#!/usr/bin/env nextflow
params.memory = "60 GB"

process GetMSList {
    input:
        val ready
        val obsid
        val level // L1
        val spectral_window // SW03
        path config_file // data_handler.toml

    
    output:
        val true

    script:
        '''
        nenudata get_ms -c ${config_file} -s ${spectral_window} ${level} ${obsid}
        '''
}

process SelectNearbySources {
    input:
        path input_model
        val fov_center_coords
        val radius
        val output_model
    
    output:
        path "${output_model}"


    shell:
        '''
        singularity exec --bind /net,/data !{params.container} editmodel -m /net/$(hostname)/$(pwd)/!{output_model} -near !{fov_center_coords} !{radius} /net/$(hostname)/$(pwd)/!{input_model}
        '''
}


process MakeClusters {
    input:
        path input_model
        val number_of_clusters
        val output_model
    
    output:
        path "${output_model}"
    
    shell:
        '''
        singularity exec --bind /net,/data !{params.container} cluster /net/$(hostname)/$(pwd)/!{input_model} /net/$(hostname)/$(pwd)/!{output_model} !{number_of_clusters}
        '''
}


process AverageDataInTime {
    input:
        val ready
        path msin
        val msout
        val timesteps_to_average
        val input_datacolumn
        val output_datacolumn

    output:
    path "${msout}"

    shell:
        '''
        #!/bin/bash
        cat >"avg.parset" <<EOL
        msin = !{msin}
        msout = !{msout}
        msout.overwrite=True
        msin.datacolumn = "!{input_datacolumn}"
        msout.datacolumn = "!{output_datacolumn}"
        steps = [avg]
        avg.type = average
        avg.timestep = !{timesteps_to_average}
        EOL
        DP3 "avg.parset"
        '''
}


///////////////////////////////////////////////////////////////////////////////////////////
// Data: L1
//
// Steps:
//      1. Bandpass calibration
//      2. AOFlagger
//      3. Averaging
///////////////////////////////////////////////////////////////////////////////////////////

// Retrieve the data from some remote cluster like DATABF to DAWN cluster
process RetrieveData {

    input:
        val ready
        val remote_host
        val obsid
        path config_file
    
    output:
        val true

    script:
        """
        nenudata retrieve ${remote_host} ${obsid} -c ${config_file} --dry_run
        """
}

// Apply bandpass calibration
process BandpassCalibration {
    input:
        val ready
        path config_file //data_handler_toml_file
        path ms

    script:
        """
        calpipe ${config_file} ${ms}
        """
}

// Convert L1 data to L2_BP
// This steps does 3 things (based on the DP3 parset specified in the config file. Can be changed though) 
// 1. Applying Bandpass solutions obtained earlier
// 2. Flagging with AOFlagger
// 3. Averaging
process ConvertL1toL2 {

    input:
        val ready
        val obsid
        val level // L2_BP
        val config_file
    
    output:
        val true

    shell:
        """
        ulimit -n 4096
        nenudata l1_to_l2 !{level} !{obsid} -c !{config_file} --l1_level 'L1' --max_concurrent 1
        """
}


///////////////////////////////////////////////////////////////////////////////////////////
// Data: L2
//
// Steps:
//      1. Calibration in the direction of NCP + A-team sources using DDECal
//      2. Subtraction of Ateams
//      3. DI correction of NCP
///////////////////////////////////////////////////////////////////////////////////////////

process ModelToolBuild {
    errorStrategy 'ignore'
    input:
        path msin
        val catalog
        val min_flux
        val radius
        val output_catalog

    output:
        path "${output_catalog}", emit: intrinsic_model

    shell:
        '''
        modeltool build  -c !{catalog} -m !{min_flux} -r !{radius} -o !{output_catalog} !{msin} > model_build.log
        '''
}


// Combines multiple models and applys beam attenuation
process ModelToolAttenuate {
    publishDir "${full_ms_path}" , mode: 'copy'
    input:
        val ready
        path full_ms_path
        val min_elevation
        val min_patch_flux
        val min_flux
        path intrinsic_model //or models
        val output_model_name

    output:
        path "${output_model_name}", emit: apparent_model
        path "${output_model_name}.txt", emit: sources_to_subtract_file

    shell:
        '''
        modeltool attenuate !{full_ms_path} !{intrinsic_model} -e !{min_elevation} -p !{min_patch_flux} -o !{output_model_name} -m !{min_flux}  > model_attenuate.log
        python3 !{projectDir}/templates/apparent_sources_left.py -i model_attenuate.log -o !{output_model_name}.txt -e 'Main'
        '''
}


// Make sourcedb sky model format given a different model catalog format
process MakeSourceDB {
    
    input:
        val ready
        val input_model
        // val sourcedb_name

    output:
        // path("${sourcedb_name}", type: 'dir')
        val true


    shell:
        '''
        rm -rf !{input_model}.sourcedb
        makesourcedb in=!{input_model} out=!{input_model}.sourcedb > make_source_db.log
        '''
}


//Convert bbs to AO format
process BBS2Model {

    input:
        val input_model
        val output_model

    output:
        path "${output_model}"


    shell:
        '''
        singularity exec --bind /net,/data !{params.container} bbs2model /net/$(hostname)/!{input_model} /net/$(hostname)/$(pwd)/!{output_model}
        '''
}


//Use editmodel tool to combine multiple AO format models
process Combine2AOModels {

    input:
        path input_model1
        path input_model2
        val output_model

    output:
        path "${output_model}"


    shell:
        '''
        singularity exec --bind /net,/data !{params.container} editmodel -m /net/$(hostname)/$(pwd)/!{output_model} /net/$(hostname)/$(pwd)/!{input_model1} /net/$(hostname)/$(pwd)/!{input_model2}
        '''
}


//Convert AO to DP3 format
process AO2DP3Model {

    input:
        path input_model
        val output_model

    output:
        path "${output_model}"


    shell:
        '''
        singularity exec --bind /net,/data !{params.container} editmodel -dppp-model /net/$(hostname)/$(pwd)/!{output_model} /net/$(hostname)/$(pwd)/!{input_model}
        python3 !{projectDir}/templates/change_patch_name.py -i /net/$(hostname)/$(pwd)/!{output_model} -x no_patch -y Main
        '''
}


// Given a DP3 DI parset run DDECAL DI calibration
process DP3Calibrate {
    publishDir "${full_ms_path}" , mode: 'copy'
    // max memory allocation
    memory "${params.memory}"
    cpus "${params.cpus}"
    // retry 4 times upon failure
    // errorStrategy 'retry'
    // maxRetries 3

    input:
        val ready
        tuple path(full_ms_path), path(sourcedb_name)
        path dp3_cal_parset_file
        val output_calibration_solutions_file //.5 extension

    output:
        path "${output_calibration_solutions_file}"

    // # sol_file_name=!{output_calibration_solutions_file%.*}
    // timestamp=$(date +%FT%T)

    shell:

        '''
        DP3 !{dp3_cal_parset_file} msin=!{full_ms_path} cal.sourcedb=!{sourcedb_name} cal.h5parm=!{output_calibration_solutions_file} > !{full_ms_path}/!{output_calibration_solutions_file}_dp3_cal.log
        '''
}


// apply DI solutions given the solutions file and a DP3 DI slutions apply parset
process ApplyDI {

    input:
        val ready
        tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file)
        path di_apply_parset_file
        val input_datacolumn
        val output_datacolumn

    output:
        path "${full_ms_path}"

    shell:
        '''
        DP3 !{di_apply_parset_file} msin=!{full_ms_path} apply.parmdb=!{calibration_solutions_file} msin.datacolumn=!{input_datacolumn} msout.datacolumn=!{output_datacolumn} > di_apply.log
        '''
}


//Subtract a sky direction(s)
//In total we subtract A-teams, 3C sources, and NCP(in 7 clusters)
process SubtractSources {

    input:
        val ready
        tuple path(full_ms_path), path(sourcedb_name), path(calibration_solutions_file), path(sources_to_subtract_file)
        path subtraction_parset
        val input_datacolumn
        val output_datacolumn

    output:
        path "${full_ms_path}"
    
    shell:
        '''
        directions_to_subtract=$(<!{sources_to_subtract_file})
        DP3 !{subtraction_parset} msin=!{full_ms_path} sub.applycal.parmdb=!{calibration_solutions_file} sub.sourcedb=!{sourcedb_name} sub.directions=${directions_to_subtract} msin.datacolumn=!{input_datacolumn} msout.datacolumn=!{output_datacolumn}> di_sub.log
        '''
}

// Collect data quality statistics
process AOqualityCollect {
    
    input:
        val ready
        path full_ms_path
        val data_column

    output:
        path "${full_ms_path}"

    shell:
        '''
        aoquality collect -d !{data_column} !{full_ms_path}
        '''
}


def readTxtIntoString (txt) {
    List tlist = file(txt).readLines()
    String tstring = tlist.collect {"${it}"}.join(" ")

    return tstring
}

def readTxtAndAppendString (txt, str) {
    List tlist = file(txt).readLines()
    String tstring = tlist.collect {"${it}" + str}.join(" ")

    return tstring
}