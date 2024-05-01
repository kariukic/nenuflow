#!/usr/bin/env nextflow
params.memory = "60 GB"

process SelectSourcesNearNCP {
    input:
        path input_model
        val output_model
        val radius
    
    output:
        path "${output_model}"


    shell:
        '''
        singularity exec --bind /net,/data !{params.container} editmodel -m /net/$(hostname)/$(pwd)/!{output_model} -near 0h00m00s 90d00m00s !{radius} /net/$(hostname)/$(pwd)/!{input_model}
        '''
}


process MakeClusters {
    input:
        path input_model
        val output_model
        val number_of_clusters
    
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


// Apply bandpass calibration
process BandpassCalibration {
    input:
        val ready
        path toml_file
        path ms

    script:
        """
        calpipe ${tomlfile} ${ms}
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
    pass
}


// Combines multiple models and applys beam attenuation
process ModelToolAttenuate {

    input:
        val ready
        val full_ms_path
        val min_elevation
        val min_patch_flux
        path intrinsic_model //or models
        val output_model_name

    output:
        path "${output_model_name}", emit: apparent_model
        path "${output_model_name}.txt", emit: sources_to_subtract_file

    shell: //TODO: parametrize the -m option
        '''
        modeltool attenuate !{full_ms_path} !{intrinsic_model} -e !{min_elevation} -p !{min_patch_flux} -o !{output_model_name} -m 0.01 > model_attenuate.log
        python3 !{projectDir}/templates/apparent_sources_left.py -i model_attenuate.log -o !{output_model_name}.txt -e 'Main'
        '''
}


// Make sourcedb sky model format given a different model catalog format
process MakeSourceDB {

    input:
        path input_model
        val sourcedb_name

    output:
        path("${sourcedb_name}", type: 'dir')


    shell:
        '''
        makesourcedb in=!{input_model} out=!{sourcedb_name} > make_source_db.log
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
    // retry 4 times upon failure
    errorStrategy 'retry'
    maxRetries 3

    input:
        path full_ms_path
        path sourcedb_name
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
        path full_ms_path
        path sourcedb_name //apparent_node$j.sourcedb
        path di_apply_parset_file
        path calibration_solutions_file //.5 extension // full path
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
        path full_ms_path
        path subtraction_parset
        path sourcedb_name
        path calibration_solutions_file //.5 extension // full path
        val sources_to_subtract_file
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
