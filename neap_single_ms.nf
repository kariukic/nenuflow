#!/usr/bin/env nextflow

// New Extension in Nançay Upgrading LOFAR (NenuFAR) Analysis Pipeline
// TODO: variables to change
//      1. Channel names
//      2. Model outputs names
//      3. Input the input params for modelattenuate at the workflow level maybe?
//      4. Combine all steps that take an intrinsic wsclean model, maybe + other models and attenuates them then does all the conversion steps to sourcedb format
//      5. Also make the di_calibration step fully parametrised in its own workflow? #Also the output solutions file


params.container = "/net/node100/data/users/lofareor/bharat/singularity/release_26-09-23/lofarpipe_26-09-23.sif"


// Copied from /net/node120/data/users/lofareor/nenufar/obs/L2_BP/20220508_200000_20220509_033100_NCP_COSMIC_DAWN/SW01.MS the took the first 240 timestamps for testing purposes
params.ms = "/net/node115/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01_16min.MS"
params.ms_avg = "SW01_16min_avg.MS"
params.make_apparent_model=true // run a 3 times selfcal loop to build a model from the data

//##############################################
//WORKFLOW: L2
//##############################################
// Summary: A-team Subtraction and DI Correction section
//
// Data: L2
//
// Steps:
//      1. model attenuate
//      2. Make sourcedp of attenuated model
//      3. DI calibration
//      4. DI solutions application
//      4. A-team subtraction
//#############################################

params.intrinsic_catalog_model = "/home/users/chege/theleap/neap/models/catalog.skymodel"
params.intrinsic_Ateam_model =  "/home/users/chege/theleap/neap/models/Ateam_LBA.skymodel"
params.intrinsic_3C_model = "/home/users/chege/theleap/neap/models/3c_sources.skymodel"
params.minimum_elevation = 0
params.minimum_patch_flux = 0

// params.attenuated_model = "apparent_sky.catalog"
// params.apparent_sourcedb = 'apparent_sky.sourcedb'

params.di_cal_ateam_parset = "/home/users/chege/theleap/neap/parsets/di_cal_apparent.parset"
params.di_cal_3c_parset="/home/users/chege/theleap/neap/parsets/di_cal_3c.parset"


params.di_calibration_solutions_file_l2_a = "di_diagonal_l2_a.h5"
params.di_calibration_solutions_file_l2_b = "di_diagonal_l2_b.h5"
params.di_calibration_solutions_file_l2_c = "di_diagonal_l2_c.h5"
params.di_calibration_solutions_file_l2_d = "di_diagonal_l2_d.h5"

params.ateams_subtraction_parset = "/home/users/chege/theleap/neap/parsets/di_ateams_subtraction.parset" //"/home/users/satyapan/NCP211212_analysis/example_parsets/ateam_di_catalog/di_sub_4.parset"
params.three_c_subtraction_parset = "/home/users/chege/theleap/neap/parsets/di_3c_subtraction.parset"

params.di_apply_parset = "/home/users/satyapan/NCP211212_analysis/example_parsets/ateam_di_catalog/di_apply.parset"


workflow {

    if params.make_apparent_model
    L2_A() //catalog

    L2_B(L2_A.out) //  apparent model 

    L2_C(L2_B.out) // apparent model 2

    L2_D(L2_C.out) // 3c subtracted

    L3(L2_D.out)
}

workflow L1 {
    main:
        retrieve_ch = RetrieveData(true, 'databf', params.obsid, params.config_file)
        ConvertL1toL2(retrieve_ch, "L2_BP", params.obsid, params.config_file)
    emit:
        ConvertL1toL2.out
}


workflow L2_A { //catalog model

    main:
        models_ch = Channel.fromPath( [params.intrinsic_catalog_model, params.intrinsic_Ateam_model] ).collect()

        model_attenuate_ch = ModelToolAttenuate ( true, params.ms, params.minimum_elevation, params.minimum_patch_flux, models_ch, 'apparent_sky_model_l2_a.catalog' )

        make_sourcedb_ch = MakeSourceDB ( model_attenuate_ch.apparent_model, 'apparent_sky_model_l2_a.sourcedb' )

        di_calibration_ch = DP3Calibrate ( params.ms, make_sourcedb_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_a )

        subtract_ateams_ch = SubtractSources( params.ms, params.ateams_subtraction_parset, make_sourcedb_ch, di_calibration_ch, model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_BP_A" ) // apply solutions after subtracting Ateams

        apply_di_ch = ApplyDI ( subtract_ateams_ch, make_sourcedb_ch, params.di_apply_parset,  di_calibration_ch, "SUBTRACTED_DATA_L2_BP_A", "CORRECTED_DATA_L2_BP_A" )

        // Image and convert output model from BBS to AO format
        wsclean_ch = WScleanImage ( apply_di_ch, "CORRECTED_DATA_L2_BP_A", "SW01_ateam_subtracted_l2_a" )

        // BBS2Model ( wsclean_ch.wsclean_model, "NCP_sky_model_l2_a.ao" )

    emit:

        WScleanImage.out.wsclean_ao_model
}


workflow L2_B { //catalog model

    take:
        wsclean_ao_model

    main:

        // Attenuate it
        ateam_model_attenuate_ch = ModelToolAttenuate ( wsclean_ao_model, params.ms, params.minimum_elevation, params.minimum_patch_flux, params.intrinsic_Ateam_model, 'ateam_sky_model_l2_b.catalog' )

        // Convert it to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_model_attenuate_ch.apparent_model, "apparent_ateam_sky_model_l2_b.ao" )

        // Grab both the aparent Ateam model and the apparent wscean model from L2_A
        // This method did not work
        // combined_l2b_model_ch = Channel.fromList ( [ateam_bbs2model_ch, wsclean_ao_model] ).collect() 

        // Combine them
        full_ao_format_model_ch = Combine2AOModels(ateam_bbs2model_ch, wsclean_ao_model, "full_sky_model_l2_b.ao") //combined_l2b_model_ch

        // Convert the combined model from dp3 format to sourceDB format
        full_dp3_format_model_ch = AO2DP3Model(full_ao_format_model_ch, "full_sky_model_l2_b.skymodel")

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB(full_dp3_format_model_ch, "full_sky_model_l2_b.sourcedb")

        // Run DI calibration
        di_calibration_l2b_ch = DP3Calibrate ( params.ms, full_sourcedb_format_model_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_b )

        // Subtract Ateam sources
        subtract_ateams_l2b_ch = SubtractSources( params.ms, params.ateams_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2b_ch, ateam_model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_BP_B" )

        // Apply the Main field solutions
        apply_di_l2b_ch = ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_BP_B", "CORRECTED_DATA_L2_BP_B" )

        // Image and convert output model from BBS to AO format
        WScleanImage ( apply_di_l2b_ch, "CORRECTED_DATA_L2_BP_B", "SW01_ateam_subtracted_l2_b" )

    emit:
        WScleanImage.out.wsclean_ao_model

}


workflow L2_C {
    take:
        wsclean_l2b_ao_model

    main:
        
         // Attenuate it
        ateam_model_attenuate_ch = ModelToolAttenuate ( wsclean_l2b_ao_model, params.ms, params.minimum_elevation, params.minimum_patch_flux, params.intrinsic_Ateam_model, 'ateam_sky_model_l2_c.catalog' )

        // Convert it to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_model_attenuate_ch.apparent_model, "apparent_ateam_sky_model_l2_c.ao" )

        // Combine them
        full_ao_format_model_ch = Combine2AOModels(ateam_bbs2model_ch, wsclean_l2b_ao_model, "full_sky_model_l2_c.ao") //combined_l2b_model_ch

        // Convert the combined model from dp3 format to sourceDB format
        full_dp3_format_model_ch = AO2DP3Model(full_ao_format_model_ch, "full_sky_model_l2_c.skymodel")

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB(full_dp3_format_model_ch, "full_sky_model_l2_c.sourcedb")

        // Run DI calibration
        di_calibration_l2b_ch = DP3Calibrate ( params.ms, full_sourcedb_format_model_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_c )

        // Subtract Ateam sources
        subtract_ateams_l2b_ch = SubtractSources ( params.ms, params.ateams_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2b_ch, ateam_model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_BP_C" )

        // Apply the Main field solutions
        apply_di_l2c_ch = ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_BP_C", "CORRECTED_DATA_L2_BP_C" )

        // Image and convert output model from BBS to AO format
        WScleanImage ( apply_di_l2c_ch, "CORRECTED_DATA_L2_BP_C", "SW01_ateam_subtracted_l2_c" )

    emit:
        WScleanImage.out.wsclean_ao_model
}


workflow L2_D {
    take:
        wsclean_l2c_ao_model
    
    main:
        // Attenuate the intrinsic 3C model. the wsclean model is just a workflow connector
        attenuate_3c_model_ch = ModelToolAttenuate ( wsclean_l2c_ao_model, params.ms, params.minimum_elevation, params.minimum_patch_flux, params.intrinsic_3C_model, 'apparent_3c_sky_model_l2_d.catalog' )

        // Convert it to AO format
        bbs2model_3c_ch = BBS2Model ( attenuate_3c_model_ch.apparent_model, "apparent_3c_sky_model_l2_d.ao" )

        // Combine it with the wsclean model from L2C
        merge_ncp_and_3c_sky_models_ch = Combine2AOModels ( bbs2model_3c_ch, wsclean_l2c_ao_model, "full_ncp_and_3c_sky_model_l2_d.ao" )

        // Convert the combined model from dp3 format to sourceDB format
        convert_model_to_dp3_format_ch = AO2DP3Model ( merge_ncp_and_3c_sky_models_ch, "full_ncp_and_3c_sky_model_l2_d.skymodel" )

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB ( convert_model_to_dp3_format_ch, "full_ncp_and_3c_sky_model_l2_d.sourcedb" )

        // Run DI calibration # this parset should change the solution interval and  minimum lambda cut
        di_calibration_l2d_ch = DP3Calibrate ( params.ms, full_sourcedb_format_model_ch, params.di_cal_3c_parset, params.di_calibration_solutions_file_l2_d )

        // Subtract 3C sources
        subtract_3c_l2d_ch = SubtractSources ( params.ms, params.three_c_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2d_ch, attenuate_3c_model_ch.sources_to_subtract_file, "CORRECTED_DATA_L2_BP_C", "SUBTRACTED_DATA_L2_BP_D" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        WScleanImage ( subtract_3c_l2d_ch, "SUBTRACTED_DATA_L2_BP_D", "SW01_3c_subtracted_l2_d" )

    emit:
        WScleanImage.out.wsclean_ao_model
}


// Skymodel commands:
// ~offringa/Software/lofartools/build/bbs2model sw03-sources.txt sw03-sources_ao.txt
// ~offringa/Software/lofartools/build/editmodel -m sw03-sources_pb_ao.txt -near 0h00m00s 90d00m00s 13.4 sw03-sources_ao.txt
// ~offringa/Software/lofartools/build/cluster sw03-sources_pb_ao.txt sw03-sources_cluster_ao.txt 7
// ~offringa/Software/lofartools/build/editmodel -dppp-model sw03-sources_cluster.txt sw03-sources_cluster_ao.txt

///////////////////////////////////////////////////////////////////////////////////////////
// Data: L3
//
// Steps:
//      1. Average to 12s
//      2. calibrate averaged MS with 8mins (40 timesteps) solution intervals
//      3. Aply solutions to unaveraged data
//      4. Subtract the NCP using these solutions 
///////////////////////////////////////////////////////////////////////////////////////////
params.dd_cal_parset = "/home/users/satyapan/NCP211212_analysis/example_parsets/ncp_sub/dd_cal.parset"
params.dd_calibration_solutions_file_l3 = "dd_diagonal_clustered_l3.h5"
params.ncp_subtraction_parset = "/home/users/chege/theleap/neap/parsets/dd_sub.parset"

workflow L3 {

    take:
        final_ncp_wsclean_ao_model
    
    main:

        // Select sources within 13.4 degrees of the NCP. Why 13.4?
        source_selection_ch = SelectSourcesNearNCP (final_ncp_wsclean_ao_model, "dd_ncp_sky_model_l3.ao", 13.4)

        // Make 7 clusters
        clustering_ch = MakeClusters(source_selection_ch, "dd_ncp_sky_model_clustered_l3.ao", 7)

        // Convert the combined model from dp3 format to sourceDB format
        convert_dd_model_to_dp3_format_ch = AO2DP3Model ( clustering_ch, "dd_ncp_sky_model_clustered_l3.skymodel" )

         // Convert the combined model from dp3 format to sourceDB format
        dd_sourcedb_format_model_ch = MakeSourceDB ( convert_dd_model_to_dp3_format_ch, "dd_ncp_sky_model_clustered_l3.sourcedb" )

        // Average the data
        time_averaging_ch = AverageDataInTime ( dd_sourcedb_format_model_ch, params.ms, params.ms_avg, 3, "SUBTRACTED_DATA_L2_BP_D", "DATA")

        // Perform DD calibration on averaged data
        dd_calibration_l3_ch = DP3Calibrate ( time_averaging_ch, dd_sourcedb_format_model_ch, params.dd_cal_parset, params.dd_calibration_solutions_file_l3 )

        // Subtract NCP using the dd solutions
        subtract_ncp_l3_ch = SubtractSources ( params.ms, params.ncp_subtraction_parset, dd_sourcedb_format_model_ch, dd_calibration_l3_ch, "/home/users/chege/theleap/neap/parsets/ncp_clusters.txt", "SUBTRACTED_DATA_L2_BP_D", "SUBTRACTED_DATA_L3" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        WScleanImage ( subtract_ncp_l3_ch, "SUBTRACTED_DATA_L3", "SW01_ncp_subtracted_l3" )

    emit:
        WScleanImage.out.wsclean_ao_model
}



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

// Retrieve the data from some remote cluster like DATABF to DAWN cluster
process RetrieveData {

    input:
        val ready
        val remote_host
        val obsid
        path config_file
    
    output:
        

    script:
        """
        nenudata retrieve ${remote_host} ${obsid} -c ${config_file}
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
        val level // L2_BP
        val config_file

    script:
        """
        nenudata l1_to_l2 ${level} ${obsid} -c ${config_file} --l1_level 'L1' --max_concurrent 1
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


//Combines multiple models and applys beam attenuation
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


//Use editmodel tool to combine multiple Ao format models
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

    input:
        path full_ms_path
        path sourcedb_name
        path dp3_cal_parset_file
        val output_calibration_solutions_file //.5 extension

    output:
        path "${output_calibration_solutions_file}"

    shell:
        '''
        DP3 !{dp3_cal_parset_file} msin=!{full_ms_path} cal.sourcedb=!{sourcedb_name} cal.h5parm=!{output_calibration_solutions_file} > dp3_cal.log
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


//WSclean image
process WScleanImage {
    publishDir "${launchDir}/results/images", pattern: "*.fits", mode: "move", overwrite: true

    input:
        path full_ms_path
        val data_column
        val image_name

    output:
        path "*.fits"
        path "${image_name}-sources.ao", emit: wsclean_ao_model
    // -niter 100000
    shell:
        '''
        wsclean -name !{image_name} -pol I -weight briggs -0.1 -data-column !{data_column} -minuv-l 20 -maxuv-l 2000 -scale 3amin -size 600 600 -make-psf -niter 100000 -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms -multiscale -no-update-model-required -join-channels -channels-out 12 -save-source-list -fit-spectral-pol 2 !{full_ms_path} > wsclean_image.log

        singularity exec --bind /net,/data !{params.container} bbs2model /net/$(hostname)/$(pwd)/!{image_name}-sources.txt /net/$(hostname)/$(pwd)/!{image_name}-sources.ao
        '''
}
