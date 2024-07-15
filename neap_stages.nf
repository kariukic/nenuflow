#!/usr/bin/env nextflow

include {
    ModelToolBuild;
    ModelToolAttenuate;
    BBS2Model;
    MakeSourceDB;
    Combine2AOModels;
    MakeClusters;
    AO2DP3Model;
    SelectNearbySources;
    AverageDataInTime
    DP3Calibrate;
    SubtractSources;
    ApplyDI;
} from './neap_processes.nf'

params.stage = null

params.ch_in = null

params.ms = null

params.l3_averaged_ms_name = null


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
//      5. A-team subtraction
//#############################################

workflow {

    if ( params.stage == "L2_A" ) {
        
        L2_A ( params.ch_in )

    }

    if ( params.stage == "L2_B" ) {

        L2_B ( params.ch_in )

    }

    if ( params.stage == "L2_C" ) {

        L2_C ( params.ch_in )

    }

    if ( params.stage == "L2_D" ) {

        L2_D ( params.ch_in )

    }

    if ( params.stage == "L3" ) {

        L3 ( params.ch_in )

    }


}


workflow L2_A {
    take:
        start_ch

    main:

        msets_channel = channel.fromPath( "${params.ms}*.MS", glob: true, checkIfExists: true, type: 'dir' ) //.collect {it}
        main_field_model_ch = ModelToolBuild (msets_channel.take(1), params.catalog, params.min_flux, params.sky_model_radius, "main_field_intrinsic_model.txt")

        ateams_ch = Channel.of( params.intrinsic_Ateam_model )

        model_attenuate_ch = ModelToolAttenuate ( true, msets_channel, params.minimum_elevation, params.minimum_patch_flux, params.min_flux, main_field_model_ch.concat( ateams_ch ).collect(), 'apparent_sky_model_l2_a.catalog' )

        catalogs_ch = msets_channel.collect { it + "/apparent_sky_model_l2_a.catalog" }

        make_sourcedb_ch = MakeSourceDB ( model_attenuate_ch.apparent_model.collect(), catalogs_ch.flatten())

        sourcedbs_ch = msets_channel.collect { it + "/apparent_sky_model_l2_a.catalog.sourcedb" }

        msets_and_sourcedbs_ch = msets_channel.flatten().merge( sourcedbs_ch.flatten() )

        di_calibration_ch = DP3Calibrate (make_sourcedb_ch.collect(), msets_and_sourcedbs_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_a )

        solutions_ch =  msets_channel.collect { it + "/${params.di_calibration_solutions_file_l2_a}" }

        sources_to_subtract_ch = msets_channel.collect { it + "/apparent_sky_model_l2_a.catalog.txt" }

        msets_sourcedbs_solutions_and_sources_to_subtract_ch = msets_and_sourcedbs_ch.merge( solutions_ch.flatten() ).merge( sources_to_subtract_ch.flatten() )

        subtract_ateams_ch = SubtractSources ( di_calibration_ch.collect(), msets_sourcedbs_solutions_and_sources_to_subtract_ch, params.ateams_subtraction_parset, "DATA", "SUBTRACTED_DATA_L2_A" ) // apply solutions after subtracting Ateams
        
        msets_sourcedbs_and_solutions_ch = msets_and_sourcedbs_ch.merge( solutions_ch.flatten() )

        ApplyDI ( subtract_ateams_ch.collect(), msets_sourcedbs_and_solutions_ch, params.di_apply_parset,  "SUBTRACTED_DATA_L2_A", "CORRECTED_DATA_L2_A" )

    emit:

        ApplyDI.out
}   



workflow L2_B { //catalog model

    take:
        wsclean_ao_model

    main:
        msets_channel = channel.fromPath( "${params.ms}*.MS", glob: true, checkIfExists: true, type: 'dir' ).collect {it}
        // Attenuate it
        ateam_model_attenuate_ch = ModelToolAttenuate ( wsclean_ao_model, msets_channel, params.minimum_elevation, params.minimum_patch_flux, params.min_flux, params.intrinsic_Ateam_model, 'ateam_sky_model_l2_b.catalog' )

        // Convert it to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_model_attenuate_ch.apparent_model, "apparent_ateam_sky_model_l2_b.ao" )

        // Combine them
        full_ao_format_model_ch = Combine2AOModels(ateam_bbs2model_ch, wsclean_ao_model, "full_sky_model_l2_b.ao") //combined_l2b_model_ch

        // Convert the combined model from dp3 format to sourceDB format
        full_dp3_format_model_ch = AO2DP3Model(full_ao_format_model_ch, "full_sky_model_l2_b.skymodel")

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB(full_dp3_format_model_ch, "full_sky_model_l2_b.sourcedb")

        // Run DI calibration
        di_calibration_l2b_ch = DP3Calibrate ( msets_channel, full_sourcedb_format_model_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_b )

        // Subtract Ateam sources
        subtract_ateams_l2b_ch = SubtractSources( msets_channel, params.ateams_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2b_ch, ateam_model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_B" )

        // Apply the Main field solutions
        ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_B", "CORRECTED_DATA_L2_B" )

    emit:

        ApplyDI.out
}

workflow L2_C {
    take:
        wsclean_l2b_ao_model

    main:
        msets_channel = channel.fromPath( "${params.ms}*.MS", glob: true, checkIfExists: true, type: 'dir' ).collect {it}

         // Attenuate it
        ateam_model_attenuate_ch = ModelToolAttenuate ( wsclean_l2b_ao_model, msets_channel, params.minimum_elevation, params.minimum_patch_flux, params.min_flux, params.intrinsic_Ateam_model, 'ateam_sky_model_l2_c.catalog' )

        // Convert it to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_model_attenuate_ch.apparent_model, "apparent_ateam_sky_model_l2_c.ao" )

        // Combine them
        full_ao_format_model_ch = Combine2AOModels(ateam_bbs2model_ch, wsclean_l2b_ao_model, "full_sky_model_l2_c.ao") //combined_l2b_model_ch

        // Convert the combined model from dp3 format to sourceDB format
        full_dp3_format_model_ch = AO2DP3Model(full_ao_format_model_ch, "full_sky_model_l2_c.skymodel")

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB(full_dp3_format_model_ch, "full_sky_model_l2_c.sourcedb")

        // Run DI calibration
        di_calibration_l2b_ch = DP3Calibrate ( msets_channel, full_sourcedb_format_model_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_c )

        // Subtract Ateam sources
        subtract_ateams_l2b_ch = SubtractSources ( msets_channel, params.ateams_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2b_ch, ateam_model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_C" )

        // Apply the Main field solutions
        ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_C", "CORRECTED_DATA_L2_C" )

        // Image and convert output model from BBS to AO format
        // WScleanImage ( apply_di_l2c_ch, "CORRECTED_DATA_L2_C", "clean_ateam_sub_l2_c" )

    emit:
        // WScleanImage.out.wsclean_ao_model
        ApplyDI.out
}


workflow L2_D {
    take:
        wsclean_l2c_ao_model
    
    main:
        msets_channel = channel.fromPath( "${params.ms}*.MS", glob: true, checkIfExists: true, type: 'dir' ).collect {it}
        // params.minimum_patch_flux = 1 here
        // Attenuate the intrinsic 3C model. the wsclean model is just a workflow connector
        attenuate_3c_model_ch = ModelToolAttenuate ( wsclean_l2c_ao_model, msets_channel, params.minimum_elevation, params.minimum_3C_patch_flux, params.min_flux, params.intrinsic_3C_model, 'apparent_3c_sky_model_l2_d.catalog' )

        // Convert it to AO format
        bbs2model_3c_ch = BBS2Model ( attenuate_3c_model_ch.apparent_model, "apparent_3c_sky_model_l2_d.ao" )

        // Combine it with the wsclean model from L2C
        merge_ncp_and_3c_sky_models_ch = Combine2AOModels ( bbs2model_3c_ch, wsclean_l2c_ao_model, "full_ncp_and_3c_sky_model_l2_d.ao" )

        // Convert the combined model from dp3 format to sourceDB format
        convert_model_to_dp3_format_ch = AO2DP3Model ( merge_ncp_and_3c_sky_models_ch, "full_ncp_and_3c_sky_model_l2_d.skymodel" )

        // Convert the combined model from dp3 format to sourceDB format
        full_sourcedb_format_model_ch = MakeSourceDB ( convert_model_to_dp3_format_ch, "full_ncp_and_3c_sky_model_l2_d.sourcedb" )

        // Run DI calibration # this parset should change the solution interval and  minimum lambda cut
        di_calibration_l2d_ch = DP3Calibrate ( msets_channel, full_sourcedb_format_model_ch, params.di_cal_3c_parset, params.di_calibration_solutions_file_l2_d )

        // Subtract 3C sources
        SubtractSources ( msets_channel, params.three_c_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2d_ch, attenuate_3c_model_ch.sources_to_subtract_file, "CORRECTED_DATA_L2_C", "SUBTRACTED_DATA_L2_D" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        WScleanImage ( subtract_3c_l2d_ch, "SUBTRACTED_DATA_L2_D", "clean_3c_sub_l2_d" )

    emit:
    
        // WScleanImage.out.wsclean_ao_model
        SubtractSources.out
}




///////////////////////////////////////////////////////////////////////////////////////////
// Data: L3
//
// Steps:
//      1. Average to 12s
//      2. calibrate averaged MS with 8mins (40 timesteps) solution intervals
//      3. Aply solutions to unaveraged data
//      4. Subtract the NCP using these solutions
///////////////////////////////////////////////////////////////////////////////////////////
// Bash Skymodel commands: (Offringa path replaced with singularity container)
// ~offringa/Software/lofartools/build/bbs2model sw03-sources.txt sw03-sources_ao.txt
// ~offringa/Software/lofartools/build/editmodel -m sw03-sources_pb_ao.txt -near 0h00m00s 90d00m00s 13.4 sw03-sources_ao.txt
// ~offringa/Software/lofartools/build/cluster sw03-sources_pb_ao.txt sw03-sources_cluster_ao.txt 7
// ~offringa/Software/lofartools/build/editmodel -dppp-model sw03-sources_cluster.txt sw03-sources_cluster_ao.txt
// params.dd_cal_parset = "/home/users/satyapan/NCP211212_analysis/example_parsets/ncp_sub/dd_cal.parset"
// params.dd_calibration_solutions_file_l3 = "dd_diagonal_clustered_l3.h5"
// params.ncp_subtraction_parset = "/home/users/chege/theleap/neap/parsets/dd_sub.parset"

workflow L3 {

    take:
        final_ncp_wsclean_ao_model
    
    main:
        msets_channel = channel.fromPath( "${params.ms}*.MS", glob: true, checkIfExists: true, type: 'dir' ).collect {it}

        // Select sources within params.sky_model_radius degrees around params.fov_center
        source_selection_ch = SelectNearbySources (final_ncp_wsclean_ao_model, params.fov_center, params.sky_model_radius, "dd_ncp_sky_model_l3.ao")

        // Make the needed number of clusters
        clustering_ch = MakeClusters(source_selection_ch, params.number_of_clusters, "dd_ncp_sky_model_clustered_l3.ao")

        // Convert the combined model from dp3 format to sourceDB format
        convert_dd_model_to_dp3_format_ch = AO2DP3Model ( clustering_ch, "dd_ncp_sky_model_clustered_l3.skymodel" )

         // Convert the combined model from dp3 format to sourceDB format
        dd_sourcedb_format_model_ch = MakeSourceDB ( convert_dd_model_to_dp3_format_ch, "dd_ncp_sky_model_clustered_l3.sourcedb" )

        // Average the data
        time_averaging_ch = AverageDataInTime ( dd_sourcedb_format_model_ch, msets_channel, params.l3_averaged_ms_name, params.ntimesteps_to_average, "SUBTRACTED_DATA_L2_D", "DATA")

        // Perform DD calibration on averaged data
        dd_calibration_l3_ch = DP3Calibrate ( time_averaging_ch, dd_sourcedb_format_model_ch, params.dd_cal_parset, params.dd_calibration_solutions_file_l3 )

        // Subtract NCP using the dd solutions
        SubtractSources ( msets_channel, params.ncp_subtraction_parset, dd_sourcedb_format_model_ch, dd_calibration_l3_ch, params.ncp_clusters, "SUBTRACTED_DATA_L2_D", "SUBTRACTED_DATA_L3" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        WScleanImage ( subtract_ncp_l3_ch, "SUBTRACTED_DATA_L3", "clean_ncp_sub_l3" )

    emit:
        // WScleanImage.out.wsclean_ao_model
        SubtractSources.out
}