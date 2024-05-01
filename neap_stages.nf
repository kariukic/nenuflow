#!/usr/bin/env nextflow

include {
    ModelToolAttenuate;
    BBS2Model;
    MakeSourceDB;
    Combine2AOModels;
    MakeClusters;
    AO2DP3Model;
    SelectSourcesNearNCP;
    AverageDataInTime
    DP3Calibrate;
    SubtractSources;
    ApplyDI;
} from './neap_processes.nf'


params.stage = null

params.ch_in = null

// Copied from /net/node120/data/users/lofareor/nenufar/obs/L2_BP/20220508_200000_20220509_033100_NCP_COSMIC_DAWN/SW01.MS the took the first 240 timestamps for testing purposes
//"/net/node[114..115]/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01_16min.MS"
params.ms = null
params.ms_avg = "SW01_16min_avg.MS"

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


workflow L2_A { // catalog model
    take:
        start_ch

    main:

        models_ch = Channel.fromPath( [params.intrinsic_catalog_model, params.intrinsic_Ateam_model] ).collect()

        model_attenuate_ch = ModelToolAttenuate ( true, params.ms, params.minimum_elevation, params.minimum_patch_flux, models_ch, 'apparent_sky_model_l2_a.catalog' )

        make_sourcedb_ch = MakeSourceDB ( model_attenuate_ch.apparent_model, 'apparent_sky_model_l2_a.sourcedb' )

        di_calibration_ch = DP3Calibrate ( params.ms, make_sourcedb_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_a )

        subtract_ateams_ch = SubtractSources ( params.ms, params.ateams_subtraction_parset, make_sourcedb_ch, di_calibration_ch, model_attenuate_ch.sources_to_subtract_file, "DATA", "SUBTRACTED_DATA_L2_BP_A" ) // apply solutions after subtracting Ateams

        ApplyDI ( subtract_ateams_ch, make_sourcedb_ch, params.di_apply_parset,  di_calibration_ch, "SUBTRACTED_DATA_L2_BP_A", "CORRECTED_DATA_L2_BP_A" )

    emit:

        ApplyDI.out
}



workflow L2_B { //catalog model

    take:
        wsclean_ao_model

    main:

        // Attenuate it
        ateam_model_attenuate_ch = ModelToolAttenuate ( wsclean_ao_model, params.ms, params.minimum_elevation, params.minimum_patch_flux, params.intrinsic_Ateam_model, 'ateam_sky_model_l2_b.catalog' )

        // Convert it to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_model_attenuate_ch.apparent_model, "apparent_ateam_sky_model_l2_b.ao" )

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
        ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_BP_B", "CORRECTED_DATA_L2_BP_B" )

    emit:

        ApplyDI.out
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
        ApplyDI ( subtract_ateams_l2b_ch, full_sourcedb_format_model_ch, params.di_apply_parset,  di_calibration_l2b_ch, "SUBTRACTED_DATA_L2_BP_C", "CORRECTED_DATA_L2_BP_C" )

        // Image and convert output model from BBS to AO format
        // WScleanImage ( apply_di_l2c_ch, "CORRECTED_DATA_L2_BP_C", "SW01_ateam_subtracted_l2_c" )

    emit:
        // WScleanImage.out.wsclean_ao_model
        ApplyDI.out
}


workflow L2_D {
    take:
        wsclean_l2c_ao_model
    
    main:
        // params.minimum_patch_flux = 1 here
        // Attenuate the intrinsic 3C model. the wsclean model is just a workflow connector
        attenuate_3c_model_ch = ModelToolAttenuate ( wsclean_l2c_ao_model, params.ms, params.minimum_elevation, 1, params.intrinsic_3C_model, 'apparent_3c_sky_model_l2_d.catalog' )

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
        SubtractSources ( params.ms, params.three_c_subtraction_parset, full_sourcedb_format_model_ch, di_calibration_l2d_ch, attenuate_3c_model_ch.sources_to_subtract_file, "CORRECTED_DATA_L2_BP_C", "SUBTRACTED_DATA_L2_BP_D" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        // WScleanImage ( subtract_3c_l2d_ch, "SUBTRACTED_DATA_L2_BP_D", "SW01_3c_subtracted_l2_d" )

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

        // Select sources within 13.4 degrees of the NCP. Why 13.4?
        source_selection_ch = SelectSourcesNearNCP (final_ncp_wsclean_ao_model, "dd_ncp_sky_model_l3.ao", 15) // changed from 13.4 just because this is a lower frequency

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
        SubtractSources ( params.ms, params.ncp_subtraction_parset, dd_sourcedb_format_model_ch, dd_calibration_l3_ch, "/home/users/chege/theleap/neap/parsets/ncp_clusters.txt", "SUBTRACTED_DATA_L2_BP_D", "SUBTRACTED_DATA_L3" )

        // Image and convert output model from BBS to AO format #TODO: REMOVE THIS PART IF YOU WANT!
        // WScleanImage ( subtract_ncp_l3_ch, "SUBTRACTED_DATA_L3", "SW01_ncp_subtracted_l3" )

    emit:
        // WScleanImage.out.wsclean_ao_model
        SubtractSources.out
}