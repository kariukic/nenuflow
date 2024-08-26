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
    AverageDataInTime;
    GetSolPerDir;
    DP3Calibrate;
    SubtractSources;
    ApplyDI;
    AOqualityCollect;
} from './neap_processes.nf'

params.stage = null

params.ch_in = null

params.ms = null


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

    if ( params.stage == "L3" ) {

        L3 ( params.ch_in )

    }


}


workflow L2_A {

    take:

        start_ch

    main:

        mset_ch = channel.fromPath( "${params.ms}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )
        main_field_model_ch = ModelToolBuild (mset_ch.take(1), params.catalog, params.min_flux, params.model_build_radius, "main_field_intrinsic_model.txt")

        ateams_ch = Channel.of( params.intrinsic_Ateam_model )

        model_attenuate_ch = ModelToolAttenuate ( true, mset_ch, params.minimum_elevation, params.minimum_patch_flux, params.min_flux, main_field_model_ch.concat( ateams_ch ).collect(), 'l2a_apparent.catalog' )

        catalogs_ch = mset_ch.collect { it + "/l2a_apparent.catalog" }

        make_sourcedb_ch = MakeSourceDB ( model_attenuate_ch.apparent_model.collect(), catalogs_ch.flatten())

        sourcedb_ch = mset_ch.collect { it + "/l2a_apparent.catalog.sourcedb" }

        mset_and_sourcedb_ch = mset_ch.flatten().merge( sourcedb_ch.flatten() )

        if ( params.sols_per_dir ) {

            mset_and_sourcedb_and_nsols_ch = mset_and_sourcedb_ch.merge( catalogs_ch.flatten() )

        }

        else {

            mset_and_sourcedb_and_nsols_ch = mset_and_sourcedb_ch.combine( channel.of( false ) )

        }

        di_calibration_ch = DP3Calibrate ( make_sourcedb_ch.collect(), mset_and_sourcedb_and_nsols_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_a, params.solint, params.maxforks.ateam )

        all_solutions_ch =  mset_ch.collect { it + "/${params.di_calibration_solutions_file_l2_a}" }

        sources_to_subtract_ch = mset_ch.collect { it + "/l2a_apparent.catalog.txt" }

        mset_sourcedb_solutions_and_sources_to_subtract_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).merge( sources_to_subtract_ch.flatten() )

        subtract_ateam_ch = SubtractSources ( di_calibration_ch.collect(), mset_sourcedb_solutions_and_sources_to_subtract_ch, params.ateams_subtraction_parset, "DATA", "SUBTRACTED_DATA_L2_A" ) // apply solutions after subtracting Ateams
        
        mset_sourcedb_and_solutions_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() )

        apply_di_ch = ApplyDI ( subtract_ateam_ch.collect(), mset_sourcedb_and_solutions_ch, params.di_apply_parset,  "SUBTRACTED_DATA_L2_A", "CORRECTED_DATA_L2_A" )

        AOqualityCollect( apply_di_ch, mset_ch, "CORRECTED_DATA_L2_A" )

    emit:

        AOqualityCollect.out
}


workflow L2_B {

    take:

        wsclean_ao_model

    main:

        mset_ch = channel.fromPath( "${params.ms}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        ateam_catalog_ch = ModelToolAttenuate ( wsclean_ao_model, mset_ch, params.minimum_elevation, params.minimum_patch_flux, params.min_flux, params.intrinsic_Ateam_model, 'l2b_apparent_ateam.catalog' )

        all_ateam_catalog_ch = mset_ch.collect { it + "/l2b_apparent_ateam.catalog" }

        // Convert apparent A-team model to AO format
        ateam_bbs2model_ch = BBS2Model ( ateam_catalog_ch.apparent_model.collect(), all_ateam_catalog_ch.flatten(), "l2b_apparent_ateam.ao" )

        all_ateam_bbs2model_ch = mset_ch.collect { it + "/l2b_apparent_ateam.ao" }

        // Combine apparent A-team model with the apparent WSclean model
        wsclean_ao_model_ch = channel.of( wsclean_ao_model )
        merge_ncp_and_ateam_ch  = Combine2AOModels( ateam_bbs2model_ch.collect(), all_ateam_bbs2model_ch.flatten().combine( wsclean_ao_model_ch ), "l2b_apparent_ncp_ateam.ao" )

        all_merge_ncp_and_ateam_ch  = mset_ch.collect { it + "/l2b_apparent_ncp_ateam.ao" }

        // Convert the combined model from dp3 format to sourceDB format
        dp3_format_ch = AO2DP3Model( merge_ncp_and_ateam_ch .collect(), all_merge_ncp_and_ateam_ch .flatten(), "l2b_apparent_ncp_ateam.skymodel" )

        all_dp3_format_ch = mset_ch.collect { it + "/l2b_apparent_ncp_ateam.skymodel" }

        // Convert the combined model from dp3 format to sourceDB format
        sourcedb_ch = MakeSourceDB( dp3_format_ch.collect(), all_dp3_format_ch.flatten() )

        all_sourcedb_ch = mset_ch.collect { it + "/l2b_apparent_ncp_ateam.skymodel.sourcedb" }

        mset_and_sourcedb_ch = mset_ch.flatten().merge( all_sourcedb_ch.flatten() )

        if ( params.sols_per_dir ) {

            mset_and_sourcedb_and_nsols_ch = mset_and_sourcedb_ch.merge( all_dp3_format_ch.flatten() )

        }

        else {
            
           mset_and_sourcedb_and_nsols_ch = mset_and_sourcedb_ch.combine( channel.of( false ) )

        }

        calibrate_ch = DP3Calibrate ( sourcedb_ch.collect(), mset_and_sourcedb_and_nsols_ch, params.di_cal_ateam_parset, params.di_calibration_solutions_file_l2_b, params.solint, params.maxforks.ateam )

        all_solutions_ch =  mset_ch.collect { it + "/${params.di_calibration_solutions_file_l2_b}" }

        sources_to_subtract_ch = mset_ch.collect { it + "/l2b_apparent_ateam.catalog.txt" }

        mset_sourcedb_solutions_and_sources_to_subtract_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).merge( sources_to_subtract_ch.flatten() )

        subtract_ateam_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_sources_to_subtract_ch, params.ateams_subtraction_parset, "DATA", "SUBTRACTED_DATA_L2_B" )

        mset_sourcedb_and_solutions_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() )

        apply_di_ch = ApplyDI ( subtract_ateam_ch.collect(), mset_sourcedb_and_solutions_ch, params.di_apply_parset,  "SUBTRACTED_DATA_L2_B", "CORRECTED_DATA_L2_B" )

        AOqualityCollect( apply_di_ch, mset_ch, "CORRECTED_DATA_L2_B" )

    emit:

        AOqualityCollect.out

}


workflow L2_C {

    take:

        wsclean_ao_model

    main:

        mset_ch = channel.fromPath( "${params.ms}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        three_c_catalog_ch = ModelToolAttenuate ( wsclean_ao_model, mset_ch, params.minimum_elevation, params.minimum_3C_patch_flux, params.min_flux, params.intrinsic_3C_model, 'l2c_apparent_3c.catalog' )

        all_3c_catalog_ch = mset_ch.collect { it + "/l2c_apparent_3c.catalog" }

        bbs2model_3c_ch = BBS2Model ( three_c_catalog_ch.apparent_model.collect(), all_3c_catalog_ch.flatten(), "l2c_apparent_3c.ao" )

        all_bbs2model_3c_ch = mset_ch.collect { it + "/l2c_apparent_3c.ao" }

        wsclean_ao_model_ch = channel.of( wsclean_ao_model )
        merge_ncp_and_3c_ch = Combine2AOModels ( bbs2model_3c_ch.collect(), all_bbs2model_3c_ch.flatten().combine( wsclean_ao_model_ch ), "l2c_apparent_ncp_3c.ao" )

        all_merge_ncp_and_3c_ch = mset_ch.collect { it + "/l2c_apparent_ncp_3c.ao" }

        dp3_format_ch = AO2DP3Model ( merge_ncp_and_3c_ch.collect(), all_merge_ncp_and_3c_ch.flatten(), "l2c_apparent_ncp_3c.skymodel" )

        all_dp3_format_ch = mset_ch.collect { it + "/l2c_apparent_ncp_3c.skymodel" }

        sourcedb_ch = MakeSourceDB( dp3_format_ch.collect(), all_dp3_format_ch.flatten() )

        all_sourcedb_ch = mset_ch.collect { it + "/l2c_apparent_ncp_3c.skymodel.sourcedb" }

        mset_and_sourcedb_ch = mset_ch.flatten().merge( all_sourcedb_ch.flatten() )

        mset_and_sourcedb_and_nsols_ch = mset_and_sourcedb_ch.combine( channel.of( false ) )

        calibrate_ch = DP3Calibrate ( sourcedb_ch.collect(), mset_and_sourcedb_and_nsols_ch, params.di_cal_3c_parset, params.di_calibration_solutions_file_l2_c, params.solint_3c, params.maxforks.three_c )

        all_solutions_ch =  mset_ch.collect { it + "/${params.di_calibration_solutions_file_l2_c}" }

        sources_to_subtract_ch = mset_ch.collect { it + "/l2c_apparent_3c.catalog.txt" }

        mset_sourcedb_solutions_and_sources_to_subtract_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).merge( sources_to_subtract_ch.flatten() )

        subtract_3c_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_sources_to_subtract_ch, params.three_c_subtraction_parset, "CORRECTED_DATA_L2_B", "SUBTRACTED_DATA_L2_C" )

        AOqualityCollect( subtract_3c_ch, mset_ch, "SUBTRACTED_DATA_L2_C" )

    emit:

        AOqualityCollect.out

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

        wsclean_ao_model
    
    main:

        mset_ch = channel.fromPath ( "${params.ms}_T*.MS", glob: true, checkIfExists: true, type: 'dir' )

        source_select_ch = SelectNearbySources ( wsclean_ao_model, params.fov_center.radec, params.sky_model_radius, "l3_ncp.ao" )

        cluster_ch = MakeClusters(source_select_ch, params.number_of_clusters, "l3_ncp_clusters.ao" )

        dp3_format_ch = AO2DP3Model ( true, cluster_ch, "l3_ncp_clusters.skymodel" )

        sourcedb_ch = MakeSourceDB ( true, dp3_format_ch )

        avgmsout_ch = mset_ch.collect { it.toString().replace(params.ms, "${params.ms}_L3") }

        mset_and_avgmsout_ch = mset_ch.flatten().merge( avgmsout_ch.flatten() )

        average_ch = AverageDataInTime ( sourcedb_ch, mset_and_avgmsout_ch, params.ntimesteps_to_average, "SUBTRACTED_DATA_L2_C", "DATA" )

        avgmsout_and_sourcedb_ch = avgmsout_ch.flatten().combine( sourcedb_ch )

        calibrate_ch = DP3Calibrate ( average_ch.collect(), avgmsout_and_sourcedb_ch, params.dd_cal_parset, params.dd_calibration_solutions_file_l3, params.solint_target, params.maxforks.target )

        all_solutions_ch =  avgmsout_ch.flatten().collect { it + "/${params.dd_calibration_solutions_file_l3}" }

        clusters_ch  = MakeDP3ClustersListFile(calibrate_ch.collect(), params.number_of_clusters, "clusters_list.txt")

        mset_and_sourcedb_ch = mset_ch.flatten().combine( sourcedb_ch )

        mset_sourcedb_solutions_and_clusters_ch = mset_and_sourcedb_ch.merge( all_solutions_ch.flatten() ).combine( clusters_ch )

        subtract_ncp_ch = SubtractSources ( calibrate_ch.collect(), mset_sourcedb_solutions_and_clusters_ch, params.ncp_subtraction_parset, "SUBTRACTED_DATA_L2_C", "SUBTRACTED_DATA_L3" )

        AOqualityCollect( subtract_ncp_ch, mset_ch, "SUBTRACTED_DATA_L3" )

    emit:

        AOqualityCollect.out

}