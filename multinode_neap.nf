#!/usr/bin/env nextflow
// New Extension in Nançay Upgrading LOFAR (NenuFAR) Analysis Pipeline

include {
    RetrieveData;
    ConvertL1toL2;
    BBS2Model;
    readTxtIntoString;
    readTxtAndAppendString;
} from "./neap_processes.nf"

// All logs for the run will be stored here
params.logs_dir = null // "/home/users/chege/theleap/neap/test/logs"

// A list of the msfiles (consecutive timechunks) that make the full Spectral window
params.mslist = null //"/home/users/chege/theleap/neap/test/sw01_56min_ms_list.txt" 

// The name of all the MSfiles on each node. All timechunks should have the same name but located on different nodes
params.ms = null// "/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01_56MinChunk.MS"    //SW01_16min.MS"

// The path to the MSfiles on each node
params.datapath = null //"/data/users/lofareor/chege/nenufar/obs/L2_BP"

// The nodes in which the msfiles are distributed
params.hosts = null //"/home/users/chege/theleap/neap/test/pssh_hosts_list.txt"

// TODO: Maybe read the ms, datapath and hosts params from the mslist param

params.obsid=null // "20231208_NT04"


import groovy.json.JsonOutput
process GetParams {
    debug true
    publishDir params.logs_dir, mode: 'copy'

    input:
        val stage

    output:
        path "${stage}_params.json", emit: params_file
        val true, emit: params_standby

    script:
        """
        echo '${JsonOutput.prettyPrint(JsonOutput.toJson(params))}' > ${stage}_params.json
        """
}


process DistributedCalibration {
    debug true

    input:
        val ready
        val ch_in
        val entry
        val params_file

    output:
        val true

    script:
        """
        pssh -v -i -h ${params.hosts} -t 0 -x "cd ${params.datapath}; bash" nextflow run ${params.serial_neap} --stage ${entry} --ch_in ${ch_in} -params-file ${params_file} > ${params.logs_dir}/${entry}_try.log 2>&1
        """
}


workflow {
    // l1_ch = Retrieve()
    // l2_ch = L1toL2( l1_ch  ) //
    // l2a_ch = Run_L2A( true )
    // l2a_ao_model = "/home/users/satyapan/NT04231208_analysis/trial_l2a_9/l2a_ateam_sub-sources_11deg.ao"
    // l2b_ch = Run_L2B( l2a_ao_model )
    l2b_ao_model = "/home/users/satyapan/NT04231208_analysis/trial_l2b_2/l2b_ateam_sub-sources_11deg.ao"
    l2c_ch = Run_L2C( l2b_ao_model )
    // l2b_ao_model = "/home/users/satyapan/NT04231208_analysis/trial_l2b_2/l2b_ateam_sub-sources.ao"
    // l3_ch = Run_L3( l2b_ao_model )

}

// TODO: Run AOFLagger for post-calibration RFI Flagging

// For pspipe list all MS in time in one line
    // steps:
    // create mslist with the observation_id as the name of the file. observationId use date_field_spectralwindow


workflow Retrieve {

    main:
        RetrieveData( true, params.remote_host, params.obsid, params.config_file )

    emit:
        RetrieveData.out

}

workflow L1toL2 {

    take:
        l1_data_ready

    main:
        ConvertL1toL2( l1_data_ready, params.obsid, "L2_BP", params.config_file )

    emit:
        ConvertL1toL2.out
}


workflow Run_L2A {

    take:

        l2_data_available

    main:

        l2a_params_ch = GetParams( 'l2a' )

        cal_l2a_ch = DistributedCalibration ( l2a_params_ch.params_standby, l2_data_available, "L2_A", l2a_params_ch.params_file)

        mses = readTxtIntoString ( params.mslist )

        solution_files = readTxtAndAppendString(params.mslist, "/${params.di_calibration_solutions_file_l2_a}")

        sols_collect_ch = H5ParmCollect( cal_l2a_ch, solution_files, "di_l2_a_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.combined_sols, mses, "aoqstats_l2a" )

        wsclean_ch = WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_A", "l2a_ateam_sub" )

        BBS2Model ( true, wsclean_ch.model, "l2a_model.ao" )


    emit:
        model = BBS2Model.out

}


workflow Run_L2B {
    take:
        wsclean_ao_model
    
    main:

        l2b_params_ch = GetParams( 'l2b' )

        cal_l2b_ch = DistributedCalibration ( l2b_params_ch.params_standby, wsclean_ao_model, "L2_B", l2b_params_ch.params_file )

        mses = readTxtIntoString ( params.mslist )

        solution_files = readTxtAndAppendString(params.mslist, "/${params.di_calibration_solutions_file_l2_b}")

        sols_collect_ch = H5ParmCollect( cal_l2b_ch, solution_files, "di_l2_b_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.combined_sols, mses, "aoqstats_l2b" )

        wsclean_ch = WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_B", "l2b_ateam_sub" )

        BBS2Model ( true, wsclean_ch.model, "l2b_model.ao" )

    emit:
        model = BBS2Model.out

}

workflow Run_L2C {
    take:
        wsclean_ao_model

    main:

        l2c_params_ch = GetParams( 'l2c' )

        cal_l2c_ch = DistributedCalibration ( l2c_params_ch.params_standby, wsclean_ao_model, "L2_C", l2c_params_ch.params_file )

        mses = readTxtIntoString ( params.mslist )

        solution_files = readTxtAndAppendString(params.mslist, "/${params.di_calibration_solutions_file_l2_c}")

        sols_collect_ch = H5ParmCollect( cal_l2c_ch, solution_files, "di_l2_c_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.combined_sols, mses, "aoqstats_l2c" )

    //     wsclean_ch = WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L2_C", "l2c_3c_sub" )

    //     BBS2Model ( true, wsclean_ch.model, "l2c_model.ao" )

    // emit:
    //     model = BBS2Model.out

    
}


workflow Run_L3 {
    take:
        wsclean_ao_model

    main:

        l3_params_ch = GetParams( 'l3' )

        cal_l3_ch = DistributedCalibration ( l3_params_ch.params_standby, wsclean_ao_model, "L3", l3_params_ch.params_file )

        mses = readTxtIntoString ( params.mslist )

        String solution_files = file(params.mslist).readLines().collect { it.toString().strip().replace(params.ms, "${params.ms}_L3") }.join("/${params.dd_calibration_solutions_file_l3} ") + "/${params.dd_calibration_solutions_file_l3}"

        sols_collect_ch = H5ParmCollect( cal_l3_ch, solution_files, "dd_l3_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.combined_sols, mses, "aoqstats_l3" )

        wsclean_ch = WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L3", "l3_main_sub" )

    //     BBS2Model ( true, wsclean_ch.model, "l3_model.ao" )

    // emit:
    //     model = BBS2Model.out

}


// WSclean image
process WScleanImage {
    label 'sing'
    publishDir "${params.datapath}/results/images", pattern: "*.fits", mode: "move", overwrite: true
    publishDir "${params.datapath}/results/images", pattern: "*.txt", mode: "copy", overwrite: true

    input:
        val ready
        val mses
        val size // 1800
        val scale // 1amin
        val spectral_pol_fit // 2
        val data_column
        val image_name

    output:
        path "*.fits"
        path "${image_name}-sources.txt", emit: model

    shell:
        '''
        wsclean -name !{image_name} -pol I -weight briggs -0.1 -data-column !{data_column} -minuv-l 20 -maxuv-l 5000 -scale !{scale} -size !{size} !{size} -make-psf -niter 100000 -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms -multiscale -no-update-model-required -join-channels -channels-out 12 -save-source-list -fit-spectral-pol !{spectral_pol_fit} !{mses} > !{params.logs_dir}/!{image_name}_wsclean_image.log
        '''
}

process AOqualityCombine {

    input:
        val ready
        val mses
        val output_name

    output:
        val true, emit: qstats

    shell:
        '''
        mkdir -p !{params.datapath}/results/aoquality
        aoquality combine !{params.datapath}/results/aoquality/!{output_name}.qs !{mses} > aoquality_combine.log
        python !{projectDir}/templates/plot_aoqstats.py -q !{params.datapath}/results/aoquality/!{output_name}.qs -o !{params.datapath}/results/aoquality/!{output_name}.png >> !{params.logs_dir}/aoquality_combine.log
        '''
}

process H5ParmCollect {

    publishDir "${params.datapath}/results/solutions/${output_name}", pattern: "*.h5", mode: "move", overwrite: true
    // publishDir "${params.datapath}/results/solutions/${output_name}", pattern: "*.png"

    input:
        val ready
        val solution_files
        val output_name

    output:
        path "${output_name}.h5", emit: combined_sols
        // path "*.png", emit: combined_sols
        // val true, emit: combined_sols

    shell:
        '''
        H5parm_collector.py !{solution_files} -o !{output_name}.h5 > h5parm_collect.log
        # soltool plot --plot_dir $(pwd) !{output_name}.h5 >> h5parm_collect.log
        '''
        

}