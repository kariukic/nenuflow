#!/usr/bin/env nextflow
// New Extension in Nançay Upgrading LOFAR (NenuFAR) Analysis Pipeline

include {
    RetrieveData;
    ConvertL1toL2;
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
//TEST DATA
// Copied from /net/node120/data/users/lofareor/nenufar/obs/L2_BP/20220508_200000_20220509_033100_NCP_COSMIC_DAWN/SW01.MS the took the first 240 timestamps for testing purposes
//"/net/node[114..115]/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01_16min.MS"

params.obsid=null // "20231208_NT04"


import groovy.json.JsonOutput
process GetParams {
    debug true
    publishDir params.logs_dir, mode: 'copy'

    input:
        val stage

    output:
        path "${stage}_params.json", emit: params_file
        val true, emit: params_json_written

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
    l2a_ch = Run_L2A( true ) //
    // l2b_ch = Run_L2B( '/home/users/chege/NT04/work/2c/b8cfd1d3ee32b16c68d8256b6e362a/ateam_subtracted_l2_a-sources.ao' ) //l2a_ch
    // l2c_ch = Run_L2C( l2b_ch )
    // l2d_ch = Run_L2D( l2c_ch )
    // l3_ch = Run_L3( l2d_ch )

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

        cal_l2a_ch = DistributedCalibration ( true, l2_data_available, "L2_A", l2a_params_ch.params_file)

        mses = readTxtIntoString ( params.mslist )
        solution_files = readTxtAndAppendString(params.mslist, "/${params.di_calibration_solutions_file_l2_a}")

        sols_collect_ch = H5ParmCollect( cal_l2a_ch, solution_files, "di_l2_a_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.comb_sols, mses, "aoqstats_l2a" )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_A", "ateam_subtracted_l2_a" )

    emit:
        model = WScleanImage.out.wsclean_ao_model
}


workflow Run_L2B {
    take:
        wsclean_ao_model
    
    main:
        l2b_params_ch = GetParams( 'l2b' )

        cal_l2b_ch = DistributedCalibration ( true, wsclean_ao_model, "L2_B", l2b_params_ch.params_file )

        mses = readTxtIntoString ( params.mslist )
        solution_files = readTxtAndAppendString(params.mslist, "/${params.di_calibration_solutions_file_l2_b}")

        sols_collect_ch = H5ParmCollect( cal_l2b_ch, solution_files, "di_l2_b_combined_solutions")

        aoq_comb_ch = AOqualityCombine( sols_collect_ch.comb_sols, mses, "aoqstats_l2b" )

        WScleanImage ( aoq_comb_ch.qstats.collect(), mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_B", "ateam_subtracted_l2_b" )

    emit:
        model = WScleanImage.out.wsclean_ao_model

}


workflow Run_L2C {

    take:
        wsclean_ao_model

    main:
        cal_l2c_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_C", params.hosts )

        // this string is used by wsclean 
        mses = readTxtIntoString (params.mslist)

        WScleanImage ( cal_l2c_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_C", "ateam_subtracted_l2_c" )

    emit:

        model = WScleanImage.out.wsclean_ao_model

}


workflow Run_L2D {
    take:
        wsclean_ao_model
    

    main:
        cal_l2d_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_D", params.hosts )

        // this string is used by wsclean 
        mses = readTxtIntoString (params.mslist)

        WScleanImage ( cal_l2d_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L2_D", "3c_subtracted_l2_d" )

    emit:
        model = WScleanImage.out.wsclean_ao_model
    
}


workflow Run_L3 {
    take:
        wsclean_ao_model

    main:
        cal_l3_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L3", params.hosts )

        // this string is used by wsclean 
        mses = readTxtIntoString (params.mslist)

        WScleanImage ( cal_l3_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L3", "ncp_subtracted_l3" )

    emit:
        model = WScleanImage.out.wsclean_ao_model

}

// WSclean image
process WScleanImage {
    label 'sing'
    publishDir "${params.datapath}/results/images", pattern: "*.fits", mode: "move", overwrite: true

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
        path "${image_name}-sources.ao", emit: wsclean_ao_model

    shell:
        '''
        wsclean -name !{image_name} -pol I -weight briggs -0.1 -data-column !{data_column} -minuv-l 20 -maxuv-l 5000 -scale !{scale} -size !{size} !{size} -make-psf -niter 100000 -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms -multiscale -no-update-model-required -join-channels -channels-out 12 -save-source-list -fit-spectral-pol !{spectral_pol_fit} !{mses} > wsclean_image.log

        bbs2model !{image_name}-sources.txt !{image_name}-sources.ao
        '''
}

process AOqualityCombine {
    publishDir "${params.datapath}/results/aoquality", mode: "copy"

    input:
        val ready
        val mses
        val output_name

    output:
        path "${output_name}.qs", type: 'dir', emit: qstats
        path "${output_name}.png"

    shell:
        '''
        aoquality combine !{output_name}.qs !{mses} > aoquality_combine.log
        python3 !{projectDir}/templates/plot_aoqstats.py -q !{output_name}.qs -o !{output_name}.png >> aoquality_combine.log
        '''
}

process H5ParmCollect {
    publishDir "${params.datapath}/results/solutions/${output_name}"
    publishDir "${params.datapath}/results/solutions/${output_name}", pattern: "*.png"

    input:
        val ready
        val solution_files
        val output_name

    output:
        path "${output_name}.h5", emit: comb_sols
        path "*.png"

    shell:
        '''
        H5parm_collector.py !{solution_files} -o !{output_name}.h5 > h5parm_collect.log
        soltool plot --plot_dir $(pwd) !{output_name}.h5 >> h5parm_collect.log
        '''
}