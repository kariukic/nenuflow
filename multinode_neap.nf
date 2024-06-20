#!/usr/bin/env nextflow
// New Extension in Nançay Upgrading LOFAR (NenuFAR) Analysis Pipeline

include {
    RetrieveData;
    ConvertL1toL2;
    readTxtIntoString;
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


process DistributedCalibration {
    debug true
    // errorStrategy 'ignore'

    input:
    val ready
    val ch_in
    val ms
    val datapath
    val serial_neap
    val entry
    path pssh_hosts

    output:
    val true

    script:
    """
    pssh -v -i -h ${pssh_hosts} -t 0 -x "cd ${datapath}; bash" nextflow run ${serial_neap} --ms ${ms} --stage ${entry} --ch_in ${ch_in} > ${params.logs_dir}/${entry}_try.log 2>&1
    """
}

//H5parm_collector.py
//soltool plot
workflow {

    // l1_ch = Retrieve()
    // l2_ch = L1toL2( true ) //l1_ch 
    l2a_ch = Run_L2A( true) // l2_ch
    l2b_ch = Run_L2B( l2a_ch )
    l2c_ch = Run_L2C( l2b_ch )
    l2d_ch = Run_L2D( l2c_ch )
    l3_ch = Run_L3( l2d_ch )

}

// Run AOFLagger for post-calibration RFI Flagging
// Run AOquality collect to  get the Aoquality statistics

// For pspipe list all MS in time in one line
    // steps:
    // create mslist with the observation_id as the name of the file. observationId use date_field_spectralwindow
    // 


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
        cal_l2a_ch = DistributedCalibration ( true, l2_data_available, params.ms, params.datapath, params.serial_neap, "L2_A", params.hosts )

        // this string is used by wsclean 
        mses = readTxtIntoString (params.mslist)

        WScleanImage ( cal_l2a_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_BP_A", "SW01_ateam_subtracted_l2_a" )

    emit:
        model = WScleanImage.out.wsclean_ao_model
}


workflow Run_L2B {
    take:
        wsclean_ao_model
    
    main:
        cal_l2b_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_B", params.hosts )

        // this string is used by wsclean 
        mses = readTxtIntoString (params.mslist)

        WScleanImage ( cal_l2b_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_BP_B", "SW01_ateam_subtracted_l2_b" )

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

        WScleanImage ( cal_l2c_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "CORRECTED_DATA_L2_BP_C", "SW01_ateam_subtracted_l2_c" )

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

        WScleanImage ( cal_l2d_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L2_BP_D", "SW01_3c_subtracted_l2_d" )

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

        WScleanImage ( cal_l3_ch, mses, params.image_size, params.image_scale, params.spectral_pol_fit, "SUBTRACTED_DATA_L3", "SW01_ncp_subtracted_l3" )

    emit:
        model = WScleanImage.out.wsclean_ao_model

}

//WSclean image
process WScleanImage {
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

        singularity exec --bind /net,/data !{params.container} bbs2model $(pwd)/!{image_name}-sources.txt $(pwd)/!{image_name}-sources.ao
        '''
}