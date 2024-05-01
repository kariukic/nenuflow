#!/usr/bin/env nextflow
// New Extension in Nançay Upgrading LOFAR (NenuFAR) Analysis Pipeline
// TODO: variables to change
//      1. Channel names
//      2. Model outputs names
//      3. Input the input params for modelattenuate at the workflow level maybe?
//      4. Combine all stages that take an intrinsic wsclean model, maybe + other models and attenuates them then does all the conversion stages to sourcedb format
//      5. Also make the di_calibration stage fully parametrised in its own workflow? #Also the output solutions file


//TODO: 
// collect and plot solutions
// Streamlit it
// Sort out results directory structure
// Maybe memory and time limits?
// Include the data splitting process

params.logs_dir = "/home/users/chege/theleap/neap/test/logs"

params.mslist = "/home/users/chege/theleap/neap/test/sw01_56min_ms_list.txt" // TODO: change this into  function that accepts a list of nodes and returns a list of MS by combining the datapath and the ms-name params
List mslist = file(params.mslist).readLines()
String mses = mslist.collect {"${it}"}.join(" ")

params.ms = "/data/users/lofareor/chege/nenufar/obs/L2_BP/SW01_56MinChunk.MS"    //SW01_16min.MS"
params.datapath = "/data/users/lofareor/chege/nenufar/obs/L2_BP"
params.hosts = "/home/users/chege/theleap/neap/test/pssh_hosts_list.txt"

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


workflow {

    l2a_ch = Run_L2A()
    l2b_ch = Run_L2B(l2a_ch)
    l2c_ch = Run_L2C(l2b_ch)
    l2d_ch = Run_L2D(l2c_ch)
    l3_ch = Run_L3(l2d_ch)

    // Run AOFLagger for post-calibration RFI Flagging
    // Run AOqulity collect to  get the Aoquality statistics

    // For pspipe list all MS in time in one line
        // steps:
        // create mslist with the observation_id as the name of the file. observationId use date_field_spectralwindow
        // 

}


workflow Run_L2A {

    main:
        cal_l2a_ch = DistributedCalibration ( true, true, params.ms, params.datapath, params.serial_neap, "L2_A", params.hosts )

        WScleanImage ( cal_l2a_ch, mses, "CORRECTED_DATA_L2_BP_A", "SW01_ateam_subtracted_l2_a" )

    emit:
        model = WScleanImage.out.wsclean_ao_model
}


workflow Run_L2B {
    take:
        wsclean_ao_model
    
    main:
        cal_l2b_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_B", params.hosts )

        WScleanImage ( cal_l2b_ch, mses, "CORRECTED_DATA_L2_BP_B", "SW01_ateam_subtracted_l2_b" )

    emit:
        model = WScleanImage.out.wsclean_ao_model

}


workflow Run_L2C {

    take:
        wsclean_ao_model

    main:
        cal_l2c_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_C", params.hosts )

        WScleanImage ( cal_l2c_ch, mses, "CORRECTED_DATA_L2_BP_C", "SW01_ateam_subtracted_l2_c" )

    emit:

        model = WScleanImage.out.wsclean_ao_model

}


workflow Run_L2D {
    take:
        wsclean_ao_model
    

    main:
        cal_l2d_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L2_D", params.hosts )

        WScleanImage ( cal_l2d_ch, mses, "SUBTRACTED_DATA_L2_BP_D", "SW01_3c_subtracted_l2_d" )

    emit:
        model = WScleanImage.out.wsclean_ao_model
    
}


workflow Run_L3 {
    take:
        wsclean_ao_model

    main:
        cal_l3_ch = DistributedCalibration ( true, wsclean_ao_model, params.ms, params.datapath, params.serial_neap, "L3", params.hosts )

        WScleanImage ( cal_l3_ch, mses, "SUBTRACTED_DATA_L3", "SW01_ncp_subtracted_l3" )

    emit:
        model = WScleanImage.out.wsclean_ao_model

}



//WSclean image
process WScleanImage {
    publishDir "${params.datapath}/results/images", pattern: "*.fits", mode: "move", overwrite: true

    input:
        val ready
        val mses
        val data_column
        val image_name

    output:
        path "*.fits"
        path "${image_name}-sources.ao", emit: wsclean_ao_model

    // # make image size and the pixel size here a parameter, fit-spectral-pol
    shell:
        '''
        # wsclean -name !{image_name} -pol I -weight briggs -0.1 -data-column !{data_column} -minuv-l 20 -maxuv-l 2000 -scale 3amin -size 600 600 -make-psf -niter 100000 -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms -multiscale -no-update-model-required -join-channels -channels-out 12 -save-source-list -fit-spectral-pol 2 !{mses} > wsclean_image.log

        wsclean -name !{image_name} -pol I -weight briggs -0.1 -data-column !{data_column} -minuv-l 20 -maxuv-l 5000 -scale 1amin -size 1800 1800 -make-psf -niter 100000 -auto-mask 3 -auto-threshold 1 -mgain 0.6 -local-rms -multiscale -no-update-model-required -join-channels -channels-out 12 -save-source-list -fit-spectral-pol 2 !{mses} > wsclean_image.log

        singularity exec --bind /net,/data !{params.container} bbs2model $(pwd)/!{image_name}-sources.txt $(pwd)/!{image_name}-sources.ao
        '''
}