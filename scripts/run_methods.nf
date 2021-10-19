// Which methods to run
params.methods = ["music", "rctd", "spotlight"]
methods_ch = Channel.value(params.methods)

process runMusic {
    container 'csangara/spade_music:latest'

    input:
    val methods from methods_ch

    output:
    file proportions

    when:
    methods.contains("music")

    script:
    """
    Rscript ../../deconvolution/music/script_nf.R
    """

}

process runSpotlight {
    container 'csangara/spade_spotlight:latest'

    input:
    val methods from methods_ch

    output:
    file proportions

    when:
    methods.contains("spotlight")

    script:
    """
    Rscript ../../deconvolution/spotlight/script_nf.R
    """

}

process runRCTD {
    container 'csangara/spade_rctd:latest'

    input:
    val methods from methods_ch

    output:
    file proportions

    when:
    methods.contains("rctd")

    script:
    """
    Rscript ../../deconvolution/rctd/script_nf.R
    """

}

