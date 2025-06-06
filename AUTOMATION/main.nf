#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.ref = ""
params.reads1 = ""
params.reads2 = ""
params.outprefix = "sample1"
params.outdir = "Results"

workflow {
    call_variants(
        file(params.ref),
        file(params.reads1),
        file(params.reads2),
        params.outprefix,
        params.outdir
    )

    run_analysis(params.outdir)
}

process call_variants {
    input:
    path ref
    path read1
    path read2
    val outprefix
    val outdir

    output:
    path "${outdir}/*"

    conda "${projectDir}/envs/variant_calling.yml"

    script:
    """
    mkdir -p ${outdir}
    bash ${projectDir}/Scripts/variant_caller.sh ${ref} ${read1} ${read2} ${outprefix}
    """
}

process run_analysis {
    input:
    val outdir

    conda "${projectDir}/envs/analysis_env.yml"

    script:
    """
    python3 ${projectDir}/Scripts/analyse_varcall.py ${outdir}
    """
}
