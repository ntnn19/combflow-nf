process PREPARE_SUBUNITS_JSON {
    tag 'prepare_subunits_json'
    publishDir params.outdir, mode: 'copy'

    conda "conda-forge::click=8.3.0"
//    container params.combfold_image

    input:
    path fasta
    val stoichiometry
    val extra_flags

    output:
    path "subunits/*.json", emit: jsons
    path "subunits/manifest.tsv", emit: manifest

    script:
    """
    prepare_subunits_json.py \\
        -f ${fasta} \\
        -o subunits \\
        -s '${stoichiometry}' \\
        ${extra_flags}
    """

    stub:
    """
    mkdir -p subunits
    echo '{"A0": {"name": "A0", "chain_names": ["A"], "start_res": 1, "sequence": "STUB"}}' > subunits/subunits_000_stub.json
    echo '{"A0": {"name": "A0", "chain_names": ["A", "B"], "start_res": 1, "sequence": "STUB"}}' > subunits/subunits_001_stub.json
    printf 'index\\tfilename\\n' > subunits/manifest.tsv
    """
}
