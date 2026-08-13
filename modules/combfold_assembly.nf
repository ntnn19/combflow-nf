process COMBFOLD_ASSEMBLY {
    tag "chunk_${chunk_id}"
    publishDir "${params.outdir}/combfold", mode: 'copy', pattern: 'parallel_chunk_*.log'
    publishDir "${params.outdir}/combfold", mode: 'copy', saveAs: { fn -> fn - ~/^out\// }

    container params.combfold_image
    cpus params.combfold_cpus
    memory params.combfold_memory
    time params.combfold_time

    input:
    tuple val(chunk_id), path(jsons, stageAs: 'jsons/*')
    path pdb_dir

    output:
    path "out/*", emit: comb_dirs
    path "parallel_chunk_${chunk_id}.log", emit: joblog

    script:
    """
    combfold_parallel.sh "${pdb_dir}" out "${params.combfold_home}" jsons/*.json
    mv parallel.log "parallel_chunk_${chunk_id}.log"
    """

    stub:
    """
    mkdir -p out/stub_comb/assembled_results
    echo "/x/out/stub_comb/assembled_results/output_clustered_0.cif 0.9" > out/stub_comb/assembled_results/confidence.txt
    echo stub > out/stub_comb/combfold.log
    echo stub > "parallel_chunk_${chunk_id}.log"
    """
}
