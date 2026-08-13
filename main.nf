#!/usr/bin/env nextflow
/*
 * combflow-nf: Nextflow port of the combflow Snakemake workflow.
 *
 * Stoichiometric screening with CombFold:
 *   1. PREPARE_SUBUNITS_JSON - expand a FASTA + stoichiometry spec into one
 *      CombFold subunits JSON per combination (Cartesian product).
 *   2. COMBFOLD_ASSEMBLY     - run CombFold assembly per combination. Jobs are
 *      chunked across Nextflow tasks; within each task the chunk's jobs run
 *      via GNU parallel with -j = nproc (CPUs of the node the task landed on).
 *   3. COLLECT_RESULTS       - merge all assemblies into a ranked summary TSV.
 */

nextflow.enable.dsl = 2

include { PREPARE_SUBUNITS_JSON } from './modules/prepare_subunits_json'
include { COMBFOLD_ASSEMBLY     } from './modules/combfold_assembly'
include { COLLECT_RESULTS       } from './modules/collect_results'

workflow {
    if (!params.fasta)         error "Missing required parameter: --fasta <path>"
    if (!params.pdb_dir)       error "Missing required parameter: --pdb_dir <path>"
    if (!params.stoichiometry) error "Missing required parameter: --stoichiometry <e.g. 'NS1:1-2,NS2A:1-2'>"

    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    // value channel: the PDB dir must be reusable by every chunk task
    // (a queue channel would be exhausted after the first chunk)
    if (!file(params.pdb_dir).isDirectory()) error "--pdb_dir is not a directory: ${params.pdb_dir}"
    pdb_dir_ch = channel.value(params.pdb_dir)

    PREPARE_SUBUNITS_JSON(fasta_ch, params.stoichiometry, params.extra_prepare_flags)

    chunk_size = params.chunk_size as int
    if (chunk_size < 1) error "--chunk_size must be >= 1, got ${params.chunk_size}"

    jsons_ch  = PREPARE_SUBUNITS_JSON.out.jsons.flatten()
    // chunk combinations; chunk id = first JSON's base name (unique per chunk)
    chunks_ch = jsons_ch.collate(chunk_size).map { files -> tuple(files[0].baseName, files) }

    COMBFOLD_ASSEMBLY(chunks_ch, pdb_dir_ch)

    COLLECT_RESULTS(COMBFOLD_ASSEMBLY.out.comb_dirs.collect())
}
