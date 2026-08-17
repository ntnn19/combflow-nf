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
include { PREPARE_FASTAS } from './modules/prepare_fastas'
include { MAP_FASTA_TO_PDB } from './modules/map_fasta_to_pdb'
include { CREATE_PDB_DIRS } from './modules/create_pdb_dirs'
include { COMBFOLD_ASSEMBLY     } from './modules/combfold_assembly'
//include { COLLECT_RESULTS       } from './modules/collect_results'

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
    jsons_ch = jsons_ch.map { files -> tuple(files.baseName, files) }
    jsons_ch.take(1).view { "jsons_ch: $it" }

    PREPARE_FASTAS(pdb_dir_ch, jsons_ch.take(1))
    
    prepare_fastas_ch  = PREPARE_FASTAS.out
    prepare_fastas_ch.take(1).view { "prepare_fastas_ch: $it" }

    MAP_FASTA_TO_PDB(pdb_dir_ch,prepare_fastas_ch.take(1))

    map_fasta_to_pdb_ch = MAP_FASTA_TO_PDB.out    
    map_fasta_to_pdb_ch.take(1).view { "map_fasta_to_pdb_ch: $it" }

    CREATE_PDB_DIRS(pdb_dir_ch,map_fasta_to_pdb_ch.take(1))

    create_pdb_dir_ch = CREATE_PDB_DIRS.out
    create_pdb_dir_ch.take(1).view { "create_pdb_dir_ch: $it" }
//    pdb_dir_per_comb_ch  = create_pdb_dir_ch.flatten()
    pdb_dir_per_comb_ch = create_pdb_dir_ch.join(jsons_ch)
    pdb_dir_per_comb_ch.take(1).view { "pdb_dir_per_comb_ch: $it" }

    // chunk combinations; chunk id = first JSON's base name (unique per chunk)
//    chunks_ch = pdb_dir_per_comb_ch.collate(chunk_size).map { files -> tuple(files[0].baseName, files) }

//    chunks_ch.take(1).view { "chunks_ch: $it" }

//    COMBFOLD_ASSEMBLY(chunks_ch, pdb_dir_ch)

//    COLLECT_RESULTS(COMBFOLD_ASSEMBLY.out.comb_dirs.collect())
}
