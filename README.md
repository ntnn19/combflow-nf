# combflow-nf

Nextflow (DSL2) port of the [combflow](combflow-main) Snakemake workflow:
stoichiometric screening of protein complexes with
[CombFold](https://github.com/dina-lab3D/CombFold).

## What it does

1. **PREPARE_SUBUNITS_JSON** — expands a FASTA of unique chain sequences and a
   per-chain copy-number range (e.g. `NS1:1-2,NS2A:1-2`) into one CombFold
   `subunits.json` per stoichiometry combination (Cartesian product), plus a
   `manifest.tsv`.
2. **COMBFOLD_ASSEMBLY** — runs the CombFold combinatorial assembly
   (`run_on_pdbs.py`) for every combination against a directory of
   AlphaFold(-Multimer) PDB predictions. Combinations are split into chunks of
   `params.chunk_size` per Nextflow task; **within each task the jobs run via
   GNU parallel with `-j $(nproc)`**, i.e. using all CPUs of the node the task
   landed on.
3. **COLLECT_RESULTS** — merges all `assembled_results/confidence.txt` files
   into `combfold_summary.tsv` (`comb`, `model`, `confidence`), ranked by
   confidence (descending).

## Requirements

- Nextflow (>= 24.04)
- The CombFold container with GNU parallel and click added. The base image
  `ntnn19/combfold` contains neither, so build the derivative first:

  ```bash
  docker build -t ntnn19/combfold-parallel:latest docker/
  docker push ntnn19/combfold-parallel:latest        # or your own registry
  ```

  On the HPC, Apptainer pulls it automatically via
  `docker://ntnn19/combfold-parallel:latest`. To use a different image name:
  `--combfold_image <image>`.

## Usage

```bash
# SLURM + Apptainer (HPC)
nextflow run main.nf -profile slurm_apptainer \
    --fasta assets/example.fasta \
    --pdb_dir /path/to/af_predictions_pdb \
    --stoichiometry 'NS1:1-2,NS2A:1-2'

# single machine, Docker
nextflow run main.nf -profile local_docker \
    --fasta assets/example.fasta --pdb_dir pdb --stoichiometry 'NS1:1-2,NS2A:1-2'
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--fasta` | — (required) | FASTA of unique chain sequences |
| `--pdb_dir` | — (required) | directory with AlphaFold PDB predictions |
| `--stoichiometry` | — (required) | per-chain copy ranges, e.g. `NS1:1-2,NS2A:1-2` |
| `--outdir` | `results` | output directory |
| `--extra_prepare_flags` | `""` | extra flags for `prepare_subunits_json.py` (e.g. `"--max-af-size 1800"`) |
| `--chunk_size` | `25` | combinations per COMBFOLD_ASSEMBLY task |
| `--combfold_cpus` | `16` | CPUs per task (SLURM `--cpus-per-task`; GNU parallel uses `nproc`) |
| `--combfold_memory` | `32.GB` | memory per assembly task |
| `--combfold_time` | `24.h` | walltime per assembly task |
| `--combfold_image` | `ntnn19/combfold-parallel:latest` | container image |
| `--combfold_home` | `/app/CombFold-master` | CombFold install path inside the image |

### Profiles

| Profile | Executor | Software |
|---|---|---|
| `slurm_apptainer` | SLURM | Apptainer (autoMounts) |
| `local_apptainer` | local | Apptainer |
| `local_docker` | local | Docker |
| `local_conda` | local | conda — **prepare step only**; assembly needs the container |

## Outputs

```
results/
├── subunits/                     # one subunits_*.json per combination + manifest.tsv
├── combfold/
│   ├── <comb>/                   # per-combination CombFold output
│   │   └── assembled_results/    # assembled models (.cif) + confidence.txt
│   ├── logs/<comb>.log           # per-combination CombFold log
│   └── parallel_chunk_*.log      # GNU parallel job logs
└── combfold_summary.tsv          # all assemblies ranked by confidence
```

## Notes

- The patched `bin/run_on_pdbs.py` must live inside the CombFold tree to find
  `CombinatorialAssembler` and `scripts/libs`. Each assembly task builds a
  symlink shadow-tree in its work dir, so this also works with read-only
  Apptainer container filesystems.
- A combination that yields no `assembled_results/confidence.txt` (e.g.
  "Could not assemble") fails the task; the failed combinations are listed in
  the task's `.command.log`.
- `nproc` inside the container reflects the SLURM CPU allocation
  (cpuset/affinity), so GNU parallel's width always matches what the node
  gave the job.
