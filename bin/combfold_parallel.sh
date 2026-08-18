#!/usr/bin/env bash
# combfold_parallel.sh <outdir> <combfold_home> <json...>
#
# Runs one CombFold assembly per subunits JSON via GNU parallel, with
# -j = nproc (the number of CPUs of the node this task landed on).
#
# Each JSON jsons/<comb>.json is assembled with its OWN PDB directory,
# staged by Nextflow next to the JSONs as pdbs/<comb>_pdbs/ (fallback:
# the single directory matching pdbs/<comb>*).
#
# The patched run_on_pdbs.py resolves "../CombinatorialAssembler" and
# "scripts/libs" relative to its own location, so a symlink shadow-tree
# mirroring the CombFold layout is created in the current directory:
#   shadow/scripts/run_on_pdbs.py      (copy of the patched script)
#   shadow/scripts/libs            ->  <combfold_home>/scripts/libs
#   shadow/CombinatorialAssembler  ->  <combfold_home>/CombinatorialAssembler
# This works with Docker and with read-only Apptainer container filesystems.
set -uo pipefail

if [ "$#" -lt 3 ]; then
    echo "usage: combfold_parallel.sh <outdir> <combfold_home> <json...>" >&2
    exit 2
fi

OUTDIR=$1
COMBFOLD_HOME=$2
shift 2
JSONS=("$@")

# --- shadow tree -----------------------------------------------------------
mkdir -p shadow/scripts
cp "$(command -v run_on_pdbs.py)" shadow/scripts/run_on_pdbs.py
ln -sfn "${COMBFOLD_HOME}/scripts/libs" shadow/scripts/libs
ln -sfn "${COMBFOLD_HOME}/CombinatorialAssembler" shadow/CombinatorialAssembler


J=$(nproc)
echo "combfold_parallel: running ${#JSONS[@]} CombFold job(s) via GNU parallel -j ${J}"

# --- per-combination PDB dir resolution ------------------------------------
# Convention: pdbs/<comb>_pdbs (as emitted by CREATE_PDB_DIRS). Fallback:
# exactly one directory matching pdbs/<comb>*. Prints the dir on stdout.
pdb_dir_for() {
    comb=$1
    if [ -d "pdbs/${comb}_pdbs" ]; then
        echo "pdbs/${comb}_pdbs"
        return 0
    fi
    matches=()
    for d in pdbs/"${comb}"*/; do
        [ -d "$d" ] && matches+=("$d")
    done
    if [ "${#matches[@]}" -eq 1 ]; then
        echo "${matches[0]%/}"
        return 0
    fi
    echo "combfold_parallel: ERROR - expected exactly one PDB dir matching pdbs/${comb}*, found ${#matches[@]}" >&2
    return 1
}

run_one() {
    json=$1
    comb=$(basename "$json" .json)
    out="${OUTDIR}/${comb}"
    pdb_dir=$(pdb_dir_for "$comb") || return 1
    rm -rf "$out"
    # run_on_pdbs.py refuses to run into a non-empty output dir, so the log
    # must NOT be created inside $out before the run - write it outside and
    # move it in afterwards.
    python shadow/scripts/run_on_pdbs.py "$json" "$pdb_dir" "$out"
    rc=$?
    return $rc
}
export -f run_one pdb_dir_for
export OUTDIR

printf '%s\n' "${JSONS[@]}" | parallel -j "$J" --joblog parallel.log run_one {}
PARALLEL_RC=$?
echo "combfold_parallel: GNU parallel exit status: ${PARALLEL_RC}"

# --- verify outputs ----------------------------------------------------------
# run_on_pdbs.py can exit 0 without producing an assembly ("Could not
# assemble"), so success is gated on the presence of confidence.txt.
FAILED=()
for json in "${JSONS[@]}"; do
    comb=$(basename "$json" .json)
    if [ ! -s "${OUTDIR}/${comb}/assembled_results/confidence.txt" ]; then
        FAILED+=("$comb")
    fi
done

if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "combfold_parallel: ERROR - ${#FAILED[@]} of ${#JSONS[@]} combination(s) produced no assembled_results/confidence.txt:" >&2
    printf '  - %s\n' "${FAILED[@]}" >&2
    echo "combfold_parallel: see ${OUTDIR}/<comb>/combfold.log for details" >&2
    exit 1
fi

echo "combfold_parallel: all ${#JSONS[@]} combination(s) assembled successfully"
