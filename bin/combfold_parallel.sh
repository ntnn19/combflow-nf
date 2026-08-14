#!/usr/bin/env bash
# combfold_parallel.sh <pdb_dir> <outdir> <combfold_home> <json...>
#
# Runs one CombFold assembly per subunits JSON via GNU parallel, with
# -j = nproc (the number of CPUs of the node this task landed on).
#
# The patched run_on_pdbs.py resolves "../CombinatorialAssembler" and
# "scripts/libs" relative to its own location, so a symlink shadow-tree
# mirroring the CombFold layout is created in the current directory:
#   shadow/scripts/run_on_pdbs.py      (copy of the patched script)
#   shadow/scripts/libs            ->  <combfold_home>/scripts/libs
#   shadow/CombinatorialAssembler  ->  <combfold_home>/CombinatorialAssembler
# This works with Docker and with read-only Apptainer container filesystems.
set -uo pipefail

if [ "$#" -lt 4 ]; then
    echo "usage: combfold_parallel.sh <pdb_dir> <outdir> <combfold_home> <json...>" >&2
    exit 2
fi

PDB_DIR=$1
OUTDIR=$2
COMBFOLD_HOME=$3
shift 3
JSONS=("$@")

# --- shadow tree -----------------------------------------------------------
mkdir -p shadow/scripts
cp "$(command -v run_on_pdbs.py)" shadow/scripts/run_on_pdbs.py
ln -sfn "${COMBFOLD_HOME}/scripts/libs" shadow/scripts/libs
ln -sfn "${COMBFOLD_HOME}/CombinatorialAssembler" shadow/CombinatorialAssembler

#mkdir -p "$OUTDIR"

J=$(nproc)
echo "combfold_parallel: running ${#JSONS[@]} CombFold job(s) via GNU parallel -j ${J}"

run_one() {
    json=$1
    comb=$(basename "$json" .json)
    out="${OUTDIR}/${comb}"
    rm -rf "$out"
#    mkdir -p "$out"
    python shadow/scripts/run_on_pdbs.py "$json" "$PDB_DIR" "$out"
}
export -f run_one
export OUTDIR PDB_DIR

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
