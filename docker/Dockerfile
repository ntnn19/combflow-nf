# combfold-parallel: ntnn19/combfold + GNU parallel + click
#
# The base image provides CombFold at /app/CombFold-master (compiled
# CombinatorialAssembler) with python 3.12, biopython, scipy and numpy.
# This derivative adds:
#   - GNU parallel  (fan-out of CombFold jobs with -j $(nproc))
#   - click         (prepare_subunits_json step)
#
# Build:   docker build -t ntnn19/combfold-parallel:latest docker/
# Push:    docker push ntnn19/combfold-parallel:latest
# Apptainer (on the HPC):
#   apptainer build combfold-parallel.sif docker://ntnn19/combfold-parallel:latest
FROM ntnn19/combfold:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends parallel \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir click==8.3.0
