nextflow run main.nf -profile local_apptainer,local_conda --fasta assets/example.fasta \
    --pdb_dir /scratch/home/nagarnat/projects/scaturro/zika_ns5/replication_complex/pdb --stoichiometry 'NS1:1-2,NS2A:1-2'
