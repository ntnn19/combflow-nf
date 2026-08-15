nextflow run main.nf -profile slurm,apptainer,conda --fasta assets/example.fasta \
    --pdb_dir /scratch/home/nagarnat/projects/scaturro/zika_ns5/replication_complex/pdb \
    --stoichiometry 'NS1:1-5,NS2A:1-5,NS2B:1-5,NS3:1-5,NS4A:1-5,NS4B:1-5,NS5:1-5' -resume
