process COLLECT_RESULTS {
    tag 'collect_results'
    publishDir params.outdir, mode: 'copy'

    input:
    path comb_dirs

    output:
    path "combfold_summary.tsv"

    script:
    """
    printf 'comb\\tmodel\\tconfidence\\n' > combfold_summary.tsv
    for d in */; do
        comb=\${d%/}
        conf="\${d}assembled_results/confidence.txt"
        if [ -f "\$conf" ]; then
            awk -v c="\$comb" '{ n=split(\$1, a, "/"); print c "\\t" a[n] "\\t" \$2 }' "\$conf" >> combfold_summary.tsv
        fi
    done
    { head -n 1 combfold_summary.tsv && tail -n +2 combfold_summary.tsv | sort -k3,3nr; } > combfold_summary.tmp
    mv combfold_summary.tmp combfold_summary.tsv
    """

    stub:
    """
    printf 'comb\\tmodel\\tconfidence\\n' > combfold_summary.tsv
    """
}
