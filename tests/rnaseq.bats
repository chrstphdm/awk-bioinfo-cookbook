#!/usr/bin/env bats
# RNA-seq recipe tests
# Run from repo root: bats docs/tests/rnaseq.bats

DATA="docs/data"

@test "HTSeq: 15 gene lines (excluding __ lines)" {
    run bash -c "awk '\$1 !~ /^__/' $DATA/htseq_counts.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "15" ]
}

@test "HTSeq: 5 summary __ lines" {
    run bash -c "awk '\$1 ~ /^__/' $DATA/htseq_counts.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "HTSeq: __no_feature count is 18432" {
    run bash -c "awk '\$1==\"__no_feature\"{print \$2}' $DATA/htseq_counts.tsv"
    [ "$status" -eq 0 ]
    [ "$output" = "18432" ]
}

@test "featureCounts: 15 gene rows (excluding headers)" {
    run bash -c "awk '/^##/{next} NR>2' $DATA/featurecounts.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "15" ]
}

@test "featureCounts: 3 sample columns" {
    run bash -c "awk '/^##/{next} NR==2{print NF-6; exit}' $DATA/featurecounts.tsv"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}

@test "featureCounts: clean matrix extraction produces header + 15 rows" {
    run bash -c "awk '/^##/{next} NR==2{printf \"gene_id\"; for(i=7;i<=NF;i++) printf \"\t\"\$i; print\"\"; next} {printf \$1; for(i=7;i<=NF;i++) printf \"\t\"\$i; print\"\"}' $DATA/featurecounts.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "16" ]
}

@test "featureCounts: filter low-count (>=10 in >=2 samples) keeps 14 genes" {
    run bash -c "awk '/^##/{next} NR==2{printf \"gene_id\"; for(i=7;i<=NF;i++) printf \"\t\"\$i; print\"\"; next} {printf \$1; for(i=7;i<=NF;i++) printf \"\t\"\$i; print\"\"}' $DATA/featurecounts.tsv | awk -v mc=10 -v ms=2 'NR==1{next} {n=0; for(i=2;i<=NF;i++) if(\$i+0>=mc) n++; if(n>=ms) k++} END{print k+0}'"
    [ "$status" -eq 0 ]
    [ "$output" = "14" ]
}
