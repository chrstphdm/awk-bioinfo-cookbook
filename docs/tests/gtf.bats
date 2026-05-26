#!/usr/bin/env bats
# GTF recipe tests
# Run from repo root: bats docs/tests/gtf.bats

DATA="docs/data"

@test "GTF: 2 gene records" {
    run bash -c "awk -F'\t' '!/^#/ && \$3==\"gene\"' $DATA/annotation.gtf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "GTF: extract gene_id with gawk match()" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"gene\"{delete a; match(\$9,/gene_id \"([^\"]+)\"/,a); print a[1]}' $DATA/annotation.gtf | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"GENE1"* ]]
    [[ "$output" == *"GENE2"* ]]
}

@test "GTF: extract gene_name" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"gene\"{delete a; match(\$9,/gene_name \"([^\"]+)\"/,a); print a[1]}' $DATA/annotation.gtf | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"GeneName1"* ]]
    [[ "$output" == *"GeneName2"* ]]
}

@test "GTF: GENE1.1 has 3 exons" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"exon\"{match(\$9,/transcript_id \"([^\"]+)\"/,a); c[a[1]]++} END{print c[\"GENE1.1\"]}' $DATA/annotation.gtf"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}

@test "GTF: GENE2.1 has 2 exons" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"exon\"{match(\$9,/transcript_id \"([^\"]+)\"/,a); c[a[1]]++} END{print c[\"GENE2.1\"]}' $DATA/annotation.gtf"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "GTF: GENE1.1 spliced length = 602 bp" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"exon\"{delete a; match(\$9,/transcript_id \"([^\"]+)\"/,a); l[a[1]]+=\$5-\$4+1} END{print l[\"GENE1.1\"]}' $DATA/annotation.gtf"
    [ "$status" -eq 0 ]
    [ "$output" = "602" ]
}

@test "GTF: detect format as GENCODE (chr prefix)" {
    run bash -c "awk -F'\t' '!/^#/{print (\$1~/^chr/ ? \"GENCODE\" : \"Ensembl\"); exit}' $DATA/annotation.gtf"
    [ "$status" -eq 0 ]
    [ "$output" = "GENCODE" ]
}
