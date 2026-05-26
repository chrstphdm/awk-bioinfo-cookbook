#!/usr/bin/env bats
# Reports & aggregation recipe tests
# Run from repo root: bats docs/tests/reports.bats

DATA="docs/data"

@test "Reports: group-by count — 2 chromosomes in VCF" {
    run bash -c "awk '!/^#/{c[\$1]++} END{for(k in c) print k}' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "Reports: min/max depth on chr1" {
    run bash -c "awk '\$1==\"chr1\"{if(NR==1||\$3<mn)mn=\$3; if(\$3>mx)mx=\$3} END{print mn, mx}' $DATA/depth.tsv"
    [ "$status" -eq 0 ]
    # min should be >= 0, max should be > 0
    min=$(echo "$output" | awk '{print $1}')
    max=$(echo "$output" | awk '{print $2}')
    (( min >= 0 ))
    (( max > 0 ))
}

@test "Reports: group-by on BED — chr1 total size" {
    run bash -c "awk '!/^#/ && \$1==\"chr1\"{total+=\$3-\$2} END{print total}' $DATA/regions.bed"
    [ "$status" -eq 0 ]
    (( output > 0 ))
}

@test "Reports: pivot-like — variant count per chrom per type" {
    run bash -c "gawk '!/^#/{t=(length(\$4)==1 && length(\$5)==1)?\"SNP\":\"INDEL\"; c[\$1][t]++} END{for(k in c) print k, c[k][\"SNP\"]+0, c[k][\"INDEL\"]+0}' $DATA/variants.vcf | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"chr1"* ]]
    [[ "$output" == *"chr2"* ]]
}
