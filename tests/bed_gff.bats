#!/usr/bin/env bats
# BED/GFF recipe tests
# Run from repo root: bats docs/tests/bed_gff.bats

DATA="docs/data"

@test "BED: 15 features total" {
    run bash -c "awk '!/^#/' $DATA/regions.bed | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "15" ]
}

@test "BED: 7 features on chr1" {
    run bash -c "awk '\$1==\"chr1\"' $DATA/regions.bed | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "7" ]
}

@test "BED: filter by chrUn keeps 2 features" {
    run bash -c "awk '\$1~/^chrUn/' $DATA/regions.bed | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "BED: merge overlapping intervals preserves count (no overlaps in test data)" {
    run bash -c "sort -k1,1 -k2,2n $DATA/regions.bed | awk 'BEGIN{OFS=\"\t\"} FNR==1{c=\$1;s=\$2;e=\$3;next} \$1==c&&\$2<=e{if(\$3>e)e=\$3;next} {print c,s,e; c=\$1;s=\$2;e=\$3} END{if(c!=\"\")print c,s,e}' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "15" ]
}

@test "GFF3: 2 gene features" {
    run bash -c "awk '!/^#/ && \$3==\"gene\"' $DATA/annotation.gff3 | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "GFF3: 5 exon features total" {
    run bash -c "awk '!/^#/ && \$3==\"exon\"' $DATA/annotation.gff3 | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "GFF3: extract ID attribute from gene lines" {
    run bash -c "awk -F'\t' '!/^#/ && \$3==\"gene\"{n=split(\$9,a,\";\"); for(i=1;i<=n;i++){gsub(/^ /,\"\",a[i]); if(a[i]~/^ID=/){sub(/ID=/,\"\",a[i]); print a[i]}}}' $DATA/annotation.gff3"
    [ "$status" -eq 0 ]
    [[ "$output" == *"gene:GENE1"* ]]
    [[ "$output" == *"gene:GENE2"* ]]
}
