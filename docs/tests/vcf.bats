#!/usr/bin/env bats
# VCF recipe tests
# Run from repo root: bats docs/tests/vcf.bats

DATA="docs/data"

@test "VCF: 20 non-header lines" {
    run bash -c "awk '!/^#/' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "20" ]
}

@test "VCF: QUAL>=30 filter keeps 14 variants" {
    run bash -c "awk '/^#/ || (\$6 != \".\" && \$6+0 >= 30)' $DATA/variants.vcf | awk '!/^#/' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "14" ]
}

@test "VCF: PASS filter keeps 13 variants" {
    run bash -c "awk '!/^#/ && \$7==\"PASS\"' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "13" ]
}

@test "VCF: 3 sample names in #CHROM header" {
    run bash -c "awk '/^#CHROM/{for(i=10;i<=NF;i++) print \$i}' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}

@test "VCF: NA12878 is a sample name" {
    run bash -c "awk '/^#CHROM/{for(i=10;i<=NF;i++) print \$i}' $DATA/variants.vcf"
    [ "$status" -eq 0 ]
    echo "$output" | grep -q "NA12878"
}

@test "VCF: 2 multi-allelic sites (ALT contains comma)" {
    run bash -c "awk '!/^#/ && \$5~/,/' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "VCF: QUAL dot handling (missing QUAL not counted as >=30)" {
    run bash -c "printf 'chr1\t100\t.\tA\tT\t.\tPASS\t.\n' | awk '\$6!=\".\" && \$6+0>=30'"
    [ "$status" -eq 0 ]
    [ -z "$output" ]
}
