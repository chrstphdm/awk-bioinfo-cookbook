#!/usr/bin/env bats
# Multi-sample recipe tests
# Run from repo root: bats docs/tests/multi_sample.bats

DATA="docs/data"

@test "Multi-sample: 5 samples in samples_list.txt" {
    run bash -c "wc -l < $DATA/samples_list.txt | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "Multi-sample: VCF missingness — NA12878 has 1 missing genotype" {
    run bash -c "awk '/^#CHROM/{for(i=10;i<=NF;i++) sn[i]=\$i; next} /^#/{next} {for(s=10;s<=NF;s++){split(\$s,g,\":\"); gt=g[1]; gsub(/\\|/,\"/\",gt); if(gt==\"./.\"||gt==\".\") m[sn[s]]++}} END{print m[\"NA12878\"]+0}' $DATA/variants.vcf"
    [ "$status" -eq 0 ]
    [ "$output" = "1" ]
}

@test "Multi-sample: VCF missingness — NA20585 has 0 missing genotypes" {
    run bash -c "awk '/^#CHROM/{for(i=10;i<=NF;i++) sn[i]=\$i; next} /^#/{next} {for(s=10;s<=NF;s++){split(\$s,g,\":\"); gt=g[1]; gsub(/\\|/,\"/\",gt); if(gt==\"./.\"||gt==\".\") m[sn[s]]++}} END{print m[\"NA20585\"]+0}' $DATA/variants.vcf"
    [ "$status" -eq 0 ]
    [ "$output" = "0" ]
}

@test "Multi-sample: metadata QC — 3 PASS, 1 FAIL, 1 WARN" {
    run bash -c "awk 'NR>1{c[\$5]++} END{print c[\"PASS\"], c[\"FAIL\"], c[\"WARN\"]}' $DATA/metadata.tsv"
    [ "$status" -eq 0 ]
    [ "$output" = "3 1 1" ]
}

@test "Multi-sample: idxstats total mapped reads" {
    run bash -c "awk '\$1!=\"*\"{m+=\$3} END{print m}' $DATA/alignments.idxstats"
    [ "$status" -eq 0 ]
    [ "$output" = "483310" ]
}
