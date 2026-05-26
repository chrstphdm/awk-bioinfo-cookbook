#!/usr/bin/env bats
# Join recipe tests
# Run from repo root: bats docs/tests/joins.bats

DATA="docs/data"

@test "Join: metadata has 5 data rows (excluding header)" {
    run bash -c "awk 'NR>1' $DATA/metadata.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "Join: inner join VCF samples with metadata returns 3 matches" {
    # NA12878, NA19238, NA20585 are in metadata; NA12890, NA18507 are not in VCF
    run bash -c "awk 'NR==FNR && NR>1{m[\$1]=1;next} /^#CHROM/{for(i=10;i<=NF;i++) if(\$i in m) c++} END{print c+0}' $DATA/metadata.tsv $DATA/variants.vcf"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}

@test "Join: left-join keeps all 5 metadata rows" {
    # Extract VCF sample names from #CHROM, then match against all metadata rows
    run bash -c "awk 'NR==FNR && /^#CHROM/{for(i=10;i<=NF;i++) vcf[\$i]=1; next} NR==FNR{next} FNR>1{print \$1, (\$1 in vcf ? \"FOUND\" : \"MISSING\")}' $DATA/variants.vcf $DATA/metadata.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "Join: NA12890 is missing from VCF data" {
    run bash -c "awk 'NR==FNR && NR>1{vcf[\$1]=1;next} NR>1 && !(\$1 in vcf){print \$1}' $DATA/variants.vcf $DATA/metadata.tsv"
    [ "$status" -eq 0 ]
    echo "$output" | grep -q "NA12890"
}

@test "Join: PASS samples in metadata" {
    run bash -c "awk 'NR>1 && \$5==\"PASS\"' $DATA/metadata.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}
