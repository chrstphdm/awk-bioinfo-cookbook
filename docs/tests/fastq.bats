#!/usr/bin/env bats
# FASTQ recipe tests
# Run from repo root: bats docs/tests/fastq.bats

DATA="docs/data"

@test "FASTQ: 20 reads total" {
    run bash -c "awk 'NR%4==1' $DATA/reads.fastq | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "20" ]
}

@test "FASTQ: reads.fastq has exactly 80 lines (20 reads × 4)" {
    run bash -c "wc -l < $DATA/reads.fastq | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "80" ]
}

@test "FASTQ: length filter >=50 keeps 13 reads" {
    run bash -c "awk 'NR%4==2 && length(\$0)>=50' $DATA/reads.fastq | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "13" ]
}

@test "FASTQ: 2 duplicate IDs (DUP.read1 appears twice)" {
    run bash -c "awk 'NR%4==1 && /^@DUP.read1/' $DATA/reads.fastq | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "FASTQ: dedup by ID keeps first occurrence only" {
    run bash -c "awk 'NR%4==1{id=substr(\$1,2); if(id in seen) {getline;getline;getline;next} seen[id]=1} 1' $DATA/reads.fastq | awk 'NR%4==1' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    # 20 reads - 2 duplicates (one occurrence each of DUP.read1 and DUP.read2 is kept) = 18
    [ "$output" = "18" ]
}

@test "FASTQ: 5 IDs in ids.txt" {
    run bash -c "wc -l < $DATA/ids.txt | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}
