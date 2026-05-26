#!/usr/bin/env bats
# SAM-derived TSV recipe tests
# Run from repo root: bats docs/tests/sam_tsv.bats

DATA="docs/data"

@test "idxstats: 6 lines (5 chroms + unmapped)" {
    run bash -c "wc -l < $DATA/alignments.idxstats | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "6" ]
}

@test "idxstats: unmapped line has * as refname" {
    run bash -c "awk '\$1==\"*\"' $DATA/alignments.idxstats | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "1" ]
}

@test "depth: 200 positions" {
    run bash -c "wc -l < $DATA/depth.tsv | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "200" ]
}

@test "depth: mean coverage chr1 is approximately 35" {
    run bash -c "awk '\$1==\"chr1\"{s+=\$3; n++} END{printf \"%d\", s/n}' $DATA/depth.tsv"
    [ "$status" -eq 0 ]
    # Mean should be around 35 (generated with gauss(35,8))
    (( output >= 30 && output <= 40 ))
}

@test "depth: chr2 has 3 zero-depth positions (523-525)" {
    run bash -c "awk '\$1==\"chr2\" && \$3==0' $DATA/depth.tsv | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "3" ]
}
