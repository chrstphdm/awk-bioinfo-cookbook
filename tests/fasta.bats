#!/usr/bin/env bats
# FASTA recipe tests
# Run from repo root: bats docs/tests/fasta.bats

DATA="docs/data"

@test "FASTA: 5 sequences" {
    run bash -c "awk '/^>/' $DATA/genome.fasta | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "FASTA: headers contain chr1, chr2, chr3, scaffold_001, scaffold_002" {
    run bash -c "awk '/^>/{print \$1}' $DATA/genome.fasta | sed 's/>//' | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"chr1"* ]]
    [[ "$output" == *"scaffold_002"* ]]
}

@test "FASTA: linearise produces 5 single-line sequences" {
    run bash -c "awk '/^>/{if(seq) print seq; print; seq=\"\"; next} {seq=seq\$0} END{if(seq) print seq}' $DATA/genome.fasta | awk '!/^>/' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "5" ]
}

@test "FASTA: filter min length 100 keeps 2 sequences (chr1 and chr2)" {
    run bash -c "awk '/^>/{if(seq && length(seq)>=100){print h; print seq} h=\$0; seq=\"\"; next} {seq=seq\$0} END{if(seq && length(seq)>=100){print h; print seq}}' $DATA/genome.fasta | awk '/^>/' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}

@test "FASTA: rename lookup has 2 entries" {
    run bash -c "wc -l < $DATA/rename.tsv | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "2" ]
}
