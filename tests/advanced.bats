#!/usr/bin/env bats
# Advanced patterns + getline recipe tests
# Run from repo root: bats docs/tests/advanced.bats

DATA="docs/data"

@test "Advanced: multi-dim array — per-chrom coverage aggregation" {
    run bash -c "gawk '{cov[\$1]+=\$3; n[\$1]++} END{for(c in cov) printf \"%s %.0f\n\",c,cov[c]/n[c]}' $DATA/depth.tsv | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"chr1"* ]]
    [[ "$output" == *"chr2"* ]]
}

@test "Advanced: match() capture — extract gene_id from GFF3" {
    run bash -c "LC_ALL=C gawk -F'\t' '!/^#/ && \$3==\"gene\"{delete a; match(\$9,/ID=([^;]+)/,a); print a[1]}' $DATA/annotation.gff3 | sort"
    [ "$status" -eq 0 ]
    [[ "$output" == *"gene:GENE1"* ]]
    [[ "$output" == *"gene:GENE2"* ]]
}

@test "Advanced: gensub — remove chr prefix" {
    run bash -c "echo 'chr1 100 200' | gawk '{c=gensub(/^chr/,\"\",1,\$1); print c, \$2, \$3}'"
    [ "$status" -eq 0 ]
    [ "$output" = "1 100 200" ]
}

@test "Advanced: force reassembly \$1=\$1 — CSV to TSV" {
    run bash -c "echo 'a,b,c' | awk 'BEGIN{FS=\",\";OFS=\"\t\"}{\$1=\$1;print}'"
    [ "$status" -eq 0 ]
    [[ "$output" == *"$(printf 'a\tb\tc')"* ]]
}

@test "Advanced: split() lookup table — filter classic HLA genes" {
    run bash -c "printf 'HLA-A\nHLA-B\nHLA-Z\n' | awk 'BEGIN{split(\"A,B,C,DRB1\",t,\",\"); for(i in t) ok[\"HLA-\"t[i]]=1} \$1 in ok'"
    [ "$status" -eq 0 ]
    [[ "$output" == *"HLA-A"* ]]
    [[ "$output" == *"HLA-B"* ]]
    [[ "$output" != *"HLA-Z"* ]]
}

@test "getline: read single value from file" {
    run bash -c "echo 42 > /tmp/test_val.txt && gawk 'BEGIN{getline v < \"/tmp/test_val.txt\"; close(\"/tmp/test_val.txt\"); print v}' && rm /tmp/test_val.txt"
    [ "$status" -eq 0 ]
    [ "$output" = "42" ]
}

@test "getline: read from command" {
    run bash -c "gawk 'BEGIN{\"echo hello\" | getline r; close(\"echo hello\"); print r}'"
    [ "$status" -eq 0 ]
    [ "$output" = "hello" ]
}

@test "getline: FASTQ block read (same result as NR%4 method)" {
    run bash -c "gawk -v ml=50 'NR%4==1{h=\$0;getline s;getline p;getline q; if(length(s)>=ml) print h}' $DATA/reads.fastq | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "13" ]
}
