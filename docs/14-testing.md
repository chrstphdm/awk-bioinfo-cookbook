# Testing AWK Scripts

AWK scripts are often written as one-liners and never tested formally. This works fine
for interactive use, but breaks silently in pipelines when input format changes, edge
cases appear, or the script is ported to a different AWK implementation.

This chapter shows how to test AWK recipes: from a simple `diff`-based check to a
full test suite with [bats](https://github.com/bats-core/bats-core), integrated into
CI.

---

## The minimal test pattern: `diff`

### Context

The simplest test: run a recipe on known input and compare the output to a known-good
expected result. A non-zero exit code means the test failed.

### Code

```bash
# Test: filter VCF variants by QUAL >= 30
# Expected: only PASS/high-QUAL variants, headers intact

awk '/^#/ || ($6 != "." && $6+0 >= 30)' docs/data/variants.vcf \
  | diff - tests/expected/vcf_qual30.txt
echo "Exit: $?"   # 0 = pass, non-zero = fail
```

```bash
# Compact pass/fail wrapper
run_test() {
    local name="$1"; shift
    if "$@" > /tmp/actual.txt 2>&1 && diff -q /tmp/actual.txt "tests/expected/${name}.txt" > /dev/null 2>&1; then
        echo "PASS: $name"
    else
        echo "FAIL: $name"
        diff "tests/expected/${name}.txt" /tmp/actual.txt
        return 1
    fi
}

run_test vcf_qual30 awk '/^#/ || ($6 != "." && $6+0 >= 30)' docs/data/variants.vcf
```

### Generate expected output files

Before writing tests, generate the expected outputs from a known-good run:

```bash
mkdir -p tests/expected

# VCF QUAL filter
awk '/^#/ || ($6 != "." && $6+0 >= 30)' docs/data/variants.vcf \
    > tests/expected/vcf_qual30.txt

# FASTQ length filter (min 50 bp)
awk 'NR%4==1{h=$0; getline s; getline p; getline q;
     if(length(s)>=50) print h"\n"s"\n+\n"q}' docs/data/reads.fastq \
    > tests/expected/fastq_min50.txt

# Join: variants with metadata
awk 'NR==FNR && NR>1{m[$1]=$4; next}
     /^#CHROM/{print $0, "population"; next}
     /^#/{print; next}
     {print $0, ($10 in m ? m[$10] : ".")}' \
    docs/data/metadata.tsv docs/data/variants.vcf \
    > tests/expected/vcf_with_population.txt
```

---

## Shell harness: multiple tests with pass/fail count

### Context

A shell script that runs a series of tests and reports a summary. Useful before
committing or deploying a script.

### Code

```bash
#!/usr/bin/env bash
# tests/run_tests.sh
# Usage: bash tests/run_tests.sh (from repo root)

PASS=0; FAIL=0
DATA="docs/data"
EXP="tests/expected"

check() {
    local name="$1"; local cmd="$2"; local expected="$3"
    actual=$(eval "$cmd" 2>&1)
    if [ "$actual" = "$(cat "$expected")" ]; then
        echo "  PASS  $name"
        ((PASS++))
    else
        echo "  FAIL  $name"
        diff <(echo "$actual") "$expected" | head -10
        ((FAIL++))
    fi
}

echo "=== AWK Cookbook Tests ==="

# --- VCF ---
check "vcf_skip_headers" \
    "awk '!/^#/' $DATA/variants.vcf | wc -l | tr -d ' '" \
    <(echo "20")

check "vcf_qual30" \
    "awk '/^#/ || (\$6 != \".\" && \$6+0 >= 30)' $DATA/variants.vcf | awk '!/^#/' | wc -l | tr -d ' '" \
    <(echo "14")

check "vcf_pass_only" \
    "awk '/^#/ || \$7==\"PASS\"' $DATA/variants.vcf | awk '!/^#/' | wc -l | tr -d ' '" \
    <(echo "12")

# --- FASTQ ---
check "fastq_read_count" \
    "awk 'NR%4==1' $DATA/reads.fastq | wc -l | tr -d ' '" \
    <(echo "20")

check "fastq_min50" \
    "awk 'NR%4==2 && length(\$0)>=50' $DATA/reads.fastq | wc -l | tr -d ' '" \
    <(echo "11")

# --- FASTA ---
check "fasta_seq_count" \
    "awk '/^>/' $DATA/genome.fasta | wc -l | tr -d ' '" \
    <(echo "5")

# --- BED ---
check "bed_chr1_features" \
    "awk '\$1==\"chr1\"' $DATA/regions.bed | wc -l | tr -d ' '" \
    <(echo "7")

# --- HTSeq ---
check "htseq_gene_count" \
    "awk '\$1 !~ /^__/' $DATA/htseq_counts.tsv | wc -l | tr -d ' '" \
    <(echo "15")

echo ""
echo "Results: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ] && exit 0 || exit 1
```

---

## bats — Bash Automated Testing System

### Context

[bats-core](https://github.com/bats-core/bats-core) (MIT licence) is a TAP-compliant
test framework for shell scripts. It provides clean `@test` syntax, `run` to capture
command output, and assertions like `assert_output`.

### Installation

```bash
# via npm (any platform)
npm install -g bats

# via Homebrew (macOS)
brew install bats-core

# via apt (Debian/Ubuntu)
sudo apt-get install bats
```

### Example test file

```bash
#!/usr/bin/env bats
# tests/vcf.bats

DATA="docs/data"

@test "VCF: correct variant count (20 non-header lines)" {
    run bash -c "awk '!/^#/' $DATA/variants.vcf | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "20" ]
}

@test "VCF: QUAL>=30 filter keeps 14 variants" {
    run bash -c "awk '/^#/ || (\$6 != \".\" && \$6+0 >= 30)' $DATA/variants.vcf | awk '!/^#/' | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "14" ]
}

@test "VCF: PASS filter keeps 12 variants" {
    run awk '/^#/ || $7=="PASS"' "$DATA/variants.vcf"
    [ "$status" -eq 0 ]
    variant_count=$(echo "$output" | awk '!/^#/' | wc -l | tr -d ' ')
    [ "$variant_count" = "12" ]
}

@test "VCF: sample names extracted from #CHROM line" {
    run bash -c "awk '/^#CHROM/{for(i=10;i<=NF;i++) print \$i}' $DATA/variants.vcf"
    [ "$status" -eq 0 ]
    [ "$(echo "$output" | wc -l | tr -d ' ')" = "3" ]
    echo "$output" | grep -q "NA12878"
    echo "$output" | grep -q "NA20585"
}
```

### Run the tests

```bash
bats tests/vcf.bats
bats tests/              # run all .bats files in the directory
```

---

## Testing edge cases

### Context

Edge cases that commonly break AWK scripts: empty files, single-line input, missing
values, and malformed records.

### Code

```bash
# Test with empty input
echo "" | awk '!/^#/{count++} END{print count+0}'   # should print 0, not error

# Test with header-only VCF (no variants)
printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" \
  | awk '!/^#/{count++} END{print count+0}'          # should print 0

# Test with a QUAL of "." (missing)
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" \
       "chr1\t100\t.\tA\tT\t.\t.\t.\n" \
  | awk '!/^#/ && $6!="." && $6+0>=30'               # should output nothing

# Test FASTQ with exactly 1 read
printf "@read1\nACGT\n+\nIIII\n" \
  | awk 'NR%4==2{count++} END{print count}'           # should print 1
```

```bash
#!/usr/bin/env bats
# tests/edge_cases.bats

@test "Empty file: VCF filter returns no output" {
    run awk '!/^#/{count++} END{print count+0}' /dev/null
    [ "$status" -eq 0 ]
    [ "$output" = "0" ]
}

@test "Missing QUAL (dot): not counted in QUAL>=30 filter" {
    input=$(printf 'chr1\t100\t.\tA\tT\t.\t.\t.\n')
    run bash -c "echo '$input' | awk '\$6!=\".\" && \$6+0>=30'"
    [ "$status" -eq 0 ]
    [ -z "$output" ]
}

@test "Single-read FASTQ: correctly counted" {
    run bash -c "printf '@r1\nACGT\n+\nIIII\n' | awk 'NR%4==2{c++} END{print c}'"
    [ "$status" -eq 0 ]
    [ "$output" = "1" ]
}
```

---

## Test POSIX compatibility across AWK implementations

### Context

Recipes marked `[gawk]` use GNU AWK extensions. Verify that POSIX recipes work with
`mawk` and BSD `awk`, and that `[gawk]` recipes correctly fail on non-gawk.

### Code

```bash
#!/usr/bin/env bash
# tests/portability.sh — run a recipe with multiple AWK implementations

recipe='!/^#/{c++} END{print c+0}'
input="docs/data/variants.vcf"
expected=20

for awk_impl in awk gawk mawk; do
    if command -v "$awk_impl" &>/dev/null; then
        result=$($awk_impl "$recipe" "$input")
        if [ "$result" = "$expected" ]; then
            echo "PASS ($awk_impl): $result"
        else
            echo "FAIL ($awk_impl): got $result, expected $expected"
        fi
    else
        echo "SKIP ($awk_impl): not installed"
    fi
done
```

```bash
# Test with Docker for isolated environments
docker run --rm -v "$(pwd):/repo" debian:bookworm-slim \
    bash -c "cd /repo && awk '!/^#/{c++}END{print c+0}' docs/data/variants.vcf"
# Should print 20 with Debian's mawk
```

---

## CI integration

Add a test job to `.github/workflows/deploy-docs.yml` that runs before the MkDocs build:

```yaml
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Install bats
        run: sudo apt-get install -y bats
      - name: Run AWK recipe tests
        run: bats tests/

  deploy:
    needs: test          # only deploy if tests pass
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: '3.x' }
      - run: pip install mkdocs-material
      - run: mkdocs gh-deploy --force
```
