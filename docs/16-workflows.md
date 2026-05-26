# End-to-end Workflows

These workflows chain multiple AWK recipes into complete pipelines. Each one
demonstrates how the individual format chapters combine in practice.

!!! tip "Test data"
    All workflows use the files in [`docs/data/`](data/README.md). Run from the repo root.

---

## Workflow 1: FASTQ QC → filter → stats report

**Problem:** You have raw FASTQ reads. You want to know how many pass a length
threshold, and generate a before/after report.

**Chapters used:** [FASTQ](02-fastq.md)

```bash
#!/usr/bin/env bash
# fastq_qc_pipeline.sh
# Usage: bash fastq_qc_pipeline.sh reads.fastq 50

FASTQ="${1:-docs/data/reads.fastq}"
MIN_LEN="${2:-50}"

echo "=== FASTQ QC Report ==="
echo "Input file : $FASTQ"
echo "Min length : $MIN_LEN"
echo ""

# Step 1: raw stats
echo "--- Raw reads ---"
awk 'NR%4==2 {
    len = length($0)
    total++
    bases += len
    if (len < min_len) short++
}
END {
    printf "Reads      : %d\n", total
    printf "Bases      : %d\n", bases
    printf "Mean len   : %.1f bp\n", bases/total
    printf "Short (<min_len) : %d\n", short+0
}' OFS='\t' \
   "min_len=$MIN_LEN" \
   "$FASTQ" 2>/dev/null || \
awk -v min_len="$MIN_LEN" 'NR%4==2 {
    len = length($0); total++; bases += len
    if (len < min_len) short++
}
END {
    printf "Reads      : %d\n", total
    printf "Bases      : %d\n", bases
    printf "Mean len   : %.1f bp\n", (total>0?bases/total:0)
    printf "Short (<"min_len"bp) : %d\n", short+0
}' "$FASTQ"

echo ""

# Step 2: filter reads by minimum length, write to temp file
FILTERED=$(mktemp /tmp/filtered_XXXXXX.fastq)
awk -v min_len="$MIN_LEN" '
NR%4==1 { h=$0; getline s; getline p; getline q
    if (length(s) >= min_len)
        printf "%s\n%s\n+\n%s\n", h, s, q
}' "$FASTQ" > "$FILTERED"

# Step 3: post-filter stats and retention rate
echo "--- Filtered reads (min ${MIN_LEN} bp) ---"
awk -v min_len="$MIN_LEN" '
NR%4==2 { total++; bases += length($0) }
END {
    printf "Reads      : %d\n", total
    printf "Bases      : %d\n", bases
    printf "Mean len   : %.1f bp\n", (total>0?bases/total:0)
}' "$FILTERED"

# Compute retention rate
raw_count=$(awk 'NR%4==1' "$FASTQ" | wc -l | tr -d ' ')
filt_count=$(awk 'NR%4==1' "$FILTERED" | wc -l | tr -d ' ')
awk -v raw="$raw_count" -v kept="$filt_count" \
    'BEGIN { printf "\nRetention  : %d/%d (%.1f%%)\n", kept, raw, kept/raw*100 }'

rm -f "$FILTERED"
```

**Expected output on `docs/data/reads.fastq` with min_len=50:**
The exact numbers depend on the test data; run the pipeline to see your results.
Key outputs: read/base counts before and after filtering, and the retention rate.

---

## Workflow 2: VCF + BED → variants annotated with region

**Problem:** Given a VCF and a BED file of regions of interest, flag which variants
fall inside a region.

**Chapters used:** [VCF](05-vcf.md), [BED/GFF](04-bed-gff.md), [Two-file Joins](07-joins.md)

```awk
# annotate_vcf_with_bed.awk
# Usage: awk -f annotate_vcf_with_bed.awk regions.bed variants.vcf

# Step 1: load BED regions into memory (first file: NR==FNR)
NR == FNR {
    if (/^#/ || NF < 3) next
    chrom = $1
    n[chrom]++
    # BED is 0-based half-open; convert to 1-based closed to match VCF POS
    r_start[chrom][n[chrom]] = $2 + 1
    r_end[chrom][n[chrom]]   = $3
    r_name[chrom][n[chrom]]  = (NF >= 4) ? $4 : "region_" n[chrom]
    next
}

# Step 2: process VCF (second file)
/^#CHROM/ { print $0 "\tREGION"; next }
/^#/      { print; next }
{
    chrom = $1; pos = $2 + 0
    region = "."
    for (r = 1; r <= n[chrom]+0; r++) {
        if (pos >= r_start[chrom][r] && pos <= r_end[chrom][r]) {
            region = r_name[chrom][r]
            break
        }
    }
    print $0 "\t" region
}
```

```bash
# Run the workflow
awk -f annotate_vcf_with_bed.awk docs/data/regions.bed docs/data/variants.vcf \
  | awk '!/^#/ && $NF != "." { print $1, $2, $4, $5, $NF }' \
  | column -t
```

!!! warning "Performance note"
    This recipe does a linear scan of all regions per variant: O(V × R).
    For a whole-genome VCF (millions of variants) and a large region file (thousands of
    intervals), this is too slow. Use `bedtools intersect` for production:
    ```bash
    bedtools intersect -a variants.vcf -b regions.bed -wa -wb | awk '...'
    ```

---

## Workflow 3: featureCounts → clean matrix → CPM → top expressed genes

**Problem:** featureCounts output needs: (1) metadata columns stripped, (2) low-count
genes removed, (3) CPM normalisation, (4) top 20 expressed genes reported.

**Chapters used:** [RNA-seq Counts](12-rnaseq.md)

```bash
#!/usr/bin/env bash
# rnaseq_pipeline.sh
COUNTS="${1:-docs/data/featurecounts.tsv}"
MIN_COUNT=10
MIN_SAMPLES=2
TOP_N=5   # use 5 for test data (only 15 genes total)

echo "=== RNA-seq Count Pipeline ==="

# Step 1: strip featureCounts metadata columns
CLEAN=$(awk '/^##/{next}
             NR==2{ printf "gene_id"; for(i=7;i<=NF;i++) printf "\t"$i; print""; next}
             {printf $1; for(i=7;i<=NF;i++) printf "\t"$i; print""}' "$COUNTS")

n_genes=$(echo "$CLEAN" | tail -n +2 | wc -l | tr -d ' ')
echo "Genes before filter: $n_genes"

# Step 2: filter low-count genes
FILTERED=$(echo "$CLEAN" | awk -v mc="$MIN_COUNT" -v ms="$MIN_SAMPLES" '
    NR==1{print; next}
    {n=0; for(i=2;i<=NF;i++) if($i+0>=mc) n++; if(n>=ms) print}')

n_kept=$(echo "$FILTERED" | tail -n +2 | wc -l | tr -d ' ')
echo "Genes after filter (≥${MIN_COUNT} counts in ≥${MIN_SAMPLES} samples): $n_kept"

# Step 3: compute CPM
CPM=$(echo "$FILTERED" | awk '
NR==1{print; next}
{for(i=2;i<=NF;i++) s[i]+=$i; r[NR]=$0; for(i=2;i<=NF;i++) v[NR][i]=$i; g[NR]=$1}
END{
    for(r=2;r<=NR;r++){
        printf "%s", g[r]
        for(i=2;i<=NF;i++) printf "\t%.2f", v[r][i]/s[i]*1e6
        print""
    }
}')

# Step 4: top N expressed genes by mean CPM
echo ""
echo "--- Top $TOP_N genes by mean CPM ---"
echo "$CPM" | awk -v top="$TOP_N" '
NR==1{n=NF-1; next}
{s=0; for(i=2;i<=NF;i++) s+=$i; print $1, s/n}' \
  | sort -k2,2rn | head -"$TOP_N" \
  | awk 'BEGIN{printf "%-12s %12s\n","gene_id","mean_CPM"}
         {printf "%-12s %12.1f\n", $1, $2}'
```

---

## Workflow 4: multi-sample VCF → cohort summary report

**Problem:** From a multi-sample VCF, produce a per-sample summary: total variants,
SNP count, INDEL count, PASS rate, and missingness rate.

**Chapters used:** [VCF](05-vcf.md), [Multi-sample Patterns](13-multi-sample.md)

```awk
# cohort_report.awk
# Usage: awk -f cohort_report.awk variants.vcf

/^#CHROM/ {
    for (i = 10; i <= NF; i++) samples[i] = $i
    n_samples = NF - 9
    next
}
/^#/ { next }

{
    n_variants++

    # Classify variant type
    is_indel = (length($4) != 1 || (split($5, alts, ",") && length(alts[1]) != 1))

    # PASS flag
    is_pass = ($7 == "PASS")

    # Find GT index in FORMAT
    n_fmt = split($9, fmt, ":")
    gt_idx = 0
    for (i = 1; i <= n_fmt; i++) if (fmt[i] == "GT") { gt_idx = i; break }

    for (s = 10; s <= NF; s++) {
        split($s, gdata, ":")
        gt = gt_idx ? gdata[gt_idx] : "./."
        gsub(/\|/, "/", gt)

        if (gt == "./." || gt == ".") {
            missing[s]++
        } else {
            total[s]++
            if (is_indel) indels[s]++
            else          snps[s]++
            if (is_pass)  pass_count[s]++
        }
    }
}

END {
    printf "%-15s %8s %6s %6s %7s %9s\n",
           "sample", "n_called", "SNPs", "INDELs", "PASS%", "missing%"
    for (s = 10; s <= 9 + n_samples; s++) {
        t    = total[s] + 0
        miss = missing[s] + 0
        printf "%-15s %8d %6d %6d %6.1f%% %8.1f%%\n",
               samples[s],
               t,
               snps[s]+0,
               indels[s]+0,
               (t > 0 ? pass_count[s]/t*100 : 0),
               (n_variants > 0 ? miss/n_variants*100 : 0)
    }
}
```

```bash
# Run and format
awk -f cohort_report.awk docs/data/variants.vcf
```

**Example output** (on `docs/data/variants.vcf`):
```
sample           n_called   SNPs INDELs   PASS%     miss%
NA12878               19     16      3   68.4%      5.0%
NA19238               19     16      3   68.4%      5.0%
NA20585               19     16      3   68.4%      0.0%
```
