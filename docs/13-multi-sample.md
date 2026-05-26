# Multi-sample Patterns

All previous chapters process one sample at a time. Real bioinformatics involves
cohorts: aggregating a metric across 50 BAM files, comparing variant counts between
populations, or flagging outlier samples before analysis.

!!! tip "Test data"
    - [`docs/data/metadata.tsv`](data/README.md) — 5 samples with batch/sex/population/QC_status
    - [`docs/data/alignments.idxstats`](data/README.md) — per-chromosome mapping stats
    - [`docs/data/variants.vcf`](data/README.md) — multi-sample VCF (3 samples)

---

## Aggregate a metric across all samples in a directory

### Context

You have one output file per sample — `NA12878.idxstats`, `NA19238.idxstats`, etc.
The shell loop iterates over files; AWK processes each one and emits a summary row.

### Code

```bash
# Compute mapped-read count per sample from samtools idxstats files
for f in *.idxstats; do
    sample=$(basename "$f" .idxstats)
    awk -v sample="$sample" \
        '$1 != "*" { mapped += $3 }
         END       { print sample, mapped }' "$f"
done | sort -k2,2rn > mapping_summary.tsv
```

```bash
# Same, but also compute unmapped rate
for f in *.idxstats; do
    sample=$(basename "$f" .idxstats)
    awk -v sample="$sample" '
        $1 != "*" { mapped   += $3 }
        $1 == "*" { unmapped  = $4 }
        END {
            total = mapped + unmapped
            printf "%s\t%d\t%d\t%.2f%%\n",
                   sample, mapped, unmapped,
                   (total > 0 ? unmapped/total*100 : 0)
        }' "$f"
done
```

### Explanation

- `basename "$f" .idxstats` strips the directory path and the `.idxstats` extension,
  leaving just the sample name. This keeps AWK code clean.
- The shell loop + single-file AWK pattern is the POSIX-compatible multi-sample
  approach. It requires no gawk extensions and works with any AWK.
- Pipe to `sort` after the loop to order by mapped count, or add a header row first
  with `echo` before the loop.

---

## Detect outlier samples by z-score

### Context

After collecting a per-sample metric (e.g. total mapped reads, variant count, mean
coverage), flag samples that are more than N standard deviations from the cohort mean.
Requires storing all values before computing statistics.

### Code

```awk
# Usage: awk -v threshold=3 -f flag_outliers.awk summary.tsv
# Input: 2-column TSV, no header: sample_id  metric_value

{ val[$1] = $2+0; sum += $2+0; count++ }

END {
    if (count == 0) exit
    mean = sum / count

    # Compute population standard deviation
    ss = 0
    for (s in val) ss += (val[s] - mean)^2
    sd = sqrt(ss / count)

    printf "%-15s %10s %8s %6s\n", "sample", "value", "z_score", "flag"
    for (s in val) {
        z = (sd > 0) ? (val[s] - mean) / sd : 0
        flag = (z > threshold || z < -threshold) ? "OUTLIER" : "OK"
        printf "%-15s %10.1f %8.2f %6s\n", s, val[s], z, flag
    }
}
```

```bash
# Example: flag samples with unusual total mapped reads
for f in *.idxstats; do
    s=$(basename "$f" .idxstats)
    awk -v s="$s" '$1!="*"{m+=$3} END{print s, m}' "$f"
done | awk -v threshold=2 -f flag_outliers.awk
```

### Explanation

- Population SD (dividing by `count`) is used here rather than sample SD (`count-1`)
  because the goal is to describe the observed cohort, not estimate a population parameter.
- `(sd > 0) ? ... : 0` prevents division by zero when all samples have identical values.
- For large cohorts, pipe the loop output to the AWK script rather than collecting into
  a file first.

---

## Compare variant counts between cohorts

### Context

Join a per-sample variant count table with a metadata table (population, batch) to
compare between groups.

### Code

```awk
# File 1: metadata.tsv (sample_id batch sex population QC_status)
# File 2: variant_counts.tsv (sample_id n_variants)
# Output: sample_id n_variants batch population QC_status

NR == FNR {
    if (NR == 1) next                    # skip header
    meta[$1] = $2 "\t" $4 "\t" $5       # batch, population, QC_status
    next
}
FNR == 1 { print $0, "batch", "population", "QC_status"; next }
{
    info = ($1 in meta) ? meta[$1] : "unknown\tunknown\tunknown"
    print $0, info
}
```

```bash
# Generate variant counts per sample from multi-sample VCF
awk '/^#CHROM/{for(i=10;i<=NF;i++) sn[i]=$i; next}
     /^#/{next}
     {for(i=10;i<=NF;i++){ split($i,g,":"); if(g[1]!="./." && g[1]!="0/0") c[sn[i]]++ }}
     END{for(s in c) print s, c[s]}' variants.vcf > variant_counts.tsv

# Then join with metadata
awk -f join_meta.awk metadata.tsv variant_counts.tsv
```

---

## Summarise QC flags across a cohort

### Context

After running per-sample QC, you have a file with one row per sample and a QC status
column. Summarise the distribution of statuses.

### Code

```awk
# Input: any TSV with QC_status in a known column
# Usage: awk -v qc_col=5 -f qc_summary.awk metadata.tsv

NR == 1 { next }    # skip header
{ total++; flag_count[$qc_col]++ }

END {
    printf "QC Summary — %d samples\n", total
    printf "%-15s %6s %8s\n", "Status", "Count", "Percent"
    for (flag in flag_count)
        printf "%-15s %6d %7.1f%%\n",
               flag, flag_count[flag], flag_count[flag]/total*100
}
```

```bash
# Direct one-liner on metadata.tsv (QC_status is column 5)
awk -v qc_col=5 'NR>1{total++;c[$qc_col]++}
                 END{for(f in c) printf "%-10s %d (%.1f%%)\n",f,c[f],c[f]/total*100}' \
    metadata.tsv
```

---

## Per-sample missingness rate from a multi-sample VCF

### Context

High genotype missingness in a sample indicates low coverage or a failed library.
Count `./. ` genotypes per sample as a fraction of total variants.

### Code

```awk
# [works with POSIX awk]
/^#CHROM/ {
    for (i = 10; i <= NF; i++) samples[i] = $i
    n_samples = NF - 9
    next
}
/^#/ { next }
{
    n_variants++
    for (s = 10; s <= NF; s++) {
        split($s, gdata, ":")
        gt = gdata[1]
        gsub(/\|/, "/", gt)                       # normalise phased
        if (gt == "./." || gt == ".") missing[s]++
    }
}
END {
    printf "%-15s %10s %12s\n", "sample", "n_missing", "miss_rate"
    for (s = 10; s <= 9 + n_samples; s++) {
        miss = missing[s] + 0
        printf "%-15s %10d %11.4f\n",
               samples[s], miss, miss / n_variants
    }
}
```

```bash
# Usage
awk -f missingness.awk variants.vcf | sort -k3,3rn
```

### Explanation

- `missing[s]` uses the column index as key, matched to `samples[s]` captured from the
  `#CHROM` line.
- Phased genotypes (`0|1`) are normalised to `0/1` before checking for `.` alleles.
- For multi-allelic missing (`.|.`), the `gsub` converts to `./.` which is then matched.

---

## Shell loop + AWK reference

### Context

The canonical pattern for processing a cohort: the shell handles the loop and filename
logic; AWK handles the data processing for each file.

### Code

```bash
# Pattern 1: one AWK call per file, aggregate in a second AWK pass
for f in data/*.idxstats; do
    sample=$(basename "$f" .idxstats)
    awk -v s="$sample" '{ ... }' "$f"
done | awk '{ aggregate }' > cohort_summary.tsv

# Pattern 2: pass all files to a single AWK invocation (gawk FILENAME/FNR)
gawk 'FNR==1 { file_idx++; sample=FILENAME; sub(/.*\//,"",sample); sub(/\.[^.]+$/,"",sample) }
      { process each line using file_idx or sample }' data/*.idxstats

# Pattern 3: use GNU parallel for large cohorts
parallel 'awk -v s={/.} -f process.awk {}' ::: data/*.idxstats | awk '{ merge }' > output.tsv
```

### When to use each pattern

| Pattern | Use when |
|---|---|
| Shell loop + pipe to AWK | Files have the same format; final aggregation is simple |
| Single AWK with multiple files | `FILENAME`/`FNR` logic is clean; all files fit comfortably in one AWK invocation |
| GNU Parallel | Cohort is large (>100 files) and each file takes non-trivial time to process |

!!! tip "Handling gzipped files"
    AWK cannot read `.gz` files directly. Use process substitution or `zcat`:

    ```bash
    # Process substitution (bash)
    awk '...' <(zcat sample.fastq.gz)

    # Or pipe via zcat in the loop
    for f in data/*.fastq.gz; do
        sample=$(basename "$f" .fastq.gz)
        zcat "$f" | awk -v s="$sample" '...'
    done
    ```
