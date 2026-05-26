# RNA-seq Count Recipes

After alignment and quantification, RNA-seq data lives in count matrices: one row per
gene, one column per sample. Two tools dominate: **HTSeq-count** and **featureCounts**.
Their output formats differ enough that separate recipes are needed for each.

!!! tip "Test data"
    - HTSeq format: [`docs/data/htseq_counts.tsv`](data/README.md) — 15 genes + `__` summary lines
    - featureCounts format: [`docs/data/featurecounts.tsv`](data/README.md) — 15 genes × 3 samples

---

## HTSeq-count format

HTSeq-count produces a 2-column TSV with no header. The last 5 lines are diagnostic
summary rows prefixed with `__`:

```
GENE1   1245
GENE2   834
...
__no_feature    18432
__ambiguous     342
__too_low_aQual 127
__not_aligned   4821
__alignment_not_unique  891
```

These `__` lines must be excluded before any downstream analysis.

---

## Remove `__` summary lines and extract counts

### Context

The simplest and most common operation: strip the diagnostic rows, keep only the gene
count rows.

### Code

```awk
# Remove __summary lines
$1 !~ /^__/ { print }
```

```bash
# Typical pipeline: filter then pass to downstream tool
awk '$1 !~ /^__/' htseq_counts.tsv | sort -k2,2rn | head -20
```

### Variants

```awk
# Capture __no_feature as a QC metric while filtering
$1 == "__no_feature"          { printf "QC no_feature: %d\n", $2 > "/dev/stderr"; next }
$1 ~ /^__/                    { next }
{ print }
```

```awk
# Report all summary stats to stderr, then output clean counts
/^__/ {
    gsub(/^__/, "", $1)
    printf "SUMMARY %-30s %d\n", $1":", $2 > "/dev/stderr"
    next
}
{ print }
```

---

## featureCounts format

featureCounts output has a 2-line header followed by data:

```
## featureCounts output, version 2.0.3
Geneid  Chr  Start  End  Strand  Length  sample1  sample2  sample3
GENE1   chr1 1001   2500 +       1500    1245     980      1102
...
```

- Line 1: `##` comment (command used)
- Line 2: column headers — first 6 are annotation metadata (`Geneid Chr Start End Strand Length`),
  columns 7+ are sample names
- Data lines: gene ID + annotation metadata + counts per sample

---

## Extract gene_id and counts from featureCounts

### Context

Strip the 5 annotation metadata columns (Chr, Start, End, Strand, Length), keeping
only gene ID and count columns — the format expected by most downstream R/Python tools.

### Code

```awk
/^##/ { next }              # skip command-line comment
NR == 2 {
    # Header line: print gene_id + sample names (cols 7+)
    printf "gene_id"
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
    next
}
{
    # Data: print gene_id + counts (cols 7+)
    printf "%s", $1
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
}
```

```bash
# One-liner version
awk '/^##/{next} NR==2{printf "gene_id"; for(i=7;i<=NF;i++) printf "\t"$i; print""; next}
     {printf $1; for(i=7;i<=NF;i++) printf "\t"$i; print""}' featurecounts.tsv
```

---

## Compute CPM (Counts Per Million)

### Context

Normalise raw counts by library size to make samples comparable. CPM = (count /
total_counts_in_sample) × 10⁶. Requires a two-pass approach: first sum columns,
then divide.

### Code

```awk
# Input: clean count matrix (gene_id + counts, with header row)
# Output: same dimensions, counts replaced by CPM values

NR == 1 {
    # Store header and number of columns
    n_samples = NF - 1
    print
    next
}
{
    # Store data and accumulate column sums
    for (i = 2; i <= NF; i++) {
        col_sum[i] += $i
        row_val[NR][i] = $i
    }
    row_gene[NR] = $1
}
END {
    for (r = 2; r <= NR; r++) {
        printf "%s", row_gene[r]
        for (i = 2; i <= NF; i++)
            printf "\t%.4f", (row_val[r][i] / col_sum[i]) * 1e6
        printf "\n"
    }
}
```

!!! warning "Memory usage"
    This recipe stores the entire matrix in memory. For matrices with >20,000 genes and
    >100 samples, consider computing column sums in a first pass (`awk '{for(i=2;i<=NF;i++) s[i]+=$i}'`),
    saving them to a file, then doing the CPM calculation in a second pass.

### Variants

```awk
# RPKM: also normalise by gene length (requires a length column or separate file)
# Input: featureCounts format (col 6 = Length, cols 7+ = counts)

/^##/ { next }
NR == 2 {
    printf "gene_id"
    for (i = 7; i <= NF; i++) { printf "\t%s", $i; col_sum[i] = 0 }
    printf "\n"
    next
}
{
    for (i = 7; i <= NF; i++) col_sum[i] += $i
    gene_len[NR] = $6
    for (i = 7; i <= NF; i++) row_val[NR][i] = $i
    row_gene[NR] = $1
}
END {
    for (r = 3; r <= NR; r++) {
        printf "%s", row_gene[r]
        for (i = 7; i <= NF; i++) {
            # RPKM = (count / (library_size_millions * gene_length_kb))
            rpkm = (row_val[r][i] / ((col_sum[i]/1e6) * (gene_len[r]/1e3)))
            printf "\t%.4f", rpkm
        }
        printf "\n"
    }
}
```

---

## Merge N HTSeq-count files into a matrix

### Context

You have one HTSeq-count file per sample (`NA12878.counts`, `NA19238.counts`, ...) and
need a single gene × sample matrix. This is one of the most practically useful AWK
recipes in an RNA-seq workflow.

### Code

```awk
# Usage: awk -f merge_htseq.awk NA12878.counts NA19238.counts NA20585.counts > matrix.tsv
# [gawk — uses FILENAME; on POSIX AWK, pass sample names via -v or use FNR==1 counter]

# Track file index
FNR == 1 {
    file_idx++
    # Strip directory and extension from filename for clean column header
    fname = FILENAME
    sub(/.*\//, "", fname)   # remove path
    sub(/\.[^.]+$/, "", fname)   # remove last extension
    sample_names[file_idx] = fname
}

$1 !~ /^__/ {
    counts[$1][file_idx] = $2
    genes[$1] = 1
}

END {
    # Header
    printf "gene_id"
    for (i = 1; i <= file_idx; i++) printf "\t%s", sample_names[i]
    printf "\n"
    # Data — one row per gene
    for (gene in genes) {
        printf "%s", gene
        for (i = 1; i <= file_idx; i++)
            printf "\t%s", ((gene, i) in counts ? counts[gene][i] : "0")
        printf "\n"
    }
}
```

!!! tip "For 100+ samples, prefer `paste`"
    AWK loads all counts into memory. For large cohorts, `paste` is faster and uses
    less memory — it just needs a header fix:

    ```bash
    # Merge counts columns with paste, then prepend gene_id header
    samples=(NA12878 NA19238 NA20585)   # expand to all samples

    # Extract gene column from first file
    awk '$1 !~ /^__/' "${samples[0]}.counts" | cut -f1 > gene_ids.txt

    # Extract count column from each file (skip __ lines)
    for s in "${samples[@]}"; do
        awk '$1 !~ /^__/ {print $2}' "${s}.counts" > "${s}.col"
    done

    # Paste together and add header
    paste gene_ids.txt $(printf '%s.col ' "${samples[@]}") \
      | awk -v hdr="gene_id\t$(printf '%s\t' "${samples[@]}" | sed 's/\t$//')" \
            'BEGIN{print hdr} {print}' > matrix.tsv
    ```

### Explanation

- `FILENAME` gives the current input file's name (gawk and most AWK implementations).
- `(gene, i) in counts` uses SUBSEP to test existence in a 2D array without creating
  an empty entry. Equivalent to `gene in counts && i in counts[gene]`.
- Genes are accumulated in `genes[$1]=1` to collect the full universe across all files —
  a gene absent in one sample will get `0` in the output.

---

## Filter low-count genes

### Context

Genes with very low counts across all samples are typically noise. A common filter:
keep genes with at least `min_count` counts in at least `min_samples` samples.

### Code

```awk
# Usage: awk -v min_count=10 -v min_samples=2 -f filter_lowcount.awk matrix.tsv
# Input: tab-separated matrix with header (gene_id + count columns)

NR == 1 { print; next }    # pass header through
{
    n_pass = 0
    for (i = 2; i <= NF; i++)
        if ($i + 0 >= min_count) n_pass++
    if (n_pass >= min_samples) print
}
```

```bash
# Quick check: how many genes survive the filter?
awk -v min_count=10 -v min_samples=2 '
    NR==1{next}
    { n=0; for(i=2;i<=NF;i++) if($i+0>=min_count) n++ }
    n>=min_samples{kept++} END{print "Kept:", kept+0}
' matrix.tsv
```

### Variants

```awk
# Report filter stats to stderr while filtering
NR == 1 { n_samples = NF - 1; print; next }
{
    total++
    n_pass = 0
    for (i = 2; i <= NF; i++) if ($i + 0 >= min_count) n_pass++
    if (n_pass >= min_samples) { kept++; print }
}
END {
    printf "Filtered: %d/%d genes kept (%.1f%%)\n",
           kept+0, total, (kept+0)/total*100 > "/dev/stderr"
}
```
