# BED / GFF Recipes

**BED** is tab-separated with 0-based, half-open coordinates: `chrom start end [name score strand ...]`.
**GFF/GFF3** is tab-separated with 1-based, closed coordinates: `seqname source feature start end score strand frame attributes`.
GTF follows the same 9-field structure as GFF.

Both formats use `#` for comment/header lines.

---

## Filter by chromosome

### Context

Extract records for a specific chromosome or a set of chromosomes.
Useful before per-chromosome analysis.

### Code

```awk
# Single chromosome
awk '$1 == "chr1"' regions.bed

# A set of chromosomes loaded from a file
# File 1: chroms.txt (one name per line)
# File 2: regions.bed
awk 'NR==FNR { keep[$1]=1; next } $1 in keep' chroms.txt regions.bed

# Exclude unplaced contigs (keep only chrN and chrX/Y/M)
awk '$1 ~ /^chr([0-9]+|[XYM])$/' regions.bed
```

### Explanation

- `$1 == "chr1"` is a **pattern** — the default action (`print`) fires when it matches.
- `$1 ~ /^chr([0-9]+|[XYM])$/` uses a regex to match standard chromosome names.
- The `NR==FNR` join loads a whitelist from a file; see [Two-file joins](07-joins.md) for the full pattern.

---

## Filter GFF by feature type

### Context

Extract only genes, or only exons, from a GFF3 file. Skip comment lines.

### Code

```awk
# Skip comments, keep only "gene" features
/^#/ { next }
$3 == "gene" { print }
```

```awk
# Keep multiple feature types
/^#/ { next }
$3 == "exon" || $3 == "CDS" { print }
```

### Explanation

- GFF column 3 is the **feature type** (`gene`, `mRNA`, `exon`, `CDS`, `UTR`, ...).
- `/^#/ { next }` discards comment lines before any field parsing happens —
  important because comment lines don't have 9 fields and will cause wrong output.

---

## Reformat coordinates (0-based to 1-based)

### Context

BED is 0-based half-open (`start` is 0-indexed, `end` is exclusive).
VCF and GFF are 1-based closed. Convert when combining sources.

### Code

```awk
# BED (0-based) → 1-based closed interval
/^#/ { print; next }
{ print $1, $2 + 1, $3, $4 }
```

```awk
# 1-based closed → BED 0-based half-open
/^#/ { print; next }
{ print $1, $2 - 1, $3, $4 }
```

### Explanation

- Adding 1 to the start converts from 0-based to 1-based; the end coordinate is the
  same in both conventions (exclusive 0-based = inclusive 1-based).
- Header/comment lines are passed through unchanged with `print; next`.

---

## Compute feature sizes from BED

### Context

What is the total covered size per chromosome? What is the mean exon length?

### Code

```awk
/^#/ { next }
{
    chrom  = $1
    size   = $3 - $2        # BED: end - start (0-based half-open)
    total[chrom] += size
    count[chrom]++
}
END {
    printf "%-10s %12s %8s %10s\n", "CHROM", "TOTAL_BP", "FEATURES", "MEAN_SIZE"
    for (c in total)
        printf "%-10s %12d %8d %10.1f\n", c, total[c], count[c], total[c]/count[c]
}
```

### Explanation

- For BED, feature size = `$3 - $2` (end minus start in 0-based half-open coordinates).
- For GFF, feature size = `$5 - $4 + 1` (end minus start + 1, 1-based closed).
- Results are stored per chromosome in associative arrays; the `END` block prints the summary.

---

## Extract GFF attribute field

### Context

GFF column 9 is a semicolon-separated list of key=value attributes:
`ID=gene:ENSG00000...; Name=BRCA2; biotype=protein_coding`.
Extract a specific attribute by name.

### Code

```awk
# Extract the "ID" attribute from GFF3
/^#/ { next }
$3 == "gene" {
    # split on ";" to get individual attributes
    n = split($9, attrs, ";")
    gene_id = ""
    for (i = 1; i <= n; i++) {
        # split each attr on "=" — trim leading spaces first
        gsub(/^ +/, "", attrs[i])
        if (split(attrs[i], kv, "=") == 2 && kv[1] == "ID") {
            gene_id = kv[2]
            break
        }
    }
    print $1, $4, $5, gene_id
}
```

### Variants

```awk
# Simpler: use match() to capture the value directly [gawk]
/^#/ { next }
$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, arr)   # GTF-style quoted attribute
    print $1, $4, $5, arr[1]
}
```

### Explanation

- GFF3 uses `key=value` pairs; GTF uses `key "value"` pairs with quotes.
- The loop approach works with any AWK; the `match()` capture approach requires gawk.
- `gsub(/^ +/, "", attrs[i])` strips leading spaces that appear after a semicolon split.

!!! note "GTF in depth"
    GTF (the format used by RNA-seq tools like STAR, featureCounts, StringTie) has a
    more complex attribute structure and a gene → transcript → exon hierarchy.
    See [GTF Annotation](11-gtf.md) for a dedicated chapter.

---

## Merge overlapping BED intervals

### Context

Before computing total coverage, overlapping intervals from the same chromosome must
be merged — otherwise a position covered by two overlapping features is counted twice.

!!! warning "Input must be sorted"
    This recipe assumes the BED file is sorted by chromosome then start position.
    Sort first: `sort -k1,1 -k2,2n regions.bed | awk -f merge_bed.awk`

### Code

```awk
# Usage: sort -k1,1 -k2,2n regions.bed | awk -f merge_bed.awk
# Merges overlapping/adjacent intervals on the same chromosome

/^#/ { print; next }

FNR == 1 {
    chrom = $1; start = $2; end = $3
    next
}

$1 == chrom && $2 <= end {
    # Overlapping or adjacent: extend current interval if needed
    if ($3 > end) end = $3
    next
}

{
    # New chromosome or non-overlapping interval: emit current, start new
    print chrom, start, end
    chrom = $1; start = $2; end = $3
}

END {
    if (chrom != "") print chrom, start, end
}
```

```bash
# With OFS for tab-separated output
sort -k1,1 -k2,2n docs/data/regions.bed \
  | awk 'BEGIN{OFS="\t"} /^#/{print;next}
         FNR==1{c=$1;s=$2;e=$3;next}
         $1==c && $2<=e{if($3>e)e=$3;next}
         {print c,s,e; c=$1;s=$2;e=$3}
         END{if(c!="")print c,s,e}'
```

### Explanation

- The key condition is `$2 <= end`: if the next interval's start is ≤ the current end,
  they overlap (BED is 0-based half-open, so adjacent intervals have `$2 == end`).
- `if ($3 > end) end = $3` extends the current merged interval if the new one reaches
  further right.
- The `END` block ensures the last accumulated interval is printed.
