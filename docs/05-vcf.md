# VCF Recipes

VCF (Variant Call Format) has a variable-length header block (lines starting with `##`
or `#CHROM`) followed by one variant per line with 8 fixed columns plus per-sample columns.

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  [FORMAT  SAMPLE...]
```

---

## Skip header lines

### Context

VCF headers contain metadata (`##`) and the sample name line (`#CHROM`).
Most AWK processing needs to skip or handle them separately.

### Code

```awk
# Skip all header lines
/^#/ { next }
{ print $1, $2, $4, $5 }    # CHROM, POS, REF, ALT
```

```awk
# Pass headers through, process variant lines
/^#/ { print; next }
{ ... }
```

```awk
# Extract sample names from the #CHROM header
/^#CHROM/ {
    for (i = 10; i <= NF; i++) samples[i] = $i
    next
}
```

### Explanation

- `##` lines are metadata (format version, contig lengths, INFO/FORMAT definitions).
- `#CHROM` is the column header line — sample names start at column 10.
- Columns 1-8 are fixed; column 9 is the FORMAT definition; columns 10+ are per-sample.

---

## Filter variants by QUAL score

### Context

Remove low-confidence calls. QUAL is column 6; `.` means missing.

### Code

```awk
# Usage: awk -v min_qual=30 -f filter_qual.awk variants.vcf

/^#/ { print; next }
$6 != "." && $6 + 0 >= min_qual { print }
```

### Explanation

- `$6 + 0` forces numeric comparison — if `$6` is `.`, `$6 + 0` evaluates to 0.
- The explicit `$6 != "."` guard is clearer and handles the missing-value case
  independently of the threshold.

---

## Filter by FILTER column

### Context

Keep only variants that passed all filters. FILTER is column 7; `PASS` means all filters
passed; `.` means no filters applied; anything else is a filter name that failed.

### Code

```awk
# Keep PASS variants only
/^#/ { print; next }
$7 == "PASS" { print }
```

```awk
# Keep PASS and unfiltered (".")
/^#/ { print; next }
$7 == "PASS" || $7 == "." { print }
```

```awk
# Exclude a specific filter
/^#/ { print; next }
$7 !~ /LowQual/ { print }
```

---

## Extract a field from the INFO column

### Context

INFO (column 8) is a semicolon-separated list of `key=value` pairs (and flags).
Extract the value of a specific key — e.g. `AF` (allele frequency), `DP` (depth).

### Code

```awk
# Extract AF from INFO
/^#/ { next }
{
    n = split($8, info_fields, ";")
    af = "."
    for (i = 1; i <= n; i++) {
        if (split(info_fields[i], kv, "=") == 2 && kv[1] == "AF") {
            af = kv[2]
            break
        }
    }
    print $1, $2, $4, $5, af
}
```

**gawk version — more concise:**

```awk
/^#/ { next }
{
    match($8, /AF=([^;]+)/, arr)
    print $1, $2, $4, $5, (arr[1] != "" ? arr[1] : ".")
}
```

### Explanation

- The loop version works with any AWK; the `match()` version requires gawk.
- `match($8, /AF=([^;]+)/, arr)` captures everything between `AF=` and the next `;`
  (or end of string) into `arr[1]`.
- Flags (INFO entries without `=`) are skipped by the `split(...) == 2` check.

---

## Count variants per chromosome

### Context

Quick summary: how many variants on each chromosome?

### Code

```awk
/^#/ { next }
{ counts[$1]++ }
END {
    for (chrom in counts)
        print chrom, counts[chrom]
}
```

```bash
# Sorted output
awk '/^#/{next} {c[$1]++} END{for(k in c) print k, c[k]}' variants.vcf \
  | sort -V -k1,1
```

### Variants

```awk
# Count by variant type (SNP vs INDEL)
/^#/ { next }
{
    is_indel = (length($4) != 1 || length($5) != 1)
    type = is_indel ? "INDEL" : "SNP"
    counts[$1][type]++
}
END {
    for (c in counts)
        print c, counts[c]["SNP"]+0, counts[c]["INDEL"]+0
}
```

### Explanation

- `counts[$1]++` uses the chromosome name as array key.
- SNP vs INDEL: a substitution is SNP if both REF and ALT are exactly 1 character long.
  Multi-allelic sites (`ALT` contains `,`) are not handled here — add a `split($5, alts, ",")`
  loop for strict classification.
