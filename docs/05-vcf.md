# VCF Recipes

VCF (Variant Call Format) has a variable-length header block (lines starting with `##`
or `#CHROM`) followed by one variant per line with 8 fixed columns plus per-sample columns.

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  [FORMAT  SAMPLE...]
```

!!! tip "Test data"
    All recipes on this page work with [`docs/data/variants.vcf`](data/README.md):
    20 variants across chr1/chr2, 3 samples (NA12878, NA19238, NA20585), FORMAT `GT:DP:GQ:AD`,
    includes SNPs, INDELs, multi-allelic sites, PASS/LowQual filters, and SnpEff ANN annotations.

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

---

## Extract per-sample genotype fields

### Context

VCF column 9 (FORMAT) defines the field order for each sample column (e.g. `GT:DP:GQ:AD`).
Columns 10+ contain per-sample data in that same colon-delimited order.
To extract `GT`, `DP`, or `GQ` you must first parse FORMAT to find each field's position,
then index into each sample column accordingly.
The FORMAT field can vary between variants — do not assume a fixed position.

### Code

```awk
# Usage: awk -f genotypes.awk variants.vcf
# Output: CHROM POS SAMPLE GT DP GQ

/^#CHROM/ {
    for (i = 10; i <= NF; i++) sample_names[i] = $i
    next
}
/^#/ { next }
{
    # Find index of GT, DP, GQ in FORMAT (col 9)
    n = split($9, fmt, ":")
    gt_idx = dp_idx = gq_idx = 0
    for (i = 1; i <= n; i++) {
        if (fmt[i] == "GT") gt_idx = i
        if (fmt[i] == "DP") dp_idx = i
        if (fmt[i] == "GQ") gq_idx = i
    }
    # Emit one row per sample
    for (s = 10; s <= NF; s++) {
        split($s, gdata, ":")
        gt = gt_idx ? gdata[gt_idx] : "."
        dp = dp_idx ? gdata[dp_idx] : "."
        gq = gq_idx ? gdata[gq_idx] : "."
        print $1, $2, sample_names[s], gt, dp, gq
    }
}
```

### Explanation

- FORMAT is re-parsed per variant because different variant types can carry different fields
  (e.g. `GT:DP` vs `GT:DP:GQ:AD`). Assuming a fixed column index is a common bug.
- `gt_idx ? gdata[gt_idx] : "."` safely returns `.` when a field is absent from FORMAT.
- `split($s, gdata, ":")` reuses the `gdata` array — AWK clears it on each `split()` call.

### Variants

```awk
# Filter: only emit samples with DP >= min_depth
# Usage: awk -v min_depth=10 -f genotypes_filter.awk variants.vcf

/^#CHROM/ { for (i=10; i<=NF; i++) sample_names[i]=$i; next }
/^#/ { next }
{
    n = split($9, fmt, ":")
    gt_idx = dp_idx = 0
    for (i=1; i<=n; i++) {
        if (fmt[i]=="GT") gt_idx=i
        if (fmt[i]=="DP") dp_idx=i
    }
    for (s=10; s<=NF; s++) {
        split($s, g, ":")
        dp = dp_idx ? g[dp_idx]+0 : 0
        if (dp >= min_depth)
            print $1, $2, sample_names[s], (gt_idx ? g[gt_idx] : "."), dp
    }
}
```

---

## Handle multi-allelic sites

### Context

The ALT field can contain multiple comma-separated alleles (`A,T` or `A,T,G`).
Many downstream tools require biallelic records. AWK can split a multi-allelic
site into one record per ALT allele, or simply flag or skip such sites.

### Code

```awk
# Split multi-allelic sites — emit one line per ALT allele
/^#/ { print; next }
{
    n_alts = split($5, alts, ",")
    if (n_alts == 1) {
        print    # biallelic: pass through unchanged
    } else {
        for (i = 1; i <= n_alts; i++) {
            $5 = alts[i]
            print
        }
    }
}
```

!!! warning "Genotype fields are not updated"
    Per-sample GT values like `1/2` reference allele indices in the original multi-allelic
    record. After splitting, these indices no longer match the new single-ALT line.
    For correct splitting including genotype remapping, use `bcftools norm -m -any`.
    The AWK recipe above is appropriate for annotation or counting tasks where GT is not needed.

### Explanation

- `split($5, alts, ",")` returns the number of alternate alleles.
- Assigning `$5 = alts[i]` modifies the field and triggers `$0` reconstruction with `OFS`.
  Ensure `OFS="\t"` is set in `BEGIN` if you need tab-separated output.

### Variants

```awk
# Count multi-allelic sites
/^#/ { next }
split($5, a, ",") > 1 { multi++ }
END { print "Multi-allelic sites:", multi+0 }
```

```awk
# Skip multi-allelic sites entirely
/^#/ { print; next }
$5 !~ /,/ { print }
```

---

## Compute allele frequency from genotypes

### Context

The INFO `AF` field is set by the variant caller and may not reflect the actual allele
frequency in a subset of samples (e.g. after subsetting a cohort). Recompute AF directly
from the per-sample genotype columns.

### Code

```awk
# Output: CHROM POS REF ALT AF_computed N_missing
/^#CHROM/ {
    for (i=10; i<=NF; i++) sample_names[i]=$i
    next
}
/^#/ { next }
{
    # Find GT field index in FORMAT
    n = split($9, fmt, ":")
    gt_idx = 0
    for (i=1; i<=n; i++) if (fmt[i]=="GT") { gt_idx=i; break }

    ref_count = alt_count = missing = 0
    for (s = 10; s <= NF; s++) {
        split($s, gdata, ":")
        gt = gt_idx ? gdata[gt_idx] : "./."
        gsub(/\|/, "/", gt)          # normalise phased (0|1) → unphased (0/1)
        n_alleles = split(gt, alleles, "/")
        for (i = 1; i <= n_alleles; i++) {
            if      (alleles[i] == ".") missing++
            else if (alleles[i] == "0") ref_count++
            else                         alt_count++
        }
    }
    total = ref_count + alt_count
    af = (total > 0) ? sprintf("%.4f", alt_count / total) : "."
    printf "%s\t%d\t%s\t%s\t%s\t%d\n", $1, $2, $4, $5, af, missing
}
```

### Explanation

- `gsub(/\|/, "/", gt)` normalises phased genotypes (`0|1`) to unphased (`0/1`) before
  splitting on `/`.
- Missing alleles (`.`) are counted separately and excluded from the AF denominator.
- Multi-allelic sites: any non-`0` allele index is counted as ALT. For per-allele AF on
  multi-allelic sites, check the allele index and count per value.

### Variants

```awk
# Report per-sample heterozygosity rate alongside AF
/^#CHROM/ { for (i=10; i<=NF; i++) sn[i]=$i; next }
/^#/ { next }
{
    n_het = n_hom_ref = n_hom_alt = n_miss = 0
    n = split($9, fmt, ":"); gt_idx=0
    for (i=1; i<=n; i++) if (fmt[i]=="GT") gt_idx=i
    for (s=10; s<=NF; s++) {
        split($s, g, ":"); gt=g[gt_idx]
        gsub(/\|/,"/",gt)
        if      (gt=="./.")  n_miss++
        else if (gt=="0/0")  n_hom_ref++
        else if (gt~/^[1-9]\/[1-9]$/ && split(gt,a,"/")==2 && a[1]==a[2]) n_hom_alt++
        else                 n_het++
    }
    printf "%s\t%d\t%d\t%d\t%d\t%d\n", $1, $2, n_hom_ref, n_het, n_hom_alt, n_miss
}
```

---

## Extract ANN annotation field (SnpEff / VEP)

### Context

After running SnpEff or VEP, functional annotations are stored in the INFO field as
`ANN=allele|effect|impact|gene_name|gene_id|...` (SnpEff) or `CSQ=...` (VEP), with
multiple transcripts separated by commas and fields within each transcript separated
by pipes. Extracting the gene name and consequence class is a common post-annotation step.

### Code

```awk
# Extract gene name and effect from SnpEff ANN field [gawk]
# Output: CHROM POS REF ALT GENE EFFECT IMPACT
/^#/ { next }
{
    delete arr
    match($8, /ANN=([^;]+)/, arr)
    if (arr[1] == "") {
        print $1, $2, $4, $5, ".", ".", "."
        next
    }
    # ANN can have multiple transcripts, comma-separated — use the first
    n = split(arr[1], transcripts, ",")
    split(transcripts[1], f, "|")
    # SnpEff ANN fields: allele|effect|impact|gene_name|gene_id|feature_type|...
    print $1, $2, $4, $5, f[4], f[2], f[3]
}
```

**POSIX version — no `match()` capture groups:**

```awk
/^#/ { next }
{
    ann = ""
    n = split($8, info_fields, ";")
    for (i = 1; i <= n; i++) {
        if (substr(info_fields[i], 1, 4) == "ANN=") {
            ann = substr(info_fields[i], 5)
            break
        }
    }
    if (ann == "") { print $1, $2, $4, $5, ".", ".", "."; next }
    # Use first transcript only
    n_tx = split(ann, transcripts, ",")
    split(transcripts[1], f, "|")
    print $1, $2, $4, $5, f[4], f[2], f[3]
}
```

### Explanation

- `match($8, /ANN=([^;]+)/, arr)` captures everything from `ANN=` to the next `;` or
  end of field. `arr[1]` holds the full ANN value.
- SnpEff can annotate multiple transcripts per variant: `ANN=T|effect1|...,T|effect2|...`.
  The recipe uses the first transcript. For the most severe consequence, sort by impact
  level (HIGH > MODERATE > LOW > MODIFIER) — typically done with `sort -t'|' -k3`.
- VEP uses `CSQ=` instead of `ANN=` and has a different pipe-separated field order defined
  in the `##CSQ` header line. Adjust the field index accordingly.

### Variants

```awk
# Count variants by impact level [gawk]
/^#/ { next }
{
    delete arr
    match($8, /ANN=([^;]+)/, arr)
    if (!arr[1]) next
    split(arr[1], tx, ",")
    split(tx[1], f, "|")
    impact[f[3]]++    # f[3] = impact: HIGH, MODERATE, LOW, MODIFIER
}
END {
    for (imp in impact)
        print imp, impact[imp]
}
```
