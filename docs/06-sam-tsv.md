# SAM-derived TSV Recipes

SAM files are large and binary-adjacent — use `samtools` to extract data, then AWK
to process the resulting TSV. Common sources:

- `samtools flagstat` → alignment summary stats
- `samtools idxstats` → per-reference read counts and lengths
- `samtools depth` → per-base coverage (3 columns: chrom, pos, depth)
- `samtools view -F 4 | cut -f1,2,3,4,5` → per-read TSV extracts

---

## Per-sample read counts from idxstats

### Context

`samtools idxstats` outputs one line per reference sequence:
`ref_name  ref_len  mapped_reads  unmapped_reads`.
Summarise total mapped and unmapped reads across all chromosomes.

### Code

```awk
# Usage: samtools idxstats sample.bam | awk -f idxstats_summary.awk

{
    mapped   += $3
    unmapped += $4
}
END {
    total = mapped + unmapped
    printf "Mapped:   %d (%.1f%%)\n", mapped,   (total > 0 ? mapped/total*100   : 0)
    printf "Unmapped: %d (%.1f%%)\n", unmapped, (total > 0 ? unmapped/total*100 : 0)
    printf "Total:    %d\n", total
}
```

### Variants

```awk
# Per-chromosome summary, skip the unmapped catch-all line ("*")
$1 != "*" {
    printf "%s\t%d\t%d\n", $1, $3, $4
}
```

---

## Filter by mapping quality

### Context

`samtools view` TSV output (via `-o` or pipe): column 5 is MAPQ.
Remove low-MAPQ reads before downstream processing.

```
read_name  flag  ref  pos  mapq  cigar  rnext  pnext  tlen  seq  qual  [tags...]
```

### Code

```awk
# Usage: samtools view -h sample.bam | awk -v min_mapq=20 -f filter_mapq.awk

/^@/ { print; next }       # pass SAM header lines through
$5 + 0 >= min_mapq { print }
```

### Explanation

- SAM header lines start with `@`; pass them through unchanged.
- Column 5 is MAPQ (0–255); 0 = unmapped or multimapper; 20 means ≥1% chance of wrong
  mapping; 30 means ≥0.1%; 60 is typically the max for BWA-MEM.
- `$5 + 0` forces numeric comparison in case of non-numeric values.

---

## Summarise coverage per region

### Context

`samtools depth` outputs chrom, position, depth — one line per base.
Compute mean coverage over each region from a BED file, or just over all positions.

### Code

```awk
# Mean depth from samtools depth output
# Usage: samtools depth -a sample.bam | awk -f mean_depth.awk

{
    bases++
    total_depth += $3
    if ($3 == 0) zero_cov++
}
END {
    if (bases > 0) {
        printf "Mean depth:        %.2f\n", total_depth / bases
        printf "Breadth (> 0x):    %.2f%%\n", (bases - zero_cov) / bases * 100
        printf "Positions:         %d\n", bases
    }
}
```

```awk
# Per-chromosome mean depth
{
    depth[$1]  += $3
    bases[$1]++
}
END {
    printf "%-10s %10s %10s\n", "CHROM", "MEAN_COV", "BASES"
    for (c in depth)
        printf "%-10s %10.2f %10d\n", c, depth[c]/bases[c], bases[c]
}
```

### Variants

```awk
# Coverage at specific thresholds (1x, 10x, 20x, 30x)
{
    bases++
    depth = $3
    if (depth >= 1)  d1x++
    if (depth >= 10) d10x++
    if (depth >= 20) d20x++
    if (depth >= 30) d30x++
}
END {
    printf "1x:  %.1f%%\n10x: %.1f%%\n20x: %.1f%%\n30x: %.1f%%\n",
        d1x/bases*100, d10x/bases*100, d20x/bases*100, d30x/bases*100
}
```

### Explanation

- `-a` flag in `samtools depth` includes zero-coverage positions; without it, uncovered
  bases are absent and breadth cannot be computed correctly.
- `(bases - zero_cov) / bases * 100` gives the fraction of positions with at least 1x.

## References

- [SAM/BAM format specification (hts-specs)](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [samtools documentation](https://www.htslib.org/doc/samtools.html)
