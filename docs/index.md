# AWK Bioinformatics Cookbook

Practical AWK recipes for bioinformatics — organised by problem, not by AWK feature.

Each recipe answers a real question:
*"How do I filter FASTQ reads by length?"*,
*"How do I extract a VCF INFO field?"*,
*"How do I join two TSV files by sample ID?"*

For every recipe: context, code, line-by-line explanation, and variants.

---

## POSIX AWK vs gawk

Most recipes work with any POSIX-compliant AWK (`awk`, `mawk`, `nawk`).
Recipes marked **[gawk]** require [GNU AWK](https://www.gnu.org/software/gawk/) and use
`ARGIND`, 3-argument `match()`, `gensub()`, `PROCINFO`, or `switch`.

```bash
awk --version    # GNU Awk x.y.z → you have gawk
```

On most Linux HPC clusters, `awk` **is** gawk.
On macOS, the default `awk` is BSD AWK — install gawk with `brew install gawk`.

---

## Recipes

| Category | Topics |
|----------|--------|
| [Fundamentals](01-fundamentals.md) | FS/OFS/RS, BEGIN/END, printf, functions, error handling |
| [FASTQ](02-fastq.md) | Block navigation (`NR%4`), length filter, dedup, read extraction |
| [FASTA](03-fasta.md) | Headers, single-line reformat, length filter, rename, N50 |
| [BED / GFF](04-bed-gff.md) | Feature filter, coordinate conversion, sizes, GFF attributes |
| [VCF](05-vcf.md) | Header handling, QUAL/FILTER, INFO field extraction |
| [SAM-derived TSV](06-sam-tsv.md) | idxstats, MAPQ filter, coverage depth |
| [Two-file Joins](07-joins.md) | `NR==FNR`, left-join, `ARGIND` [gawk], dynamic FS |
| [Reports & Aggregation](08-reports.md) | Group-by, min/max/mean, pivot, sorted output [gawk] |
| [Advanced Patterns](09-advanced.md) | Multi-dim arrays, `match()` capture [gawk], `gensub()`, `$1=$1`, `delete` |
