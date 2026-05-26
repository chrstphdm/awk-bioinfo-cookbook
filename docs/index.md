# AWK Bioinformatics Cookbook

Practical AWK recipes for bioinformatics — organised by problem, not by AWK feature.

Each recipe answers a real question:
*"How do I filter FASTQ reads by length?"*,
*"How do I extract a VCF INFO field?"*,
*"How do I join two TSV files by sample ID?"*

For every recipe: context, code, line-by-line explanation, and variants.

!!! tip "New here?"
    - **Python or R user?** Start with [Why AWK?](00-why-awk.md) — includes a decision
      matrix, side-by-side comparisons, and honest advice on when to use each tool.
    - **AWK user?** Jump straight to [Fundamentals](01-fundamentals.md) or any format chapter.
    - **Nextflow user?** See [AWK in Nextflow](15-nextflow.md) for quoting rules and module patterns.

---

## Test Data

All recipes use synthetic files in [`docs/data/`](data/README.md).
Clone the repo and run recipes directly:

```bash
git clone https://github.com/chrstphdm/awk-bioinfo-cookbook
cd awk-bioinfo-cookbook/docs/data

# Quick sanity checks
awk 'NR%4==1' reads.fastq | wc -l        # → 20 reads
awk '!/^#/' variants.vcf | wc -l         # → 20 variants
awk '!/^#/' regions.bed | wc -l          # → 15 BED features
```

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
See [Why AWK? → Implementations](00-why-awk.md#awk-implementations-which-one-do-you-have)
for a full comparison table.

---

## Recipes

| Category | Topics |
|----------|--------|
| [Why AWK?](00-why-awk.md) | Implementations (gawk/mawk/nawk/bioawk), AWK vs Python vs R, Parquet |
| [Fundamentals](01-fundamentals.md) | FS/OFS/RS, BEGIN/END, printf, functions, error handling |
| [FASTQ](02-fastq.md) | Block navigation (`NR%4`), length filter, dedup, read extraction |
| [FASTA](03-fasta.md) | Headers, single-line reformat, length filter, rename, N50 |
| [BED / GFF](04-bed-gff.md) | Feature filter, coordinate conversion, sizes, GFF attributes, interval merge |
| [VCF](05-vcf.md) | Header handling, QUAL/FILTER, INFO/ANN extraction, genotypes, multi-allelic |
| [SAM-derived TSV](06-sam-tsv.md) | idxstats, MAPQ filter, coverage depth |
| [Two-file Joins](07-joins.md) | `NR==FNR`, left-join, `ARGIND` [gawk], dynamic FS |
| [Reports & Aggregation](08-reports.md) | Group-by, min/max/mean, pivot, sorted output [gawk] |
| [Advanced Patterns](09-advanced.md) | Multi-dim arrays, `match()` capture [gawk], `gensub()`, `$1=$1`, `delete` |
| [getline](10-getline.md) | Read from file/command, coprocess [gawk], common pitfalls |
| [Testing AWK Scripts](10-testing.md) | diff-based tests, bats, CI integration, POSIX compatibility |
| [GTF Annotation](11-gtf.md) | gene_id/gene_name, exon counts, transcript lengths, UTRs, Ensembl vs GENCODE |
| [RNA-seq Counts](12-rnaseq.md) | HTSeq-count, featureCounts, CPM, merge N files, low-count filter |
| [Multi-sample Patterns](13-multi-sample.md) | Aggregate cohorts, outlier detection, missingness, QC summary |
| [AWK in Nextflow](15-nextflow.md) | `script:` blocks, quoting rules, params, nf-core modules, debugging |
| [Workflows](16-workflows.md) | FASTQ QC, VCF+BED annotation, RNA-seq pipeline, cohort report |
