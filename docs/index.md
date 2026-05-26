# AWK Bioinformatics Cookbook

## About this cookbook

This cookbook is the result of nearly 20 years of bioinformatics practice — from
diagnostic sequencing pipelines in clinical labs to large-scale population genomics
projects. I'm [Christophe Demay](https://chrstphdm.github.io/), a senior
bioinformatics consultant and workflow engineer. Over the years, I've maintained
pipelines across 12+ clinical areas, built HLA typing workflows, and processed
datasets of 40,000+ whole genomes and much more.

AWK has been a constant throughout all of this. Not because it's trendy — but because
when you need to filter a 50 GB VCF at 3 AM on an HPC node, AWK is already there,
it starts in milliseconds, and the one-liner you write is often faster than the Python
script you'd spend 15 minutes setting up.

**Why this guide exists:** I got tired of re-Googling the same AWK patterns, copying
snippets between Slack threads, and explaining the `NR==FNR` idiom to every new team
member. This cookbook consolidates what I've learned — and what the official AWK
documentation doesn't teach — into one place, organised by the problems you actually
face in bioinformatics.

**What makes it different:**

- **Problem-driven, not feature-driven** — recipes answer real questions, not "here is
  how `gsub` works"
- **Battle-tested patterns** — from clinical NGS, HLA typing, RNA-seq, and population
  genomics pipelines
- **Honest about AWK's limits** — a clear decision matrix tells you when Python, R, or
  a specialised tool is the better choice
- **Executable** — every recipe runs on the included synthetic test data
- **Nextflow-ready** — dedicated chapter on AWK in Nextflow modules, quoting rules,
  `templates/` directory patterns

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
| [Fundamentals](01-fundamentals.md) | FS/OFS/RS, BEGIN/END, printf, `.awk` script files, functions, error handling |
| [FASTQ](02-fastq.md) | Block navigation (`NR%4`), length filter, dedup, read extraction |
| [FASTA](03-fasta.md) | Headers, single-line reformat, length filter, rename, N50 |
| [BED / GFF](04-bed-gff.md) | Feature filter, coordinate conversion, sizes, GFF attributes, interval merge |
| [VCF](05-vcf.md) | Header handling, QUAL/FILTER, INFO/ANN extraction, genotypes, multi-allelic |
| [SAM-derived TSV](06-sam-tsv.md) | idxstats, MAPQ filter, coverage depth |
| [Two-file Joins](07-joins.md) | `NR==FNR`, left-join, `ARGIND` [gawk], dynamic FS |
| [Reports & Aggregation](08-reports.md) | Group-by, min/max/mean, pivot, sorted output [gawk] |
| [Advanced Patterns](09-advanced.md) | Multi-dim arrays, `match()` capture [gawk], `gensub()`, `$1=$1`, `delete` |
| [getline](10-getline.md) | Read from file/command, coprocess [gawk], common pitfalls |
| [GTF Annotation](11-gtf.md) | gene_id/gene_name, exon counts, transcript lengths, UTRs, Ensembl vs GENCODE |
| [RNA-seq Counts](12-rnaseq.md) | HTSeq-count, featureCounts, CPM, merge N files, low-count filter |
| [Multi-sample Patterns](13-multi-sample.md) | Aggregate cohorts, outlier detection, missingness, QC summary |
| [Testing AWK Scripts](14-testing.md) | diff-based tests, bats, CI integration, POSIX compatibility |
| [AWK in Nextflow](15-nextflow.md) | `script:` blocks, quoting rules, params, nf-core modules, debugging |
| [Workflows](16-workflows.md) | FASTQ QC, VCF+BED annotation, RNA-seq pipeline, cohort report |
