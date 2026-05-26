# AWK Bioinformatics Cookbook

Practical AWK recipes for bioinformatics — organised by problem, not by AWK feature.

Each recipe answers a real question you'd type into a search engine:
*"How do I filter FASTQ reads by length?"*,
*"How do I extract a VCF INFO field?"*,
*"How do I compute allele frequency from genotypes?"*,
*"How do I merge HTSeq count files into a matrix?"*

For every recipe you get: context, code, line-by-line explanation, and variants.

**Browse the full cookbook:** [chrstphdm.github.io/awk-bioinfo-cookbook](https://chrstphdm.github.io/awk-bioinfo-cookbook)

---

## Quick start

```bash
git clone https://github.com/chrstphdm/awk-bioinfo-cookbook
cd awk-bioinfo-cookbook/docs/data

# Filter a VCF by QUAL ≥ 30, keeping headers
awk '/^#/ || ($6 != "." && $6+0 >= 30)' variants.vcf | head

# Mean coverage per chromosome from samtools depth output
awk '{ sum[$1]+=$3; n[$1]++ } END { for(c in sum) printf "%s\t%.1f\n",c,sum[c]/n[c] }' depth.tsv
```

All recipes use synthetic test data in `docs/data/` — clone and run immediately.

---

## POSIX AWK vs gawk

Most recipes work with any POSIX-compliant AWK (`awk`, `mawk`, `nawk`).
Recipes marked **[gawk]** require [GNU AWK](https://www.gnu.org/software/gawk/).
See [Why AWK?](docs/00-why-awk.md) for a full comparison of implementations
(gawk, mawk, nawk, bioawk) and an honest decision matrix: AWK vs Python vs R.

```bash
awk --version          # GNU Awk x.y.z  → gawk
```

On most Linux HPC clusters, `awk` is gawk. On macOS, install gawk with `brew install gawk`.

---

## Recipes

### [Why AWK?](docs/00-why-awk.md)
AWK implementations compared, decision matrix (AWK vs Python vs R), Parquet/columnar formats, when AWK becomes the wrong tool.

### [1. Fundamentals](docs/01-fundamentals.md)
FS/OFS/RS, BEGIN/END, printf, `.awk` script files, functions, error handling.

### [2. FASTQ](docs/02-fastq.md)
Block navigation (`NR%4`), length filter, dedup, read extraction.

### [3. FASTA](docs/03-fasta.md)
Headers, single-line reformat, length filter, rename, N50.

### [4. BED / GFF](docs/04-bed-gff.md)
Feature filter, coordinate conversion, sizes, GFF attributes, interval merge.

### [5. VCF](docs/05-vcf.md)
Header handling, QUAL/FILTER, INFO extraction, per-sample genotypes, multi-allelic sites, allele frequency, SnpEff ANN parsing.

### [6. SAM-derived TSV](docs/06-sam-tsv.md)
idxstats, MAPQ filter, coverage depth and thresholds.

### [7. Two-file Joins](docs/07-joins.md)
`NR==FNR`, left-join, `ARGIND` **[gawk]**, dynamic FS.

### [8. Reports & Aggregation](docs/08-reports.md)
Group-by, min/max/mean, pivot, sorted output **[gawk]**.

### [9. Advanced Patterns](docs/09-advanced.md)
Multi-dim arrays, `match()` capture **[gawk]**, `gensub()`, `$1=$1`, `delete`, `switch/case`.

### [10. getline](docs/10-getline.md)
Read from file/command, coprocess **[gawk]**, common pitfalls.

### [Testing AWK Scripts](docs/10-testing.md)
diff-based tests, bats integration, CI, POSIX compatibility testing.

### [11. GTF Annotation](docs/11-gtf.md)
gene_id/gene_name extraction, exon counts, transcript lengths, UTRs, Ensembl vs GENCODE.

### [12. RNA-seq Counts](docs/12-rnaseq.md)
HTSeq-count, featureCounts, CPM normalisation, merge N files into matrix, low-count filter.

### [13. Multi-sample Patterns](docs/13-multi-sample.md)
Cohort aggregation, outlier detection, QC summary, per-sample missingness.

### [AWK in Nextflow](docs/15-nextflow.md)
`script:` blocks, quoting rules, params, nf-core module patterns, debugging in `work/`.

### [Workflows](docs/16-workflows.md)
End-to-end pipelines: FASTQ QC, VCF+BED annotation, RNA-seq count pipeline, cohort VCF report.

### [Test Data](docs/data/README.md)
15 synthetic files (FASTQ, FASTA, VCF, BED, GTF, GFF3, depth, counts, metadata) — all inter-consistent.

---

## Test suite

```bash
# Install bats (Bash Automated Testing System)
brew install bats-core   # macOS
# or: sudo apt-get install bats   # Debian/Ubuntu

# Run all tests
bats docs/tests/
```

---

## Contributing

Open an issue or PR. New recipes should follow the format in any existing `.md` file.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — use freely, attribution appreciated.
