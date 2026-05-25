# AWK Bioinformatics Cookbook

Practical AWK recipes for bioinformatics — organised by problem, not by AWK feature.

Each recipe answers a real question you'd type into a search engine:
*"How do I filter FASTQ reads by length?"*,
*"How do I extract a VCF field?"*,
*"How do I join two TSV files by sample ID?"*

For every recipe you get: context, code, line-by-line explanation, and variants.

---

## POSIX AWK vs gawk

Most recipes work with any POSIX-compliant AWK (`awk`, `mawk`, `nawk`).
Recipes marked **[gawk]** require [GNU AWK](https://www.gnu.org/software/gawk/) and use
features like `ARGIND`, the 3-argument `match()`, `gensub()`, or `PROCINFO`.

Check which AWK you have:

```bash
awk --version          # GNU Awk x.y.z  → gawk
awk 'BEGIN{print ARGV[0]}' /dev/null   # mawk, nawk, etc.
```

On most Linux HPC clusters, `awk` is gawk. On macOS, the default `awk` is BSD AWK —
install gawk with `brew install gawk` and call it as `gawk`.

---

## Recipes

### 1. Fundamentals
- [Field and record separators](docs/01-fundamentals.md#field-and-record-separators)
- [BEGIN and END blocks](docs/01-fundamentals.md#begin-and-end-blocks)
- [Built-in variables cheatsheet](docs/01-fundamentals.md#built-in-variables)
- [Arithmetic and string operations](docs/01-fundamentals.md#arithmetic-and-string-operations)
- [Ternary operator](docs/01-fundamentals.md#ternary-operator)
- [printf for formatted output](docs/01-fundamentals.md#printf-for-formatted-output)
- [Writing to multiple output files](docs/01-fundamentals.md#writing-to-multiple-output-files)
- [Custom functions](docs/01-fundamentals.md#custom-functions)
- [Error handling: flags and stderr](docs/01-fundamentals.md#error-handling)

### 2. FASTQ
- [Navigate a FASTQ block with NR%4](docs/02-fastq.md#navigate-a-fastq-block)
- [Filter reads by minimum length](docs/02-fastq.md#filter-reads-by-minimum-length)
- [Deduplicate reads by ID](docs/02-fastq.md#deduplicate-reads-by-id)
- [Count reads and basic stats](docs/02-fastq.md#count-reads-and-basic-stats)
- [Extract a subset of reads by ID list](docs/02-fastq.md#extract-reads-by-id-list)

### 3. FASTA
- [Print sequence headers only](docs/03-fasta.md#print-headers-only)
- [Filter sequences by minimum length](docs/03-fasta.md#filter-by-minimum-length)
- [Reformat multi-line FASTA to single-line](docs/03-fasta.md#reformat-to-single-line)
- [Rename headers from a lookup table](docs/03-fasta.md#rename-headers-from-lookup-table)
- [Count sequences and total bases](docs/03-fasta.md#count-sequences-and-total-bases)

### 4. BED / GFF
- [Filter BED features by chromosome](docs/04-bed-gff.md#filter-by-chromosome)
- [Filter GFF features by type](docs/04-bed-gff.md#filter-gff-by-type)
- [Reformat coordinates (0-based to 1-based)](docs/04-bed-gff.md#reformat-coordinates)
- [Compute feature sizes from BED](docs/04-bed-gff.md#compute-feature-sizes)
- [Extract GFF attribute field](docs/04-bed-gff.md#extract-gff-attribute)

### 5. VCF
- [Skip header lines](docs/05-vcf.md#skip-header-lines)
- [Filter variants by QUAL score](docs/05-vcf.md#filter-by-qual)
- [Filter by FILTER column value](docs/05-vcf.md#filter-by-filter-column)
- [Extract a field from the INFO column](docs/05-vcf.md#extract-info-field)
- [Count variants per chromosome](docs/05-vcf.md#count-per-chromosome)

### 6. SAM-derived TSV
- [Compute per-sample read counts from flagstat-style TSV](docs/06-sam-tsv.md#per-sample-read-counts)
- [Filter by mapping quality](docs/06-sam-tsv.md#filter-by-mapping-quality)
- [Summarise coverage per region](docs/06-sam-tsv.md#coverage-per-region)

### 7. Two-file Joins
- [Join two TSV files on a shared key (NR==FNR)](docs/07-joins.md#nrfnr-join)
- [Left-join: keep all rows from file 1](docs/07-joins.md#left-join)
- [Multi-file processing with ARGIND](docs/07-joins.md#argind-multi-file) **[gawk]**
- [Dynamic FS per input file](docs/07-joins.md#dynamic-fs-per-file)

### 8. Reports and Aggregation
- [Group-by and count](docs/08-reports.md#group-by-count)
- [Compute min, max, mean per group](docs/08-reports.md#min-max-mean)
- [Pivot: rows to columns](docs/08-reports.md#pivot)
- [Sorted output with PROCINFO](docs/08-reports.md#sorted-output) **[gawk]**
- [Print a matrix with ORS switching](docs/08-reports.md#matrix-output)

### 9. Advanced Patterns
- [Multi-dimensional arrays](docs/09-advanced.md#multi-dimensional-arrays)
- [match() with capture groups](docs/09-advanced.md#match-capture-groups) **[gawk]**
- [gensub() for in-place substitution](docs/09-advanced.md#gensub) **[gawk]**
- [Force field reassembly with $1=$1](docs/09-advanced.md#force-reassembly)
- [delete array for memory management](docs/09-advanced.md#delete-array)
- [split() to initialise lookup tables](docs/09-advanced.md#split-lookup-table)
- [switch/case](docs/09-advanced.md#switch-case) **[gawk]**

---

## Contributing

Open an issue or PR. New recipes should follow the format in any existing `.md` file.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — use freely, attribution appreciated.
