# Test Data

Synthetic files used throughout the cookbook recipes. All files are self-consistent: they share the same sample IDs (`NA12878`, `NA19238`, `NA20585`), chromosomes (`chr1`, `chr2`, `chr3`), and gene names (`GENE1`, `GENE2`).

!!! note "Synthetic data — not from real samples"
    The sample identifiers (`NA12878`, `NA19238`, `NA20585`) are used as realistic-looking
    labels only. All sequences, variant calls, depth values, and counts are **entirely
    fabricated** and do not originate from the 1000 Genomes Project or any other study.
    These files are released under the same [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
    licence as the rest of the cookbook.

Clone the repo and run recipes from the `docs/data/` directory:

```bash
git clone https://github.com/chrstphdm/awk-bioinfo-cookbook
cd awk-bioinfo-cookbook/docs/data
awk 'NR%4==1' reads.fastq | wc -l   # → 20 reads
```

---

## Files

### `reads.fastq`

**Format:** FASTQ (4 lines per read)  
**Used in:** Chapter 02 (FASTQ recipes)  
**Content:** 20 reads with varying lengths (15–100 bp). Includes:
- 8 reads ≥ 50 bp (used in length-filter recipe)
- 3 duplicate read IDs: `DUP.read1` (×2) and `DUP.read2` (×2)
- 2 short reads < 30 bp: `SHORT.read1`, `SHORT.read2`
- Illumina-style headers with instrument/lane metadata after space

---

### `ids.txt`

**Format:** Plain text, one read ID per line (no `@`)  
**Used in:** Chapter 02 — extract reads by ID list  
**Content:** 5 read IDs to extract from `reads.fastq`

---

### `genome.fasta`

**Format:** FASTA, multi-line sequences (wrapped at 60 chars)  
**Used in:** Chapter 03 (FASTA recipes)  
**Content:** 5 sequences:
- `chr1` (500 bp) — contains an N-stretch at positions 121–130
- `chr2` (200 bp)
- `chr3` (50 bp) — very short, useful for length-filter
- `scaffold_001` (80 bp)
- `scaffold_002` (30 bp)

---

### `rename.tsv`

**Format:** 2-column TSV (old_name → new_name), no header  
**Used in:** Chapter 03 — rename FASTA headers from lookup table  
**Content:** Mapping for the two scaffold sequences

---

### `regions.bed`

**Format:** BED6 (chrom, start, end, name, score, strand), 0-based half-open  
**Used in:** Chapter 04 (BED/GFF recipes), Chapter 16 (workflows)  
**Content:** 15 features across `chr1`, `chr2`, `chr3`, and `chrUn_001`

---

### `annotation.gff3`

**Format:** GFF3, tab-separated, 1-based closed coordinates  
**Used in:** Chapter 04 (BED/GFF recipes)  
**Content:** Complete gene models for `GENE1` (chr1, +) and `GENE2` (chr2, -):
- gene → mRNA → exon + CDS features
- 5' and 3' UTR features
- GFF3-style attributes: `ID=...; Name=...; Parent=...`

---

### `annotation.gtf`

**Format:** GTF (Ensembl-style), tab-separated, 1-based closed coordinates  
**Used in:** Chapter 11 (GTF recipes)  
**Content:** Same gene models as `annotation.gff3` but in GTF format:
- gene, transcript, exon, CDS, UTR feature types
- Quoted attributes: `gene_id "GENE1"; transcript_id "GENE1.1";`
- `gene_biotype` attribute (Ensembl convention; GENCODE uses `gene_type`)

---

### `variants.vcf`

**Format:** VCF 4.2, 3 samples (NA12878, NA19238, NA20585)  
**Used in:** Chapter 05 (VCF recipes), Chapter 16 (workflows)  
**Content:** 20 variants across `chr1` and `chr2`:
- 8 SNPs, 4 INDELs (insertion and deletion), 2 multi-allelic sites
- QUAL values: range 14–92, 2 records with `QUAL=.`
- FILTER values: 12× PASS, 5× LowQual, 1× LowDepth, 1× `./.:`
- FORMAT: `GT:DP:GQ:AD` — includes `0/0`, `0/1`, `1/1`, `./.` genotypes
- INFO: `AF`, `DP`, `ANN` (SnpEff-style annotation for functional impact)

---

### `alignments.idxstats`

**Format:** `samtools idxstats` output (4 columns: ref, length, mapped, unmapped)  
**Used in:** Chapter 06 (SAM-derived TSV recipes)  
**Content:** 5 chromosomes + `*` unmapped catch-all

---

### `depth.tsv`

**Format:** `samtools depth` output (3 columns: chrom, pos, depth)  
**Used in:** Chapter 06 (coverage recipes)  
**Content:** 200 positions (chr1:1000–1099 and chr2:500–599). Mean depth ~30x. Three zero-depth positions at chr2:523–525 for testing coverage thresholds.

---

### `metadata.tsv`

**Format:** TSV with header (sample_id, batch, sex, population, QC_status)  
**Used in:** Chapter 07 (joins), Chapter 13 (multi-sample patterns)  
**Content:** 5 samples — 3 match samples in `variants.vcf` (NA12878, NA19238, NA20585), 2 extra samples with no variant data (NA12890, NA18507) to demonstrate left-join behaviour.

---

### `htseq_counts.tsv`

**Format:** HTSeq-count output (2 columns: gene_id, count), no header  
**Used in:** Chapter 12 (RNA-seq recipes)  
**Content:** 15 genes + 5 `__` summary lines at the bottom (`__no_feature`, `__ambiguous`, `__too_low_aQual`, `__not_aligned`, `__alignment_not_unique`)

---

### `featurecounts.tsv`

**Format:** featureCounts output (2-line header + data)  
**Used in:** Chapter 12 (RNA-seq recipes)  
**Content:** 15 genes × 3 samples (NA12878, NA19238, NA20585). Columns 1–6 are featureCounts metadata (Geneid, Chr, Start, End, Strand, Length); columns 7+ are counts per sample.

---

### `samples_list.txt`

**Format:** Plain text, one sample ID per line  
**Used in:** Chapter 13 (multi-sample patterns)  
**Content:** 5 sample IDs for shell loop examples

---

## Quick sanity checks

```bash
# FASTQ: should print 20 headers
awk 'NR%4==1' reads.fastq | wc -l

# FASTA: should print 5 sequence names
awk '/^>/{print $1}' genome.fasta

# BED: should print 15 features
awk '!/^#/' regions.bed | wc -l

# VCF: should print 20 variants (non-header lines)
awk '!/^#/' variants.vcf | wc -l

# HTSeq: should print 15 gene lines (skip __ lines)
awk '$1 !~ /^__/' htseq_counts.tsv | wc -l
```
