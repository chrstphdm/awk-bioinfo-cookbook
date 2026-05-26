# GTF Annotation Recipes

GTF (Gene Transfer Format) is the annotation format consumed by RNA-seq aligners
(STAR, HISAT2) and quantification tools (featureCounts, StringTie). It follows the
same 9-column structure as GFF but uses a different attribute syntax and a stricter
feature hierarchy.

!!! tip "Test data"
    All recipes work with [`docs/data/annotation.gtf`](data/README.md):
    2 genes (GENE1 on chr1+, GENE2 on chr2-), each with transcript, exon, CDS, and UTR features.

---

## GTF format primer

### Columns

| # | Name | Content |
|---|---|---|
| 1 | seqname | Chromosome / scaffold name |
| 2 | source | Annotation source (e.g. `ensembl`, `havana`) |
| 3 | feature | Feature type: `gene`, `transcript`, `exon`, `CDS`, `UTR`, `start_codon`, `stop_codon` |
| 4 | start | 1-based start coordinate (closed) |
| 5 | end | 1-based end coordinate (closed) |
| 6 | score | `.` for most annotation files |
| 7 | strand | `+` or `-` |
| 8 | frame | `0`, `1`, or `2` for CDS; `.` otherwise |
| 9 | attributes | Semicolon-separated `key "value";` pairs |

### Attribute format: GTF vs GFF3

GTF attributes use **quoted values**:
```
gene_id "GENE1"; transcript_id "GENE1.1"; gene_name "GeneName1";
```

GFF3 attributes use **unquoted key=value**:
```
ID=gene:GENE1;Name=GENE1;Parent=transcript:GENE1.1
```

This difference matters for parsing — a GFF3 parser will not work on GTF attributes.

### Feature hierarchy

```
gene
└── transcript
    ├── exon        (one per exon; exon_number attribute)
    ├── CDS         (coding sequence, has frame field)
    ├── UTR         (Ensembl) or five_prime_utr / three_prime_utr (GENCODE)
    ├── start_codon
    └── stop_codon
```

`gene_id` is present on **every** row. `transcript_id` is present on transcript-level
and below features. `exon_number` is present on exon and CDS rows.

### Ensembl vs GENCODE variants

Both use GTF but differ in attribute names and chromosome naming:

| Attribute | Ensembl | GENCODE |
|---|---|---|
| Gene biotype | `gene_biotype "protein_coding"` | `gene_type "protein_coding"` |
| Chr prefix | No (`1`, `X`) | Yes (`chr1`, `chrX`) |
| Gene version | `gene_version "5"` | Not always present |
| `gene_name` | Usually present | Always present |

!!! warning "Always use `-F'\\t'` and `LC_ALL=C` with GTF files"
    GTF is **tab-separated**, but column 9 (attributes) contains spaces: `gene_id "GENE1"; gene_name "GeneName1";`.
    If you use AWK's default field separator (whitespace), `$9` will only contain `gene_id`
    instead of the entire attribute string. Always set `FS="\t"` explicitly.

    Additionally, some locales (notably UTF-8 on macOS) cause gawk regex to fail on
    attribute strings containing `;` and `"`. Set `LC_ALL=C` for reliable matching:

    ```bash
    LC_ALL=C awk -F'\t' '...' annotation.gtf
    # or: LC_ALL=C awk 'BEGIN{FS="\t"} ...' annotation.gtf
    ```
    See [Why AWK? → Implementations](00-why-awk.md#awk-implementations-which-one-do-you-have) for details.

---

## Extract gene records: gene_id and gene_name

### Context

The gene-level summary line contains all gene-level attributes. Extract gene
coordinates and identifiers for a gene list or BED file.

### Code

```awk
# [gawk] — uses 3-argument match() for clean capture
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "gene" {
    delete arr
    match($9, /gene_id "([^"]+)"/, arr);   gene_id   = arr[1]
    match($9, /gene_name "([^"]+)"/, arr); gene_name = (arr[1] != "" ? arr[1] : ".")
    printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, $4, $5, gene_id, gene_name, $7
}
```

**POSIX version — split on `;`, then parse each attribute:**

```awk
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "gene" {
    gene_id = gene_name = "."
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", attrs[i])                # trim leading/trailing spaces
        if (match(attrs[i], /^gene_id "([^"]+)"/, a))   gene_id   = a[1]
        if (match(attrs[i], /^gene_name "([^"]+)"/, a)) gene_name = a[1]
    }
    printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, $4, $5, gene_id, gene_name, $7
}
```

For plain POSIX AWK (no 3-arg match), use `index()` + `substr()` to extract the
quoted value:

```awk
BEGIN { FS = "\t" }

function get_attr(attr_str, key,    pos, val) {
    pos = index(attr_str, key " \"")
    if (pos == 0) return "."
    pos += length(key) + 2              # skip key + space + opening quote
    val = substr(attr_str, pos)
    sub(/".*/, "", val)                 # trim from closing quote onwards
    return val
}

/^#/ { next }
$3 == "gene" {
    gene_id   = get_attr($9, "gene_id")
    gene_name = get_attr($9, "gene_name")
    printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, $4, $5, gene_id, gene_name, $7
}
```

### Explanation

- `delete arr` before reusing `arr` in `match()` prevents stale captures from a longer
  previous match polluting shorter subsequent matches. See [Advanced Patterns](09-advanced.md#match-with-capture-groups-gawk).
- `gsub(/^ +/, "", attrs[i])` strips leading spaces after semicolon splits — GTF
  attributes are typically written as `gene_id "X"; transcript_id "Y";` where each
  attribute after the first is preceded by a space.

---

## Count exons per transcript

### Context

How many exons does each transcript have? A quick QC check before computing lengths.

### Code

```awk
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "exon" {
    match($9, /transcript_id "([^"]+)"/, arr)   # [gawk]
    if (arr[1] != "") exon_count[arr[1]]++
}
END {
    for (tx in exon_count)
        print tx, exon_count[tx]
}
```

```bash
# Sorted by exon count descending
awk -F'\t' '/^#/{next} $3=="exon"{ match($9,/transcript_id "([^"]+)"/,a); c[a[1]]++ }
            END{ for(t in c) print t, c[t] }' annotation.gtf \
  | sort -k2,2rn
```

---

## Compute transcript lengths

### Context

Transcript length = **sum of exon lengths**, not `gene_end - gene_start`.
The gene span includes introns; the spliced transcript does not.

### Code

```awk
# [gawk]
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "exon" {
    delete arr
    match($9, /transcript_id "([^"]+)"/, arr)
    tx = arr[1]
    # GTF is 1-based closed: length = end - start + 1
    tx_length[tx] += $5 - $4 + 1
    tx_gene[tx]    = tx_gene[tx] == "" ? get_gene_id($9) : tx_gene[tx]
}

function get_gene_id(attrs,    a) {
    delete a
    match(attrs, /gene_id "([^"]+)"/, a)
    return a[1]
}

END {
    print "transcript_id\tgene_id\tspliced_length"
    for (tx in tx_length)
        printf "%s\t%s\t%d\n", tx, tx_gene[tx], tx_length[tx]
}
```

!!! note "Exon length vs gene span"
    For GENE1 in the test data, the gene spans chr1:1001–2500 (1500 bp), but its
    three exons total 200+200+200 = 600 bp of spliced sequence. AWK on exon rows
    gives the correct spliced length; AWK on the gene row gives the genomic span.

### Variants

```awk
# Report both spliced length and genomic span for comparison
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, a)
    gene_span[a[1]] = $5 - $4 + 1
}
$3 == "exon" {
    match($9, /gene_id "([^"]+)"/, a)
    exon_total[a[1]] += $5 - $4 + 1
}
END {
    printf "%-12s %12s %14s %10s\n", "gene_id", "genomic_span", "spliced_length", "pct_exonic"
    for (g in gene_span)
        printf "%-12s %12d %14d %9.1f%%\n",
               g, gene_span[g], exon_total[g]+0,
               (exon_total[g]+0) / gene_span[g] * 100
}
```

---

## Build a gene → transcript → exon hierarchy

### Context

Multi-dimensional arrays allow you to represent the full gene model in memory.
Useful for reporting or reformatting annotation data.

### Code

```awk
# [gawk] — builds a 3-level hierarchy
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "exon" {
    delete arr
    match($9, /gene_id "([^"]+)"/, arr);       gid = arr[1]
    match($9, /transcript_id "([^"]+)"/, arr); tid = arr[1]
    match($9, /exon_number "([^"]+)"/, arr);   en  = arr[1] + 0
    exons[gid][tid][en] = $4 "\t" $5 "\t" $7    # store start, end, strand
}
END {
    print "gene_id\ttranscript_id\texon_number\tstart\tend\tstrand"
    for (gid in exons)
        for (tid in exons[gid])
            for (en in exons[gid][tid])
                printf "%s\t%s\t%d\t%s\n", gid, tid, en, exons[gid][tid][en]
}
```

### Explanation

- `exons[gid][tid][en]` is a 3-level nested array (gawk only). In POSIX AWK, use
  `SUBSEP`: `exons[gid SUBSEP tid SUBSEP en]`.
- `delete arr` before each `match()` call prevents stale group captures.
- `en + 0` converts `exon_number` to an integer so iteration is in numeric order
  (in gawk, you can use `PROCINFO["sorted_in"] = "@ind_num_asc"` in the `for` loop).

---

## UTR positions (5' and 3', strand-aware)

### Context

5' and 3' UTRs are defined relative to the transcript's direction of transcription,
not absolute genomic position. On a `+` strand gene, the 5' UTR is at lower
coordinates; on a `-` strand gene, it is at higher coordinates.

### Code

```awk
# [gawk] — handles both Ensembl (UTR) and GENCODE (five_prime_utr/three_prime_utr)
BEGIN { FS = "\t" }
/^#/ { next }
$3 ~ /^(UTR|five_prime_utr|three_prime_utr)$/ {
    delete arr
    match($9, /transcript_id "([^"]+)"/, arr); tid = arr[1]
    strand = $7

    if ($3 == "five_prime_utr" || $3 == "three_prime_utr") {
        utr_type = $3                          # GENCODE: type is explicit
    } else {
        # Ensembl: infer from strand and relative position
        # (requires knowing the CDS boundaries — simplified heuristic here)
        utr_type = (strand == "+") ? "five_prime_utr" : "three_prime_utr"
    }
    print tid, $1, $4, $5, strand, utr_type
}
```

!!! note "Ensembl vs GENCODE UTR annotation"
    Ensembl uses the generic `UTR` feature type for all UTRs, requiring strand and CDS
    boundary information to distinguish 5' from 3'. GENCODE explicitly labels them
    `five_prime_utr` and `three_prime_utr`. When writing portable scripts, check which
    variant you have before parsing UTR records.

---

## Detect GTF format variant

### Context

A quick heuristic to determine whether a GTF uses Ensembl or GENCODE conventions —
useful at the start of a pipeline to set the right attribute name (`gene_biotype` vs
`gene_type`).

### Code

```awk
BEGIN { FS = "\t" }
/^#/ { next }
FNR <= 5 {
    # Chromosome naming: GENCODE uses "chr" prefix
    if ($1 ~ /^chr/) {
        format = "GENCODE"
    } else {
        format = "Ensembl"
    }
    # Attribute style: check for gene_type (GENCODE) vs gene_biotype (Ensembl)
    if ($9 ~ /gene_type/) {
        attr_biotype = "gene_type"
    } else {
        attr_biotype = "gene_biotype"
    }
    print "Detected:", format, "| biotype attribute:", attr_biotype > "/dev/stderr"
    exit
}
```

```bash
# One-liner version
awk -F'\t' '!/^#/ && FNR<=3 { print ($1~/^chr/ ? "GENCODE" : "Ensembl"); exit }' annotation.gtf
```

### Variants

```awk
# Report biotype distribution (works for both Ensembl and GENCODE)
BEGIN { FS = "\t" }
/^#/ { next }
$3 == "gene" {
    # Try gene_type first (GENCODE), then gene_biotype (Ensembl)
    delete arr
    if (match($9, /gene_type "([^"]+)"/, arr) || match($9, /gene_biotype "([^"]+)"/, arr))
        biotypes[arr[1]]++
}
END {
    for (bt in biotypes) print biotypes[bt], bt
}
```
