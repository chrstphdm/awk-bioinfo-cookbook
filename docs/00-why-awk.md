# Why AWK?

If you already use Python or R for bioinformatics, here is a concise answer to
"should I bother learning AWK?" — followed by an honest look at where it stops
being the right tool.

---

## The 30-second argument

AWK is available on every Unix system without installation, starts in milliseconds,
processes files in a single pass with constant memory, and composes naturally in
shell pipelines. For filtering, reformatting, and summarising tab-separated data —
which is most of what bioinformatics involves — it is often the fastest tool to both
write and run.

```bash
# Filter a VCF by QUAL ≥ 30, keeping headers — one line, no imports, no venv
awk '/^#/ || ($6 != "." && $6+0 >= 30)' variants.vcf > filtered.vcf
```

The equivalent in Python requires `import sys`, opening the file, a loop, a type
conversion, and a print statement — roughly 10 lines for 1 line of AWK.

---

## AWK implementations: which one do you have?

Not all `awk` commands are equal. The cookbook marks gawk-specific recipes with
**[gawk]** — here is what each implementation supports:

| Implementation | Where you find it | Strengths | Limitations |
|---|---|---|---|
| **gawk** (GNU AWK) | Linux (most distros), `brew install gawk` on macOS | Full feature set: 3-arg `match()`, `gensub()`, `ARGIND`, `PROCINFO`, coprocesses, `FPAT` | Slightly slower than mawk on simple one-liners |
| **mawk** | Ubuntu/Debian default `awk`, many HPC clusters | ~2–3× faster than gawk on large files with simple patterns | No 3-arg `match()`, no coprocess, limited regex |
| **nawk / BWK awk** | macOS default `awk`, Solaris | The "One True AWK" — reference implementation | Minimal features, rarely updated |
| **bioawk** | Optional install (`brew install bioawk`) | Format-aware: `-c bed`, `-c sam`, `-c fastx`; built-in `revcomp()`, `gc()`, `trimq()` | Fork of nawk, not actively maintained; less portable |
| **POSIX awk** | Guaranteed on any POSIX system | Maximum portability | Subset of gawk; no extended features |

**Detect your implementation:**

```bash
awk --version 2>&1 | head -1
# GNU Awk 5.x.x → you have gawk
# mawk 1.x.x → you have mawk
# (no output or error) → likely BSD/nawk

# gawk-specific check:
awk 'BEGIN { print PROCINFO["version"] }'
# prints version string if gawk, empty if not
```

**Practical advice:**

- On most Linux HPC clusters, `awk` is gawk. You can use `[gawk]` recipes directly.
- On macOS, the default `awk` is BSD awk (nawk). Install gawk: `brew install gawk`,
  then call it as `gawk`.
- **Set `LC_ALL=C` for reliable regex matching.** Some locales (notably UTF-8 on macOS)
  cause gawk regex to fail on strings containing special characters like `;` and `|`,
  which are common in VCF INFO and GFF/GTF attribute fields. Setting `LC_ALL=C` forces
  byte-level matching and is also significantly faster on large files:

```bash
# Recommended: prepend LC_ALL=C for any AWK processing of bioinformatics files
LC_ALL=C awk '...' variants.vcf

# Or export it once in your session / .bashrc
export LC_ALL=C
```

- If you need POSIX portability (e.g. scripts shared across heterogeneous clusters),
  avoid `[gawk]` features or add a runtime check:

```bash
if awk 'BEGIN{exit (length(PROCINFO)==0)}'; then
    AWK=awk
else
    AWK=gawk
fi
```

---

## Decision matrix: AWK vs Python vs R

| Situation | Best tool | Why |
|---|---|---|
| Filter/reformat a large TSV in a pipeline | **AWK** | Zero startup, O(1) memory, streams |
| Quick column arithmetic on plain TSV | **AWK** | Faster to write than pandas |
| Join two TSV files by a shared key | **AWK** | `NR==FNR` pattern, no library needed |
| One-liner in a Nextflow `script:` block | **AWK** | Shell-native, no dependency |
| Glue step between two bioinformatics tools | **AWK** | Composes in pipes |
| Parse JSON, XML, or YAML | **Python** | AWK has no structured format support |
| Complex multi-step logic across many files | **Python** | Functions, classes, libraries |
| HTTP requests, database queries | **Python** | AWK has no networking |
| Differential expression (DESeq2, edgeR) | **R** | Bioconductor ecosystem, no AWK equivalent |
| Statistical modelling, PCA, clustering | **R / Python** | scipy, numpy, Bioconductor |
| Visualisation (ggplot2, seaborn) | **R / Python** | AWK produces text, not plots |
| Interactive data exploration | **R / Python** | RStudio, Jupyter — AWK is batch-only |

**The key insight:** AWK and R are not competing for the same tasks. AWK sits earlier
in the pipeline — it filters, cleans, and reshapes data into a form that R or Python
can consume efficiently. They are collaborators, not alternatives.

---

## Side-by-side comparison

Three real tasks, compared honestly.

### Task 1: filter a VCF by QUAL ≥ 30, keep headers

=== "AWK"

    ```bash
    awk '/^#/ || ($6 != "." && $6+0 >= 30)' variants.vcf
    ```
    One line. Works on a 100 GB file with no memory overhead.

=== "Python"

    ```python
    import sys
    min_qual = 30
    with open("variants.vcf") as f:
        for line in f:
            if line.startswith("#"):
                sys.stdout.write(line)
            else:
                fields = line.split("\t")
                if fields[5] != "." and float(fields[5]) >= min_qual:
                    sys.stdout.write(line)
    ```
    10 lines. Equally memory-efficient (reads line by line), but requires a script file
    or an awkward `python3 -c` invocation.

=== "R"

    ```r
    # R is not the right tool here — it would load the full VCF into memory.
    # Use VariantAnnotation for VCF in R; do not parse VCF with read.table().
    ```

**Winner: AWK** — fastest to write, identical performance.

---

### Task 2: mean coverage depth per chromosome from `samtools depth`

=== "AWK"

    ```bash
    awk '{ sum[$1] += $3; count[$1]++ }
         END { for (c in sum) printf "%s\t%.2f\n", c, sum[c]/count[c] }' depth.tsv
    ```
    Streams the file once, O(1) memory per chromosome.

=== "Python"

    ```python
    import pandas as pd
    df = pd.read_csv("depth.tsv", sep="\t", header=None, names=["chrom","pos","depth"])
    print(df.groupby("chrom")["depth"].mean().reset_index().to_string(index=False))
    ```
    4 lines including pandas import. Loads the entire file into memory — problematic
    for whole-genome depth files (can be 10–50 GB).

=== "R"

    ```r
    depth <- read.table("depth.tsv", col.names=c("chrom","pos","depth"))
    aggregate(depth ~ chrom, data=depth, FUN=mean)
    ```
    Also loads the full file. Fine for small regions, not for whole-genome.

**Winner: AWK** — streaming is the decisive advantage for large depth files.

---

### Task 3: PCA on a genotype matrix

=== "AWK"

    ```awk
    # AWK cannot do PCA — it has no matrix algebra or eigendecomposition.
    # Use AWK to prepare the input matrix, then hand off to R or Python.
    awk 'NR>1 { for (i=7; i<=NF; i++) printf "%s%s", $i, (i<NF?"\t":"\n") }' \
        featurecounts.tsv > counts_only.tsv
    ```

=== "R"

    ```r
    counts <- read.table("counts_only.tsv", header=FALSE)
    pca <- prcomp(t(counts), scale.=TRUE)
    biplot(pca)
    ```

=== "Python"

    ```python
    import pandas as pd
    from sklearn.decomposition import PCA
    counts = pd.read_csv("counts_only.tsv", sep="\t", header=None)
    pca = PCA(n_components=2).fit_transform(counts.T)
    ```

**Winner: R / Python** — AWK has no linear algebra. Use AWK to extract the matrix,
then hand it to the right tool.

---

## AWK and columnar formats (Parquet, Arrow)

AWK reads plain text. It cannot open Parquet or Arrow binary files directly.

This matters because large bioinformatics datasets are increasingly stored as Parquet:
gnomAD (v4+), UK Biobank, and many public cohort studies use it because Parquet
compresses data in internal row groups rather than compressing the file as a whole —
making it both dense *and* splittable for distributed computing.

**The hybrid pipeline:** convert to TSV on the fly, pipe into AWK:

```bash
# Using DuckDB (fast, free, single binary)
duckdb -c "COPY (SELECT chrom, pos, ref, alt, af FROM 'gnomad.parquet' WHERE af > 0.01)
           TO '/dev/stdout' (FORMAT CSV, DELIMITER '\t', HEADER false)" \
  | awk '$4 != "." && length($3)==1 && length($4)==1'   # keep SNPs only
```

```bash
# Using parquet-tools + Python one-liner
parquet-tools cat --format jsonl gnomad.parquet \
  | python3 -c "import sys,json; [print('\t'.join(str(r[k]) for k in ['chrom','pos','ref','alt'])) for r in (json.loads(l) for l in sys.stdin)]" \
  | awk '...'
```

**When this is worth doing:** when you need a quick filter or aggregation on a Parquet
dataset and do not want to write a Python/R script. DuckDB handles the format, AWK
handles the logic. For repeated or complex queries, use DuckDB or Python directly.

---

## When AWK becomes the wrong tool

Be honest about the boundary:

| Problem | Why AWK struggles | Better choice |
|---|---|---|
| CSV with quoted fields (e.g. `"field,with,commas"`) | `FS=","` breaks on embedded commas; `FPAT` (gawk) helps but is fragile | Python `csv` module, `csvkit` |
| JSON / XML / YAML input | No parser — you'd need fragile regex hacks | `jq`, Python |
| HTTP requests, API calls | No networking | Python `requests`, `curl` |
| Parallel processing across many files | AWK is single-threaded | GNU Parallel + AWK, Python multiprocessing |
| Interactive debugging | No REPL, no debugger | Python (pdb), R (RStudio) |
| Complex multi-file logic with state | Quickly becomes unreadable | Python |
| Files > available RAM with group-by on high-cardinality keys | Hash table exhaustion | Sort first, then AWK; or Python |

!!! tip "Installation"
    === "macOS"
        ```bash
        brew install gawk
        # Then use 'gawk' explicitly, or alias: alias awk=gawk
        ```
    === "Linux (Debian/Ubuntu)"
        ```bash
        sudo apt-get install gawk
        ```
    === "HPC cluster"
        ```bash
        module load gawk   # or: which awk && awk --version
        # On most HPC systems, awk is already gawk
        ```
