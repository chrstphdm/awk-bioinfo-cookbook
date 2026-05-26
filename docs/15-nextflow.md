# AWK in Nextflow Pipelines

AWK is a natural fit for Nextflow `script:` blocks: it is always available on the
execution node, has no installation dependencies, and composes cleanly in shell
pipelines. It appears frequently in production Nextflow workflows as a glue step
between bioinformatics tools.

!!! warning "Examples not tested in a live pipeline"
    The AWK code in these examples has been validated independently on test data.
    However, the Nextflow process definitions (escaping, interpolation, output declarations)
    have **not** been executed in a `nextflow run`. If you find an escaping error or a
    Nextflow syntax issue, please [open an issue](https://github.com/chrstphdm/awk-bioinfo-cookbook/issues).

!!! note "Licence and attribution"
    All code examples in this chapter are original and written for this cookbook.
    Nextflow DSL2 syntax is described according to the [Nextflow documentation](https://www.nextflow.io/docs/latest/)
    (Apache 2.0 licence). nf-core module patterns are described based on the
    [nf-core module template](https://nf-co.re/docs/contributing/modules) (MIT licence);
    no code is copied from nf-core repositories.

---

## AWK in a `script:` block — the basics

### Context

A Nextflow `script:` block runs as a shell script on the worker node. AWK can be used
directly, just like in a terminal.

### Code

```nextflow
process FILTER_VARIANTS {
    input:
    path vcf

    output:
    path "filtered.vcf"

    script:
    """
    awk '/^#/ || (\$6 != "." && \$6+0 >= ${params.min_qual})' ${vcf} > filtered.vcf
    """
}
```

```nextflow
process COMPUTE_COVERAGE_STATS {
    input:
    path depth_file
    val  sample_id

    output:
    path "${sample_id}.coverage_stats.tsv"

    script:
    """
    awk '{ sum[\$1] += \$3; count[\$1]++ }
         END { for (c in sum) printf "%s\\t%s\\t%.2f\\n", "${sample_id}", c, sum[c]/count[c] }' \\
        ${depth_file} > ${sample_id}.coverage_stats.tsv
    """
}
```

---

## The quoting problem — and how to solve it

### Context

AWK uses single quotes in the shell: `awk '{ print $1 }'`. But Nextflow's default
`script:` blocks are single-quoted Groovy strings. This causes a conflict: the `'`
inside the AWK program ends the Groovy string prematurely.

### The three approaches

=== "Triple-quoted strings `\"\"\"` (recommended)"

    Use triple double-quotes for the script block. AWK single quotes work normally,
    but `$` variables require escaping: `\$1` for AWK fields, `${var}` for Nextflow
    variables.

    ```nextflow
    script:
    """
    awk '/^#/ || \$7 == "PASS"' ${vcf} > pass_only.vcf
    """
    ```

    - AWK single quotes: work without escaping ✓
    - AWK field variables: escape as `\$1`, `\$NF` ✓
    - Nextflow variables: `${vcf}`, `${params.min_qual}` ✓
    - **Recommended for AWK-heavy processes**

=== "Single-quoted strings `'...'` (avoid)"

    Single-quoted Groovy strings do not interpolate `${}`, and AWK's `'` terminates
    the string. You must escape every AWK single quote as `\'`, which quickly becomes
    unreadable.

    ```nextflow
    // Avoid: requires escaping AWK single quotes as \'
    script:
    'awk \'/^#/ || $7 == "PASS"\' ' + vcf + ' > pass_only.vcf'
    ```

=== "External `.awk` file (recommended for >5 lines)"

    Place the AWK program in a `.awk` file in the pipeline's `bin/` directory and
    call it with `-f`. This **eliminates all quoting issues** — no `\$` escaping,
    no triple-quoted strings:

    ```nextflow
    // bin/filter_pass.awk is committed to the pipeline repo
    script:
    """
    awk -v min_qual=${params.min_qual} \
        -f ${projectDir}/bin/filter_pass.awk \
        ${vcf} > pass_only.vcf
    """
    ```

    - Zero quoting conflicts — AWK code is plain AWK, Nextflow code is plain Nextflow
    - Clean separation of AWK logic and Nextflow plumbing
    - Easier to test the AWK script independently (`awk -f bin/filter_pass.awk test.vcf`)
    - Version-controlled alongside the pipeline
    - Reusable across multiple processes
    - **Recommended when the AWK script is longer than ~5 lines**

    See [Fundamentals → AWK script files](01-fundamentals.md#awk-script-files-awk) for
    naming conventions (`.awk` vs `.gawk`), shebang lines, and the `-f` flag.

---

## Passing Nextflow parameters to AWK

### Context

AWK's `-v` flag passes shell variables into the AWK namespace. This is the correct
way to pass Nextflow `params.*` values without shell injection risks.

### Code

```nextflow
process FILTER_BY_QUAL {
    input:  path vcf
    output: path "filtered.vcf"

    script:
    def min_qual = params.min_qual ?: 30     // Groovy default value
    """
    awk -v min_qual=${min_qual} \\
        '/^#/ || (\$6 != "." && \$6+0 >= min_qual)' \\
        ${vcf} > filtered.vcf
    """
}
```

```nextflow
process EXTRACT_SAMPLE_COUNTS {
    input:
    path   vcf
    val    sample_name

    output:
    tuple val(sample_name), path("${sample_name}.counts.tsv")

    script:
    """
    awk -v sample="${sample_name}" '
        /^#CHROM/ { for(i=10;i<=NF;i++) if(\$i==sample) col=i }
        !/^#/ && col { split(\$col,g,":"); if(g[1]!="./." && g[1]!="0/0") n++ }
        END { print sample, n+0 }
    ' ${vcf} > ${sample_name}.counts.tsv
    """
}
```

!!! warning "Do not interpolate user input directly into AWK programs"
    Never do: `awk '{ if ($1 == "${params.chrom}") print }' file`
    The `${params.chrom}` value is interpolated directly into the AWK source — if it
    contains `'` or other shell metacharacters, it can break the script or cause
    unexpected behaviour.

    **Correct:** use `-v`: `awk -v chrom="${params.chrom}" '$1 == chrom' file`

---

## AWK as a format conversion step

### Context

AWK is ideal for lightweight format transformations between tools — converting output
from one tool into the input format expected by the next.

### Code

```nextflow
// Convert featureCounts output to clean gene × sample matrix (strip metadata cols)
process CLEAN_FEATURECOUNTS {
    input:  path counts
    output: path "counts_matrix.tsv"

    script:
    """
    awk '/^##/{next}
         NR==2{ printf "gene_id"; for(i=7;i<=NF;i++) printf "\\t%s", \$i; print ""; next }
         { printf "%s", \$1; for(i=7;i<=NF;i++) printf "\\t%s", \$i; print "" }' \\
        ${counts} > counts_matrix.tsv
    """
}
```

```nextflow
// Aggregate samtools idxstats across samples
process SUMMARISE_MAPPING {
    input:  tuple val(sample), path(idxstats)
    output: path "${sample}.mapping.tsv"

    script:
    """
    awk -v s="${sample}" \\
        '\$1 != "*" { mapped += \$3 }
         END        { print s, mapped }' \\
        ${idxstats} > ${sample}.mapping.tsv
    """
}

// Gather all per-sample summaries into a cohort table
process MERGE_MAPPING_STATS {
    input:  path "stats_*.tsv"
    output: path "cohort_mapping.tsv"

    script:
    """
    echo -e "sample\\tmapped_reads" > cohort_mapping.tsv
    cat stats_*.tsv | sort -k2,2rn >> cohort_mapping.tsv
    """
}
```

---

## AWK in nf-core modules — structure and `templates/`

### Context

nf-core modules follow a standard structure. For short AWK commands (1–5 lines),
inline `script:` blocks are standard. For longer AWK programs, nf-core uses a
`templates/` directory inside the module — this is the **current recommended approach**
for non-trivial scripts.

All examples below are original and illustrative — not copied from nf-core repositories.

### Module with inline AWK (short scripts)

```
modules/local/vcf_qual_filter/
├── main.nf
├── meta.yml
└── environment.yml
```

```nextflow
// modules/local/vcf_qual_filter/main.nf
process VCF_QUAL_FILTER {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: vcf
    path "versions.yml"                    , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def min_qual = task.ext.args   ?: "30"
    """
    awk -v min_qual=${min_qual} \\
        '/^#/ || (\$6 != "." && \$6+0 >= min_qual)' \\
        ${vcf} > ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

### Module with `templates/` directory (complex scripts)

When the AWK logic is > 5 lines, extract it into a `templates/` file. The `template`
directive loads the script from a file relative to the module's directory:

```
modules/local/vcf_cohort_report/
├── main.nf
├── meta.yml
├── environment.yml
└── templates/
    └── cohort_report.sh
```

```nextflow
// modules/local/vcf_cohort_report/main.nf
process VCF_COHORT_REPORT {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.report.tsv"), emit: report
    path "versions.yml"                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'cohort_report.sh'
}
```

```bash
#!/usr/bin/env bash
# templates/cohort_report.sh
# Nextflow variables (${prefix}, ${vcf}) are interpolated before execution.

awk '
/^#CHROM/ {
    for (i=10; i<=NF; i++) samples[i] = \$i
    n_samples = NF - 9; next
}
/^#/ { next }
{
    n_variants++
    is_pass = (\$7 == "PASS")
    n_fmt = split(\$9, fmt, ":")
    gt_idx = 0
    for (i=1; i<=n_fmt; i++) if (fmt[i]=="GT") { gt_idx=i; break }
    for (s=10; s<=NF; s++) {
        split(\$s, g, ":")
        gt = gt_idx ? g[gt_idx] : "./."
        if (gt == "./.") { missing[s]++ }
        else { total[s]++; if (is_pass) pass_count[s]++ }
    }
}
END {
    for (s=10; s<=9+n_samples; s++)
        printf "%s\\t%d\\t%.1f%%\\t%.1f%%\\n",
            samples[s], total[s]+0,
            (total[s]>0 ? pass_count[s]/total[s]*100 : 0),
            (n_variants>0 ? (missing[s]+0)/n_variants*100 : 0)
}' "${vcf}" > "${prefix}.report.tsv"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: \$(awk --version 2>&1 | head -1)
END_VERSIONS
```

### How `template` and `moduleDir` work

| Variable / directive | What it does |
|---|---|
| `template 'script.sh'` | Loads `templates/script.sh` relative to the module's `main.nf` |
| `moduleDir` | Resolves to the directory containing the current module's `main.nf` |
| `${moduleDir}/environment.yml` | Standard nf-core pattern for conda environment |

- Inside a template file, Nextflow variables (`${prefix}`, `${vcf}`, `${task.cpus}`)
  are interpolated **before** the script runs — exactly like an inline `script:` block.
- AWK `$` field variables must still be escaped as `\$` inside templates (same as inline).
- The template file must have a shebang (`#!/usr/bin/env bash`).

### When to use each approach

| Approach | When | AWK suitability |
|---|---|---|
| Inline `script:` | 1–5 line AWK one-liners | Most nf-core modules do this |
| `templates/` directory | Complex scripts (> 5 lines) | Recommended — clean, testable, version-controlled alongside module |
| `${projectDir}/bin/` + `-f` | Pipeline-specific scripts, not for shared modules | Legacy — still works for pipeline-level scripts |

### Explanation

- `task.ext.prefix` and `task.ext.args` allow callers to customise the output prefix
  and AWK threshold without modifying the module.
- The `versions.yml` block is the nf-core convention for recording tool versions. For
  AWK, `awk --version` (gawk) or `awk -W version` (mawk) provides the version string.
- `label 'process_low'` marks this as a lightweight step that does not need large CPU
  or memory allocations.

---

## `exec:` vs `script:` — when to use each

| Block | Runs where | Use when |
|---|---|---|
| `script:` | Worker node (shell) | AWK, samtools, any CLI tool |
| `exec:` | Nextflow head process (Groovy) | Pure data manipulation, no external tools |

AWK must always be in `script:`. `exec:` is for Groovy string manipulation and data
routing — it does not have a shell environment.

```nextflow
// CORRECT: AWK in script:
process FILTER {
    script:
    """
    awk '!/^#/{print}' ${vcf}
    """
}

// WRONG: AWK in exec: — exec: is Groovy, not shell
process FILTER_BROKEN {
    exec:
    // This will fail — there is no awk here
    "awk '!/^#/{print}' ${vcf}"
}
```

---

## Debugging AWK inside Nextflow

### Context

When an AWK step fails or produces unexpected output inside a Nextflow run, the first
step is to isolate the AWK script from the pipeline.

### Approach

```bash
# 1. Find the work directory for the failing task
nextflow run main.nf ... 2>&1 | grep "work/"
# Or check: .nextflow.log

# 2. Navigate to the work directory
cd work/ab/cdef1234...

# 3. The .command.sh file contains the exact script Nextflow ran
cat .command.sh

# 4. Run it directly to reproduce the error
bash .command.sh

# 5. Isolate the AWK command and test it manually
awk -v min_qual=30 '/^#/ || ($6 != "." && $6+0 >= min_qual)' input.vcf | head
```

```bash
# Add debug output to stderr without breaking the pipeline
awk -v min_qual=30 '
    /^#/ { print; next }
    {
        passes = ($6 != "." && $6+0 >= min_qual)
        if (NR <= 5) printf "DEBUG: QUAL=%s passes=%d\n", $6, passes > "/dev/stderr"
        if (passes) print
    }' input.vcf > filtered.vcf
```

```bash
# Use nextflow -resume to rerun only the failed step after fixing
nextflow run main.nf ... -resume
```

!!! tip "Check which AWK is on the execution node"
    Container images may have `mawk` instead of `gawk`. If your recipe uses
    `[gawk]` features, add a guard:

    ```bash
    if ! awk 'BEGIN{exit (length(PROCINFO)==0)}' 2>/dev/null; then
        echo "ERROR: gawk required, found: $(awk --version 2>&1 | head -1)" >&2
        exit 1
    fi
    ```
