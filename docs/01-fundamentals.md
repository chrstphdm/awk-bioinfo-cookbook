# Fundamentals

---

## Field and record separators

### Context

AWK splits each input line into **fields** (`$1`, `$2`, ..., `$NF`) using the **field separator**
(`FS`). A complete line is a **record**; records are separated by the **record separator** (`RS`).
Getting these right is the first thing to do with any new file format.

### Code

```awk
# TSV input, tab-separated output
awk 'BEGIN { FS="\t"; OFS="\t" } { print $1, $3 }' file.tsv

# CSV input
awk 'BEGIN { FS="," } { print $2 }' file.csv

# Multiple possible separators (regex FS)
awk 'BEGIN { FS="[,;|]" } { print $1, $2 }' mixed.txt

# Set from the command line (no BEGIN block needed)
awk -F'\t' '{ print $1 }' file.tsv
```

### Explanation

| Variable | Default | Meaning |
|----------|---------|---------|
| `FS`     | `" "`   | Input field separator. A single space means "split on any whitespace and ignore leading/trailing whitespace." |
| `OFS`    | `" "`   | Output field separator, used when you print comma-separated fields. |
| `RS`     | `"\n"`  | Input record separator (one line = one record by default). |
| `ORS`    | `"\n"`  | Output record separator, printed after each `print` statement. |
| `NF`     | —       | Number of fields in the current record. |
| `NR`     | —       | Total number of records read so far (across all files). |
| `FNR`    | —       | Record number within the current file (resets at each new file). |

Setting `FS` in `BEGIN` applies from the very first record. Setting it mid-script only takes
effect on the **next** record.

### Variants

```awk
# Change OFS to reformat output on the fly
awk 'BEGIN { FS=","; OFS="\t" } { $1=$1; print }' file.csv
#   $1=$1 forces AWK to rebuild $0 with the new OFS — see the "Force reassembly" recipe.

# Multi-character FS
awk 'BEGIN { FS=" :: " } { print $2 }' annotations.txt
```

---

## BEGIN and END blocks

### Context

`BEGIN` runs once before any input is read. `END` runs once after all input is processed.
Use `BEGIN` to set separators, initialise variables, or validate arguments.
Use `END` to print summaries, flush output, or handle errors.

### Code

```awk
BEGIN {
    FS = "\t"
    OFS = "\t"
    total = 0
}
NR > 1 {          # skip header
    total += $3   # accumulate column 3
}
END {
    print "Total:", total
    print "Records processed:", NR - 1
}
```

### Explanation

- `BEGIN` and `END` are **pattern blocks** with no associated input record — `$0`, `NR`, `NF` are
  undefined (or zero) inside `BEGIN`.
- Multiple `BEGIN` and `END` blocks are legal and execute in order.
- `exit` in the middle of the script jumps straight to `END`.

### Variants

```awk
# Validate a mandatory -v argument
BEGIN {
    if (sample_id == "") {
        print "ERROR: -v sample_id=... is required" > "/dev/stderr"
        exit 1
    }
}
```

---

## Built-in variables

### Cheatsheet

| Variable   | Meaning |
|------------|---------|
| `NR`       | Total records read (across all files) |
| `FNR`      | Records read in the current file |
| `NF`       | Number of fields in the current record |
| `FS`       | Input field separator |
| `OFS`      | Output field separator |
| `RS`       | Input record separator |
| `ORS`      | Output record separator |
| `FILENAME` | Name of the current input file |
| `ARGC`     | Number of command-line arguments |
| `ARGV`     | Array of command-line arguments |
| `SUBSEP`   | Separator used in multi-dimensional array keys (`\034`) |

**gawk-specific:**

| Variable            | Meaning |
|---------------------|---------|
| `ARGIND`            | Index of the current file in ARGV (1-based) |
| `PROCINFO["sorted_in"]` | Controls array iteration order |
| `PROCINFO["version"]`   | gawk version string |

---

## Arithmetic and string operations

### Context

AWK types are dynamic — a variable is a number or a string depending on context.
Knowing when each applies avoids silent bugs.

### Code

```awk
# Arithmetic
{ sum += $3 }
END { print sum / NR }

# String concatenation (juxtaposition)
{ key = $1 "_" $2 }

# String comparison (lexicographic)
$3 > "PASS" { print }

# Numeric comparison
$5 > 30 { print }

# Force numeric context
{ val = $4 + 0 }   # ensures numeric comparison even if field looks like a string

# String length
{ if (length($2) > 20) print }

# Substring
{ prefix = substr($1, 1, 4) }

# Index / search
{ if (index($1, "chr") == 1) print }   # starts with "chr"

# sub() — replace first match
{ sub(/old/, "new", $2); print }

# gsub() — replace all matches
{ gsub(/,/, "\t"); print }   # CSV to TSV
```

---

## Ternary operator

### Context

AWK supports C-style ternary expressions. Useful for in-line conditional values without
a full `if/else` block.

### Code

```awk
# Print "PASS" or "FAIL" based on a quality score
{ label = ($5 >= 30) ? "PASS" : "FAIL"; print $1, label }

# Handle missing fields gracefully
{ allele2 = ($8 != "") ? $8 : $7 }

# Nested ternary (use sparingly — readability drops fast)
{ status = ($3 > 100) ? "HIGH" : ($3 > 10) ? "MED" : "LOW" }
```

### Explanation

```
(condition) ? value_if_true : value_if_false
```

The ternary is evaluated as an expression, so it works inside assignments, `print`, and
function arguments.

---

## printf for formatted output

### Context

`print` adds `ORS` at the end and `OFS` between arguments.
`printf` gives you full control — no automatic newline, C-style format strings.

### Code

```awk
# Column-aligned output
END {
    printf "%-20s %8s %6s\n", "SAMPLE", "READS", "BASES"
    for (s in counts)
        printf "%-20s %8d %6.2f\n", s, counts[s]["reads"], counts[s]["bases"] / 1e9
}

# Percentage
{ printf "%s\t%.2f%%\n", $1, ($2 / total) * 100 }

# Build a line incrementally, then terminate it
{
    ORS = "\t"
    for (i = 1; i <= NF; i++) printf "%s\t", $i
    printf "\n"
}
```

### Format specifiers

| Specifier | Meaning |
|-----------|---------|
| `%s`      | String |
| `%d`      | Integer |
| `%f`      | Float |
| `%e`      | Scientific notation |
| `%g`      | Shorter of `%f` or `%e` |
| `%-20s`   | Left-aligned, 20-char wide |
| `%08d`    | Zero-padded, 8-char wide |
| `%.4f`    | 4 decimal places |

---

## Writing to multiple output files

### Context

AWK can write to arbitrary files with `>` (truncate) or `>>` (append).
Pass file paths as `-v` variables so the script stays generic.

### Code

```awk
# Route records to different files based on a field value
{
    if ($4 == "PASS")
        print > passed_file
    else
        print > failed_file
}
```

```bash
awk -v passed_file=pass.tsv -v failed_file=fail.tsv -f route.awk input.tsv
```

### Explanation

- `>` truncates on the **first write** then appends. AWK keeps the file handle open for the
  duration of the run — no overhead from repeated opens.
- `>>` always appends (useful when the output file pre-exists).
- Use `close(filename)` to flush and close a file mid-script (e.g. if you need to re-read it).
- `/dev/stderr` is a valid target: `print "ERROR: ..." > "/dev/stderr"`.

### Variants

```awk
# One output file per chromosome
{ print > ("chr_" $1 ".bed") }

# Close files after writing to avoid hitting the OS open-file limit
{
    outfile = "sample_" $1 ".tsv"
    print > outfile
    close(outfile)
}
```

---

## Custom functions

### Context

AWK supports named functions for reusable logic. Arguments are passed by value (scalars)
or by reference (arrays). Functions can be defined anywhere in the script.

### Code

```awk
# Function: normalise a genomic coordinate string
function normalise_allele(allele,    fields, n) {
    # local variables declared after the last named param, separated by extra spaces
    n = split(allele, fields, ":")
    return fields[1] ":" fields[2]   # return 2-field resolution
}

{
    allele_2f = normalise_allele($6)
    print $1, allele_2f
}
```

### Explanation

- Function signature: `function name(param1, param2,    local1, local2)`
  Extra spaces before local variable names are a convention — AWK treats them as
  regular parameters that happen to be empty at call time.
- Arrays passed as arguments are passed **by reference** — modifications inside the
  function affect the caller's array.
- Scalar arguments are passed **by value**.
- No `return` is required; functions without `return` return an empty string.

---

## Error handling

### Context

AWK scripts used in pipelines need predictable error behaviour: print a clear message,
exit with a non-zero code, and route errors to stderr so they don't pollute stdout.

### Code

```awk
BEGIN {
    if (input_file == "") {
        print "ERROR: -v input_file=... is required" > "/dev/stderr"
        exit 1
    }
}

{
    if (NF < 5) {
        print "ERROR: line " NR " has only " NF " fields (expected >= 5)" > "/dev/stderr"
        exit 1
    }
}

END {
    # If exit was called with a non-zero code, END still runs.
    # Use a flag to avoid printing a "success" message on error.
    if (_exit_error) exit _exit_error
    print "Done. Processed " NR " records."
}
```

### Pattern: accumulated error flag

When you want to **collect all errors** before aborting (useful for input validation):

```awk
BEGIN {
    errors = ""
    if (sample_id == "")  errors = errors "ERROR: sample_id is required\n"
    if (cohort == "")     errors = errors "ERROR: cohort is required\n"
    if (errors != "") {
        printf "%s", errors > "/dev/stderr"
        exit 1
    }
}
```

### Explanation

- `exit N` in any block (including `BEGIN`) jumps to `END` and sets the process exit code to N.
- Inside `END`, calling `exit` again exits immediately without re-entering `END`.
- Writing to `"/dev/stderr"` keeps error messages out of stdout/pipes.
