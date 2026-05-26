# Advanced Patterns

These recipes cover gawk-specific features and idiomatic AWK patterns that appear in
production bioinformatics scripts. Most require gawk unless marked POSIX.

---

## Multi-dimensional arrays

### Context

Standard AWK arrays are one-dimensional. True multi-dimensional arrays (gawk) use
`array[key1][key2]` syntax. POSIX AWK fakes it with `array[key1 SUBSEP key2]` — same
semantics, uglier syntax.

### Code

**gawk (recommended):**

```awk
# Store per-sample, per-chromosome coverage
{
    sample = $1; chrom = $2; depth = $3
    coverage[sample][chrom] += depth
    count[sample][chrom]++
}
END {
    for (s in coverage)
        for (c in coverage[s])
            printf "%s\t%s\t%.2f\n", s, c, coverage[s][c] / count[s][c]
}
```

**POSIX AWK (any awk):**

```awk
{
    key = $1 SUBSEP $2
    coverage[key] += $3
    count[key]++
}
END {
    for (key in coverage) {
        split(key, parts, SUBSEP)
        printf "%s\t%s\t%.2f\n", parts[1], parts[2], coverage[key] / count[key]
    }
}
```

### Explanation

- gawk's `array[a][b]` creates a true sub-array for each value of `a`.
  You can call `length(array[a])` to count the number of keys at the second level.
- POSIX `SUBSEP` (default `\034`, a non-printable character) creates a composite key.
  To check existence: `if ((a, b) in array)` is the POSIX equivalent of `if (b in array[a])`.
- Three-level nesting (`data[gene][typer][resolution]`) is common for HLA typing data,
  pipeline status tracking, and any multi-dimensional pivot structure.

---

## match() with capture groups **[gawk]**

### Context

Extract structured parts of a string using a single regex — no split, no substr arithmetic.
The 3-argument form of `match()` fills an array with the captured groups.

### Code

```awk
# Extract sample ID, UUID, and date from a filename:
# format: SAMPLENAME_xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx_20240115.vcf.gz

{
    if (match($0, /^(.+)_([0-9a-fA-F-]+)_([0-9]{8})\./, arr)) {
        sample_id = arr[1]
        uuid      = arr[2]
        date      = arr[3]
        print sample_id, uuid, date
    }
}
```

```awk
# Parse a Nextflow log line to extract process name, ID, and status
/INFO.*process >/ {
    if (match($0, /\[(..\/.{6})\] (\w+) process > (.+) \((.+)\)/, arr)) {
        hash    = arr[1]
        status  = arr[2]
        process = arr[3]
        nxf_id  = arr[4]
        counts[process][status]++
    }
}
```

### Explanation

- `match(string, /regex/, array)` — 3-argument form is gawk-only.
  - `arr[0]` = the full match.
  - `arr[1]`, `arr[2]`, ... = capture group 1, 2, ... (parenthesised groups).
- After a successful match, `RSTART` and `RLENGTH` are also set (as in 2-arg match).
- **Always `delete arr` before reusing the array in a loop** — stale captures from a
  longer previous match will pollute shorter subsequent matches.

```awk
# Delete before reuse
for (i = 1; i <= NF; i++) {
    delete arr
    if (match($i, /pattern(.+)/, arr))
        print arr[1]
}
```

---

## gensub() for in-place substitution **[gawk]**

### Context

`sub()` and `gsub()` modify in place and return the number of replacements.
`gensub()` returns the modified string without changing the original —
useful when you need both the original and modified value, or when working on a
field without modifying `$0`.

### Code

```awk
# Normalise chromosome names: remove "chr" prefix if present
{
    chrom = gensub(/^chr/, "", 1, $1)    # "g" for global, "1" for first-only
    print chrom, $2, $3
}
```

```awk
# Replace only the first occurrence
{ new = gensub(/foo/, "bar", 1, $2); print $1, new }

# Replace all occurrences
{ new = gensub(/foo/, "bar", "g", $2); print $1, new }

# Back-reference in replacement: swap two colon-separated fields
{ new = gensub(/^([^:]+):([^:]+)/, "\\2:\\1", 1, $3); print $1, new }
```

### Explanation

- Signature: `gensub(regexp, replacement, how [, target])`
  - `how`: `"g"` = global, `1` = first occurrence, `2` = second, etc.
  - `target`: defaults to `$0` if omitted.
- Back-references in `replacement`: `\\1`, `\\2` refer to capture groups in the regex.
- `gsub()` POSIX alternative: `tmp = $3; gsub(/^chr/, "", tmp)` — modifies `tmp` in place.

---

## Force field reassembly with $1=$1

### Context

Changing `OFS` does not automatically reformat `$0` — AWK only rebuilds `$0` from
fields when you modify a field. The idiom `$1=$1` (or `$1=$1""`) triggers that rebuild
without actually changing the content.

### Code

```awk
# Convert a CSV file to TSV
BEGIN { FS=","; OFS="\t" }
{ $1=$1; print }
```

```awk
# Trim leading/trailing spaces from all fields
{
    for (i=1; i<=NF; i++) gsub(/^ +| +$/, "", $i)
    $1=$1    # rebuild $0 with OFS
    print
}
```

### Explanation

- `$1=$1` assigns field 1 to itself — a no-op for the data, but it marks the record as
  "modified," causing AWK to reconstruct `$0` by joining all fields with `OFS`.
- After this, `print` (with no arguments, which prints `$0`) outputs the reformatted line.
- Alternative: `print $1, $2, ..., $NF` — explicit but tedious for many fields.

---

## delete array for memory management

### Context

gawk's `delete` removes an entry or an entire array from memory. In long-running scripts
processing millions of records, not deleting accumulators can exhaust memory.

### Code

```awk
# Free per-sample data after printing it
END {
    for (sample in data) {
        print_report(sample, data[sample])
        delete data[sample]    # free memory for this sample
    }
}
```

```awk
# Reset a capture array before reuse
{
    delete arr
    if (match($0, /pattern(.+)/, arr))
        process(arr[1])
}
```

```awk
# Delete an entire array at once (gawk)
END { delete counts }

# POSIX equivalent: loop and delete each key
END { for (k in counts) delete counts[k] }
```

### Explanation

- `delete array[key]` removes one entry.
- `delete array` (gawk) removes all entries and the array itself in one operation.
- Critical pattern: always `delete arr` before calling `match()` with a capture array
  in a loop. If the previous match captured 3 groups and the new one only captures 1,
  `arr[2]` and `arr[3]` will still hold stale values.

---

## split() to initialise lookup tables

### Context

A list of valid values is known at script-writing time (gene names, filter tags,
status codes). Load it into an associative array for O(1) membership tests.

### Code

```awk
BEGIN {
    classic_genes = "A,B,C,DRA,DRB1,DPA1,DPB1,DQA1,DQB1"
    split(classic_genes, tmp, ",")
    for (i = 1; i <= length(tmp); i++)
        is_classic[tmp[i]] = 1
}
{
    gene = $2
    if (gene in is_classic)
        print "CLASSIC:", $0
}
```

### Explanation

- `split(str, arr, sep)` fills `arr[1]`, `arr[2]`, ... and returns the number of parts.
- The subsequent `for` loop converts the indexed array into a lookup array (`is_classic`).
- `gene in is_classic` is O(1) and does not modify the array (unlike `is_classic[gene]`
  which creates an empty entry).
- Useful when the list is short and stable. For large or dynamic lists, load from a file
  with `NR==FNR`.

---

## switch/case **[gawk]**

### Context

Multiple exclusive conditions on the same variable. More readable than a chain of
`if/else if` when the cases are values of a single variable.

### Code

```awk
# Classify HLA allele resolution based on the number of fields
{
    n_fields = split($6, f, ":")
    switch (n_fields) {
        case 4:
            resolution = "4-field"
            allele_2f  = f[1] ":" f[2]
            allele_3f  = f[1] ":" f[2] ":" f[3]
            break
        case 3:
            resolution = "3-field"
            allele_2f  = f[1] ":" f[2]
            allele_3f  = $6
            break
        case 2:
            resolution = "2-field"
            allele_2f  = $6
            allele_3f  = $6
            break
        default:
            resolution = "unknown"
            allele_2f  = $6
            allele_3f  = $6
    }
    print $1, allele_2f, allele_3f, resolution
}
```

### Explanation

- `switch` is gawk-specific. POSIX equivalent: chain of `if/else if`.
- `break` is required to prevent fall-through (same as C).
- Cases can be strings, numbers, or regex: `case /^chr/: ...`

---

!!! note "See also: getline"
    The next natural step after mastering these patterns is `getline` — AWK's facility
    for reading from a file or command mid-script. It is essential for bidirectional
    pipelines, reading companion files without loading them into memory, and coprocesses.
    See [getline](10-getline.md).
