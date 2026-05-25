# Two-file Joins

Joining two files on a shared key is one of the most common AWK tasks in bioinformatics:
annotate variants with gene names, filter a VCF by a sample list, enrich a TSV with
metadata. AWK handles this without writing to disk.

---

## NR==FNR join (inner join)

### Context

Match records in file 2 against a lookup table built from file 1.
Only print records that appear in both files.

### Code

```awk
# File 1: metadata.tsv  — sample_id  batch  sex  population
# File 2: results.tsv   — sample_id  coverage  n_variants

NR == FNR {
    meta[$1]["batch"]  = $2
    meta[$1]["sex"]    = $3
    meta[$1]["pop"]    = $4
    next
}
$1 in meta {
    print $1, meta[$1]["batch"], meta[$1]["sex"], meta[$1]["pop"], $2, $3
}
```

```bash
awk -F'\t' -f join.awk metadata.tsv results.tsv
```

### Explanation

- `NR == FNR` is the canonical two-file idiom. `NR` counts all records seen so far across
  all files; `FNR` counts records within the current file. They are equal **only while
  reading the first file**.
- `next` in the first-file block skips to the next record, so the second block never
  fires during file 1 processing.
- `$1 in meta` checks key existence without creating an entry (safer than `meta[$1] != ""`
  which creates an empty entry as a side effect).

### Common mistakes

```awk
# WRONG: this fires for every line in both files
{ if (NR == FNR) ... }

# RIGHT: use NR==FNR as the rule pattern
NR == FNR { ... }
```

---

## Left-join: keep all rows from file 1

### Context

Annotate every record from file 2 with metadata from file 1, even if no match exists.
Output a default value when the key is missing.

### Code

```awk
NR == FNR {
    annot[$1] = $2
    next
}
{
    val = ($1 in annot) ? annot[$1] : "NA"
    print $0, val
}
```

### Explanation

- The ternary `($1 in annot) ? annot[$1] : "NA"` avoids creating spurious empty entries
  in the `annot` array for unmatched keys.
- Every record from file 2 is printed — this is the left-join behaviour.

---

## Multi-file processing with ARGIND **[gawk]**

### Context

When you have more than two files, or when each file needs different processing logic,
`ARGIND` is cleaner than tracking file transitions manually.
`ARGIND` is the 1-based index of the current file in the argument list.

### Code

```awk
#!/usr/bin/gawk -f
# Usage: gawk -f multi_file.awk cohort.tsv delivered.tsv qc_flags.tsv

ARGIND == 1 {
    # Load cohort membership
    cohort[$1] = $2
    next
}
ARGIND == 2 {
    # Load delivered samples
    delivered[$1] = $0
    next
}
ARGIND == 3 {
    # Process QC flags, join with previous data
    sample = $1
    if (sample in cohort && sample in delivered)
        print sample, cohort[sample], $2
}
```

### Explanation

- Unlike `NR==FNR`, `ARGIND` works for any number of files and resets to the correct
  index at each file boundary.
- `FNR` still resets to 1 at each new file — combine `ARGIND` with `FNR != 1` to skip
  per-file headers.
- `ARGIND` is gawk-specific. POSIX equivalent: `FNR == 1 { ++file_id }`.

### POSIX equivalent (no gawk)

```awk
FNR == 1 { ++file_id }

file_id == 1 { cohort[$1] = $2; next }
file_id == 2 { delivered[$1] = $0; next }
file_id == 3 {
    if ($1 in cohort && $1 in delivered)
        print $1, cohort[$1], $2
}
```

---

## Dynamic FS per input file

### Context

Different files in the same AWK run use different separators (e.g. a CSV metadata file
and a space-delimited PLINK output). Set `FS` differently for each file.

### Code

```awk
# File 1: comma-separated   File 2: space-separated

FNR == 1 { ++file_id }

file_id == 1 {
    FS = ","
    # Note: FS change takes effect on the NEXT record.
    # For the very first line, set FS in BEGIN or use a workaround:
    if (FNR == 1) { FS = ","; $0 = $0 }   # re-split current line
    data[$1] = $2
    next
}
file_id == 2 {
    FS = " "
    if ($1 in data) print $1, data[$1], $3
}
```

### Practical note

Changing `FS` mid-record does not re-split the current `$0`. To re-split immediately,
reassign `$0 = $0` after changing `FS`, or set `FS` in `BEGIN` if it applies to the
first file. For clarity, prefer `ARGIND`-based processing (gawk) which makes per-file
logic explicit.
