# FASTA Recipes

FASTA is a two-part format per sequence: a header line starting with `>`, followed by
one or more sequence lines. Multi-line FASTA (sequence wrapped at 60 or 80 chars) is
common in reference genomes; single-line FASTA is easier to process.

The key AWK pattern is: detect `>` to start a new record, accumulate sequence lines,
act when the next `>` arrives or at `END`.

---

## Print headers only

### Context

List all sequence names in a FASTA file — for a quick inventory, or to build an ID list.

### Code

```awk
/^>/ { print substr($0, 2) }
```

Or, keeping the `>`:

```awk
/^>/ { print }
```

### Explanation

- `/^>/` is a pattern that matches lines starting with `>`.
- `substr($0, 2)` strips the leading `>` character.
- The block is executed for every matching line — no state needed.

---

## Reformat multi-line FASTA to single-line

### Context

Most AWK processing assumes one sequence = one line. Convert a wrapped FASTA
(sequence split across multiple lines) to single-line format first.

### Code

```awk
/^>/ {
    if (seq != "")
        print header "\n" seq
    header = $0
    seq = ""
    next
}
{ seq = seq $0 }
END {
    if (seq != "") print header "\n" seq
}
```

### Explanation

- When a new `>` header is found, print the **previous** record (if any), then reset.
- Non-header lines are concatenated onto `seq` without separator — `seq = seq $0`.
- The `END` block prints the last record (there is no following `>` to trigger the flush).

### Variants

```awk
# One-liner version using RS and ORS trick (gawk / mawk)
awk 'BEGIN { RS=">"; OFS="\n" } NR > 1 {
    split($0, lines, "\n")
    header = ">" lines[1]
    seq = ""
    for (i=2; i<=length(lines); i++) seq = seq lines[i]
    print header "\n" seq
}' input.fasta
```

---

## Filter sequences by minimum length

### Context

Remove short contigs from an assembly or short sequences from a database before
running a search.

### Code

```awk
# Usage: awk -v min_len=500 -f filter_fasta.awk assembly.fasta

/^>/ {
    if (header != "" && length(seq) >= min_len)
        print header "\n" seq
    header = $0
    seq = ""
    next
}
{ seq = seq $0 }
END {
    if (header != "" && length(seq) >= min_len)
        print header "\n" seq
}
```

### Explanation

- The flush-on-next-header pattern (same as the single-line reformatter above) is
  extended with a `length(seq) >= min_len` guard.
- `min_len` defaults to 0 if not set with `-v` — all records pass through.

### Variants

```awk
# Keep a count of kept vs discarded sequences
/^>/ {
    if (header != "") {
        total++
        if (length(seq) >= min_len) { kept++; print header "\n" seq }
    }
    header = $0; seq = ""; next
}
{ seq = seq $0 }
END {
    if (header != "") {
        total++
        if (length(seq) >= min_len) { kept++; print header "\n" seq }
    }
    print "Kept " kept "/" total " sequences" > "/dev/stderr"
}
```

---

## Rename headers from a lookup table

### Context

Rename sequence IDs in a FASTA file using a two-column TSV: old name → new name.
Common after assembly, before database submission.

### Code

```awk
# File 1: rename.tsv — two columns: old_id  new_id
# File 2: sequences.fasta

NR == FNR {
    names[$1] = $2
    next
}
/^>/ {
    old_id = substr($0, 2)
    if (old_id in names)
        print ">" names[old_id]
    else
        print $0    # keep original if no mapping found
    next
}
{ print }
```

```bash
awk -f rename_fasta.awk rename.tsv sequences.fasta
```

### Explanation

- `NR == FNR` loads the lookup table from the first file.
- For each header line in the FASTA, the old ID (after stripping `>`) is looked up in
  the `names` array. Unknown IDs pass through unchanged.
- Non-header lines (`{ print }`) are printed as-is.

---

## Count sequences and total bases

### Context

How many sequences in this FASTA? What is the total assembly size? N50?

### Code

```awk
/^>/ {
    if (seq != "") {
        count++
        n = length(seq)
        total += n
        lens[count] = n
    }
    seq = ""
    next
}
{ seq = seq $0 }
END {
    if (seq != "") {
        count++
        n = length(seq)
        total += n
        lens[count] = n
    }
    # Sort lengths descending to compute N50
    n50 = compute_n50(lens, count, total)
    printf "Sequences:   %d\n", count
    printf "Total bases: %d\n", total
    printf "N50:         %d\n", n50
}

function compute_n50(lens, n, total,    i, cumul, sorted) {
    # Bubble sort descending (fine for < 10k sequences; use sort | awk for large assemblies)
    for (i = 1; i <= n; i++) sorted[i] = lens[i]
    for (i = 1; i <= n; i++)
        for (j = i+1; j <= n; j++)
            if (sorted[j] > sorted[i]) { tmp=sorted[i]; sorted[i]=sorted[j]; sorted[j]=tmp }
    cumul = 0
    for (i = 1; i <= n; i++) {
        cumul += sorted[i]
        if (cumul >= total / 2) return sorted[i]
    }
}
```

### Explanation

- Sequence lengths are stored in `lens[count]` for N50 computation.
- `compute_n50` sorts lengths descending and walks until the cumulative sum exceeds
  half the total — the length at that point is the N50.
- For very large assemblies (millions of sequences), the bubble sort is too slow.
  In practice: pipe through `sort -rn` and use a simpler AWK accumulator.

### Variants

```bash
# Faster N50 for large assemblies: get lengths first, then sort externally
awk '/^>/ { if (seq!="") print length(seq); seq=""; next } { seq=seq$0 } END { if(seq!="") print length(seq) }' assembly.fasta \
  | sort -rn \
  | awk -v total=$(grep -v ">" assembly.fasta | tr -d '\n' | wc -c) \
        '{ cumul+=$1; if(cumul >= total/2) { print $1; exit } }'
```
