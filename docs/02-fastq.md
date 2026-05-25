# FASTQ Recipes

FASTQ is a 4-line format: `@header`, sequence, `+`, quality scores.
The core AWK trick is `NR % 4` to identify which line you are on within each read block.

```
NR % 4 == 1  →  header line    (@READID ...)
NR % 4 == 2  →  sequence
NR % 4 == 3  →  separator      (+)
NR % 4 == 0  →  quality scores
```

---

## Navigate a FASTQ block

### Context

You need to do different things on each of the 4 lines of a FASTQ record —
store the header, inspect the sequence, then act on both together at the quality line.

### Code

```awk
{
    if      (NR % 4 == 1) header = $0
    else if (NR % 4 == 2) seq    = $0
    else if (NR % 4 == 3) next           # skip the "+" line
    else {
        # NR % 4 == 0: quality line — we have all 4 fields now
        qual = $0
        # ... process the read here ...
        print header "\n" seq "\n+\n" qual
    }
}
```

### Explanation

- `NR % 4` gives 1, 2, 3, 0 for the successive lines of each read.
- `next` skips to the next record without executing the rest of the block —
  handy to silently discard the `+` line.
- At `NR % 4 == 0` (quality line) all four lines are available via variables.

### Variants

```awk
# Print only sequence lines (every 4th line starting at line 2)
NR % 4 == 2 { print }

# Print only headers
NR % 4 == 1 { print substr($0, 2) }   # strip the leading @
```

---

## Filter reads by minimum length

### Context

Short reads degrade alignment quality. Remove reads below a length threshold
before mapping.

### Code

```awk
# Usage: awk -v min_len=50 -f filter_length.awk reads.fastq

{
    if      (NR % 4 == 1) header = $0
    else if (NR % 4 == 2) seq    = $0
    else if (NR % 4 == 3) next
    else {
        if (length(seq) >= min_len)
            print header "\n" seq "\n+\n" $0
    }
}
```

### Explanation

- `length(seq)` returns the number of characters in the sequence string.
- The quality line (`$0` at `NR % 4 == 0`) is always the same length as the sequence —
  both get printed or both get discarded.
- `min_len` is passed via `-v min_len=50` on the command line.

### Variants

```awk
# Also keep a count of filtered vs kept reads
{
    if      (NR % 4 == 1) header = $0
    else if (NR % 4 == 2) seq    = $0
    else if (NR % 4 == 3) next
    else {
        total++
        if (length(seq) >= min_len) {
            kept++
            print header "\n" seq "\n+\n" $0
        }
    }
}
END {
    print "Kept:", kept, "/", total > "/dev/stderr"
}
```

---

## Deduplicate reads by ID

### Context

Duplicate reads — same read ID appearing more than once — can occur after merging FASTQ
files. Keep only the first occurrence.

### Code

```awk
{
    if      (NR % 4 == 1) {
        # Extract read ID (everything up to the first space)
        split($0, parts, " ")
        read_id = substr(parts[1], 2)    # strip leading @
        is_dup = (read_id in seen)
        seen[read_id]++
        header = $0
    }
    else if (NR % 4 == 2) seq  = $0
    else if (NR % 4 == 3) next
    else {
        if (!is_dup)
            print header "\n" seq "\n+\n" $0
    }
}
```

### Explanation

- `seen[read_id]++` uses the `!seen[key]++` deduplication pattern:
  the first time a key appears, `seen[key]` is 0 (falsy), so `!seen[key]` is true.
  After the increment it becomes 1 (truthy), so subsequent occurrences are flagged.
- `in seen` checks key existence without creating an entry — used here to set a flag
  for use on the quality line.
- IDs are extracted from the header by splitting on spaces — Illumina headers have
  additional metadata after the first space.

---

## Count reads and basic stats

### Context

Quick sanity check: how many reads, what is the total base count, average read length?

### Code

```awk
NR % 4 == 2 {
    reads++
    bases += length($0)
}
END {
    if (reads > 0)
        printf "Reads: %d\nTotal bases: %d\nAvg length: %.1f\n", reads, bases, bases/reads
}
```

### Explanation

- Only the sequence line is inspected — headers and quality lines are ignored.
- `bases += length($0)` accumulates total nucleotides.
- `printf` with `%.1f` prints one decimal place for the average.

---

## Extract reads by ID list

### Context

You have a list of read IDs (one per line) and want to extract exactly those reads
from a FASTQ file. Classic two-file join, adapted for FASTQ.

### Code

```awk
# File 1: ids.txt — one read ID per line (no @)
# File 2: reads.fastq

NR == FNR {
    wanted[$1] = 1
    next
}
{
    if      (FNR % 4 == 1) {
        split($0, parts, " ")
        read_id = substr(parts[1], 2)
        in_wanted = (read_id in wanted)
        header = $0
    }
    else if (FNR % 4 == 2) seq  = $0
    else if (FNR % 4 == 3) next
    else {
        if (in_wanted)
            print header "\n" seq "\n+\n" $0
    }
}
```

```bash
awk -f extract_ids.awk ids.txt reads.fastq
```

### Explanation

- `NR == FNR` is true only while reading the **first file** (ids.txt) — at that point
  `NR` and `FNR` are the same because no file boundary has been crossed yet.
  `next` skips to the next record so no FASTQ processing happens for the ID file.
- Once AWK moves to the second file (reads.fastq), `NR > FNR`, so the FASTQ processing
  block takes over.
- `FNR % 4` (not `NR % 4`) is used for the FASTQ block because `NR` continues counting
  from the end of the first file.

### Variants

```awk
# Report missing IDs at the end
NR == FNR {
    wanted[$1] = 0   # 0 = not found yet
    next
}
FNR % 4 == 1 {
    split($0, parts, " ")
    read_id = substr(parts[1], 2)
    if (read_id in wanted) wanted[read_id] = 1
}
END {
    for (id in wanted)
        if (wanted[id] == 0)
            print "MISSING:", id > "/dev/stderr"
}
```
