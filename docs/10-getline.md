# getline — Reading Outside the Main Loop

`getline` lets AWK read input from a file or command mid-script, outside the normal
one-record-at-a-time flow. It is essential for: looking up data from a companion file
during processing, capturing the output of an external command, and reading multiple
lines of a single logical record at once.

`getline` returns **1** on success, **0** on EOF, **-1** on error.
**Always check the return value.**

---

## Read from a file mid-script

### Context

You are processing a file record by record and need to read one piece of data from a
companion file during processing — but the companion file is too large to load into
memory upfront with `NR==FNR`.

### Code

```awk
# Read a single-value file once in BEGIN (e.g. a genome size)
# File genome.size contains one line: "2789820082"
BEGIN {
    if ((getline genome_size < "genome.size") <= 0) {
        print "ERROR: cannot read genome.size" > "/dev/stderr"
        exit 1
    }
    close("genome.size")
}
{
    # Use genome_size in calculations
    coverage = $3 / genome_size
    printf "%s\t%.4f\n", $1, coverage
}
```

```awk
# Read a sorted companion file line-by-line alongside the main file
# Both files must be sorted on the same key (chrom + pos)
# Main file: depth.tsv (chrom pos depth)
# Companion: regions.bed (chrom start end name)

BEGIN {
    companion = "regions.bed"
    # Prime the first line
    if ((getline cline < companion) > 0) {
        split(cline, cf, "\t")
    }
}
{
    chrom = $1; pos = $2
    # Advance companion until it catches up
    while (cf[1] < chrom || (cf[1] == chrom && cf[3]+0 < pos)) {
        if ((getline cline < companion) <= 0) break
        split(cline, cf, "\t")
    }
    in_region = (cf[1] == chrom && cf[2]+0 <= pos && pos <= cf[3]+0)
    print $0, (in_region ? cf[4] : ".")
}
END { close(companion) }
```

### Explanation

- `getline var < "file"` reads the next line from `file` into `var` **without touching `$0` or `NF`**.
- The file handle stays open between records — each call reads the next line sequentially.
  This is stateful: AWK remembers where it left off.
- `close("file")` is optional here but is a good habit; it releases the file descriptor.
- For random-access lookups, the `NR==FNR` join pattern (see [Two-file Joins](07-joins.md))
  is simpler. Use `getline` when both files are sorted on the same key and you want to
  stream them together without loading either into memory.

---

## Read from a command

### Context

You need the current date, the output of `wc -l`, or any other shell command result
inside an AWK script — without exiting the script.

### Code

```awk
# Capture today's date in BEGIN
BEGIN {
    "date +%Y-%m-%d" | getline today
    close("date +%Y-%m-%d")
    print "# Report generated:", today
}
{ print }
```

```awk
# Count lines in a file named by a field, add the count as a new column
{
    cmd = "wc -l < " $2
    cmd | getline line_count
    close(cmd)
    print $1, $2, line_count + 0
}
```

```awk
# Sort a set of values collected from the input, then print them in order
{ values[NR] = $1 }
END {
    # Build a sort command from the collected values
    for (i = 1; i <= NR; i++) printf "%s\n", values[i] | "sort -n"
    close("sort -n")
}
```

!!! warning "Always call `close()` after reading from a command"
    If you omit `close(cmd)`, the pipe stays open and subsequent calls to the same
    command string return the **next line of the same process's output** — not a fresh
    invocation. This is the single most common `getline` bug.

    ```awk
    # BUG: close() is missing — second call gets line 2 of the same wc process
    { "wc -l < " $2 | getline n; print $1, n }

    # CORRECT: close() ensures a fresh invocation per record
    { cmd = "wc -l < " $2; cmd | getline n; close(cmd); print $1, n }
    ```

### Explanation

- `cmd | getline var` executes `cmd` as a shell command, reads one line of its stdout
  into `var`. Successive calls without `close()` read successive lines of the same
  process — useful for multi-line output, but usually not what you want per record.
- The command string is the pipe handle. `close("date +%Y-%m-%d")` must match the
  exact string used to open the pipe.
- Avoid spawning a new process per input line for performance-sensitive scripts:
  prefer loading data into an array in `BEGIN` when the data is small.

---

## Read the next record early (FASTQ alternative)

### Context

An alternative to the `NR % 4` idiom for FASTQ: explicitly consume all four lines of
a read block when you encounter the header line. This makes the block structure
explicit in the code.

### Code

```awk
# Usage: awk -v min_len=50 -f fastq_getline.awk reads.fastq
NR % 4 == 1 {
    header = $0
    getline seq
    getline plus_line
    getline qual
    if (length(seq) >= min_len)
        printf "%s\n%s\n%s\n%s\n", header, seq, plus_line, qual
}
```

### Explanation

- After the three `getline` calls, `$0` holds the quality line (the last one read).
- When the main loop advances with `next` or at end of the rule block, it reads the
  **next** line after the quality line — which is the next `@` header.
- The `NR % 4 == 1` guard ensures this only fires on header lines.

!!! note "NR % 4 vs getline — which to prefer?"
    The `NR % 4` idiom is more robust: it does not depend on the record structure being
    exactly 4 lines, and it handles malformed FASTQ more gracefully (a missing quality
    line does not offset all subsequent reads). Use `getline` when you need all four
    variables (header, seq, plus, qual) simultaneously in one rule — it reads more
    naturally. See [FASTQ Recipes](02-fastq.md) for the standard `NR % 4` approach.

---

## Coprocess — bidirectional communication **[gawk]**

### Context

A coprocess is a child program that stays running for the lifetime of the AWK script.
You can send lines to its stdin and read responses from its stdout, one line at a time.
Useful when calling an external tool per record would be too slow to spawn as a new
process each time.

### Code

```awk
# Toy example: uppercase each sequence using tr as a coprocess [gawk]
BEGIN {
    coprocess = "tr a-z A-Z"
}
/^>/ { print; next }    # pass headers through
{
    print $0 |& coprocess          # send sequence to coprocess stdin
    coprocess |& getline result    # read response from coprocess stdout
    print result
}
END {
    close(coprocess)
}
```

```awk
# Practical example: look up gene names via a Python helper script [gawk]
# The helper reads gene IDs on stdin and writes "id\tname" on stdout
BEGIN {
    lookup = "python3 lookup_genes.py"
}
/^#/ { print; next }
{
    print $3 |& lookup                 # send gene ID
    lookup |& getline response         # read back "id\tname"
    split(response, r, "\t")
    print $1, $2, r[2], $4, $5        # replace col 3 with gene name
}
END { close(lookup) }
```

### Explanation

- `print data |& proc` sends `data` to the coprocess's stdin.
- `proc |& getline var` reads one line from the coprocess's stdout into `var`.
- **Deadlock risk:** if the coprocess buffers output (e.g. Python's default stdout
  buffering), it may never send data back until it exits. Use line-buffered mode
  (`python3 -u` or `sys.stdout.flush()`) in the child process.
- `close(coprocess)` sends EOF to the coprocess stdin and waits for it to exit.
- `|&` is gawk-specific. POSIX AWK has no coprocess facility.

---

## Common getline pitfalls

!!! warning "Five mistakes to avoid"

    **1. Not checking the return value**

    ```awk
    # BUG: if the file is empty, getline returns 0 and var is ""
    getline var < "config.txt"
    print "Config:", var    # silently prints "Config: "

    # CORRECT
    if ((getline var < "config.txt") <= 0) {
        print "ERROR: config.txt missing or empty" > "/dev/stderr"; exit 1
    }
    ```

    **2. Forgetting `close()` on commands**

    See the warning in [Read from a command](#read-from-a-command) above.

    **3. `getline` without a variable replaces `$0` and resets `NF`**

    ```awk
    # BUG: getline without var replaces $0 — the original line is gone
    { getline; print $1 }    # prints field 1 of the NEXT line, not the current

    # CORRECT: use a variable to keep $0 intact
    { getline next_line; print $1, next_line }
    ```

    **4. Relative paths depend on the working directory, not the script location**

    ```awk
    # This opens "genome.size" relative to wherever the script is run from
    getline size < "genome.size"
    # Pass the path as a variable instead:
    # awk -v ref_dir=/data/ref '{ getline size < (ref_dir "/genome.size") }' ...
    ```

    **5. Using getline inside a loop without a termination condition**

    ```awk
    # BUG: infinite loop if getline never returns 0
    while (getline line < "big.file") { process(line) }
    # This is actually correct — the while loop terminates on EOF (return 0)
    # The bug is omitting the > 0 check which conflates error (-1) with EOF (0):
    while ((getline line < "big.file") > 0) { process(line) }   # CORRECT
    ```
