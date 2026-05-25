# Reports and Aggregation

AWK is a natural fit for computing summary statistics from tabular data: counts per
group, min/max/mean, pivot tables. No Python, no R, no temporary files.

---

## Group-by and count

### Context

How many variants per sample? How many reads per chromosome? Any "count by key" question.

### Code

```awk
# Count records per value of column 1
{ counts[$1]++ }
END {
    for (key in counts)
        print key, counts[key]
}
```

```bash
# Sorted output, most frequent first
awk '{c[$1]++} END{for(k in c) print c[k], k}' data.tsv | sort -rn
```

### Explanation

- `counts[$1]++` increments the counter for key `$1`. First access auto-initialises to 0.
- The `END` block iterates over all keys; order is arbitrary (hash order).
- Pipe to `sort` for deterministic output.

### Variants

```awk
# Count unique values in column 2, grouped by column 1
{ seen[$1][$2] = 1 }
END {
    for (group in seen)
        print group, length(seen[group])
}
```

---

## Compute min, max, mean per group

### Context

Per-sample coverage statistics, per-gene expression range, per-chromosome variant density.

### Code

```awk
# File: coverage.tsv  columns: sample_id  chrom  depth
NR > 1 {
    sample = $1
    depth  = $3 + 0

    count[sample]++
    sum[sample]   += depth

    if (!(sample in min_val) || depth < min_val[sample]) min_val[sample] = depth
    if (!(sample in max_val) || depth > max_val[sample]) max_val[sample] = depth
}
END {
    printf "%-15s %8s %8s %10s %8s\n", "SAMPLE", "MIN", "MAX", "MEAN", "N"
    for (s in count)
        printf "%-15s %8.1f %8.1f %10.2f %8d\n",
            s, min_val[s], max_val[s], sum[s]/count[s], count[s]
}
```

### Explanation

- `!(sample in min_val)` is true on the first record for that sample, so the initial
  value is set directly rather than comparing against an arbitrary sentinel.
- `depth + 0` forces numeric interpretation.

---

## Pivot: rows to columns

### Context

Input: long-format table with `sample`, `metric`, `value`.
Output: wide-format table with one column per metric.

```
# Input                    # Output
S1  coverage  45.2         SAMPLE  coverage  n_variants
S1  n_variants 12345       S1      45.2      12345
S2  coverage  32.1         S2      32.1       9876
S2  n_variants  9876
```

### Code

```awk
NR > 1 {
    data[$1][$2] = $3
    if (!($2 in metrics)) {
        metrics[$2] = ++n_metrics
    }
}
END {
    # Print header
    printf "SAMPLE"
    for (i = 1; i <= n_metrics; i++) printf "\t%s", metric_names[i]
    printf "\n"

    for (sample in data) {
        printf "%s", sample
        for (i = 1; i <= n_metrics; i++)
            printf "\t%s", (metric_names[i] in data[sample] ? data[sample][metric_names[i]] : "NA")
        printf "\n"
    }
}
# Build ordered metric list
NR > 1 {
    if (!($2 in metrics)) {
        n_metrics++
        metrics[$2] = n_metrics
        metric_names[n_metrics] = $2
    }
}
```

### Simpler version (fixed metrics known in advance)

```awk
NR > 1 { data[$1][$2] = $3 }
END {
    print "SAMPLE\tcoverage\tn_variants"
    for (s in data)
        print s "\t" data[s]["coverage"] "\t" data[s]["n_variants"]
}
```

---

## Sorted output with PROCINFO **[gawk]**

### Context

AWK's `for (key in array)` iterates in arbitrary (hash) order. In gawk, `PROCINFO["sorted_in"]`
controls the iteration order globally or per-array.

### Code

```awk
END {
    PROCINFO["sorted_in"] = "@ind_str_asc"    # alphabetical by key
    for (chrom in counts)
        print chrom, counts[chrom]
}
```

### Sort modes

| Value | Order |
|-------|-------|
| `@ind_str_asc`  | Keys, string ascending (alphabetical) |
| `@ind_str_desc` | Keys, string descending |
| `@ind_num_asc`  | Keys, numeric ascending |
| `@ind_num_desc` | Keys, numeric descending |
| `@val_str_asc`  | Values, string ascending |
| `@val_num_desc` | Values, numeric descending (largest first) |

### Explanation

- `PROCINFO["sorted_in"]` is set globally here; it applies to all subsequent `for` loops.
- Reset to `""` to restore hash-order iteration for specific arrays.
- Requires gawk. POSIX alternative: collect keys into an array, sort externally with
  `sort`, and process the sorted output in a second AWK pass.

---

## Print a matrix with ORS switching

### Context

Output a gene × typer × resolution matrix — rows are genes, columns are metrics.
Switching `ORS` lets you print comma-separated values on the same line, then a newline
at the end of each row.

### Code

```awk
# data[gene][metric] is populated earlier in the script
END {
    metrics_list = "cov_mean,cov_min,cov_max"
    split(metrics_list, metrics_arr, ",")

    # Header row
    ORS = ","
    print "GENE"
    for (i = 1; i < length(metrics_arr); i++) print metrics_arr[i]
    ORS = "\n"
    print metrics_arr[length(metrics_arr)]

    # Data rows
    for (gene in data) {
        ORS = ","
        print gene
        for (i = 1; i < length(metrics_arr); i++)
            print (metrics_arr[i] in data[gene] ? data[gene][metrics_arr[i]] : "NA")
        ORS = "\n"
        print (metrics_arr[length(metrics_arr)] in data[gene] ? data[gene][metrics_arr[length(metrics_arr)]] : "NA")
    }
}
```

### Explanation

- `ORS` is switched to `","` before printing column values, then back to `"\n"` for
  the last column so each row terminates properly.
- This avoids trailing commas without post-processing.
- The same effect can be achieved with `printf` which gives finer control without
  ORS manipulation.
