# Contributing

Thanks for considering a contribution. This guide is meant to stay focused, practical,
and tested — here is how to keep it that way.

## Recipe format

Every recipe follows the same structure:

```markdown
## Recipe title

### Context

One paragraph: what problem does this solve, and when would you use it.

### Code

\```awk
# Working AWK code with comments
\```

### Explanation

- Bullet points explaining non-obvious lines.

### Variants

Optional: alternative approaches, edge cases, one-liner versions.
```

**Rules:**

- Code must work on the test data in `docs/data/`.
- Mark gawk-specific recipes with **[gawk]** in the heading or above the code block.
- Provide a POSIX alternative when feasible.
- Use `BEGIN { FS = "\t" }` explicitly for tab-separated formats (GTF, BED, GFF, VCF).
- Add `LC_ALL=C` in bash examples when regex matches fields with `;` or `|`.

## Adding a test

Every recipe with non-trivial logic should have a corresponding test in `docs/tests/`.

Tests use [bats-core](https://github.com/bats-core/bats-core):

```bash
#!/usr/bin/env bats
# docs/tests/my_chapter.bats

DATA="docs/data"

@test "Description of what is being tested" {
    run bash -c "awk '...' $DATA/input_file | wc -l | tr -d ' '"
    [ "$status" -eq 0 ]
    [ "$output" = "expected_value" ]
}
```

Run all tests from the repo root:

```bash
bats docs/tests/
```

Use `gawk` explicitly in tests that require gawk features.

## Commit conventions

```
feat: add new recipe or chapter
fix:  correct an existing recipe or test
docs: update README, index, or non-recipe content
```

Keep commits atomic: one logical change per commit.

## Test data

If your recipe needs new test data, add it to `docs/data/` and document it in
`docs/data/README.md`. Keep files small (< 50 KB), synthetic, and consistent with
existing sample IDs (`NA12878`, `NA19238`, `NA20585`) and chromosomes (`chr1`, `chr2`).

## What not to do

- Don't add recipes that duplicate existing tools without adding value (e.g. don't
  rewrite `bedtools intersect` in AWK unless the recipe teaches something).
- Don't add Python/R code — this is an AWK cookbook. Reference other tools in prose
  when AWK is the wrong choice.
- Don't commit without running `bats docs/tests/` first.

## Questions?

Open an [issue](https://github.com/chrstphdm/awk-bioinfo-cookbook/issues).
