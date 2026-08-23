# transientDTU

[![R package check](https://github.com/kohze/transientDTU/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kohze/transientDTU/actions/workflows/R-CMD-check.yaml)
[![Bioconductor devel](https://github.com/kohze/transientDTU/actions/workflows/bioconductor-devel.yaml/badge.svg)](https://github.com/kohze/transientDTU/actions/workflows/bioconductor-devel.yaml)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](LICENSE.md)

`transientDTU` detects explicit, bounded **diverge-reconverge transcript-usage
episodes** in ordered multi-group studies. It is a post-inference decision
layer: upstream tools estimate differential transcript usage; `transientDTU`
turns the complete stage-specific evidence into auditable episodes and a stable
candidate ranking.

## Why use it?

Pairwise DTU result tables answer whether usage differs in individual
contrasts. They do not directly encode the more specific temporal question:

> Does one group separate coherently from the required alternatives, remain
> separated for one or more consecutive stages, and reconverge on both sides?

`transientDTU` makes that rule explicit and reusable. It provides:

- method-agnostic long-table input;
- preflight schema, coverage, and replicate-cell auditing;
- any number of ordered stages and two or more groups;
- component-q, effect, comparator-coherence, and flank requirements;
- configurable episode length, flank depth, and irregular-time annotations;
- explicit missing-data and boundary behavior;
- complete, quantile, or median replicate-level robustness checks;
- reciprocal exchange annotation;
- deterministic one-row-per-gene ranking;
- joint threshold-grid stability analysis;
- known-truth simulations and operating-characteristic summaries; and
- native `SummarizedExperiment` interoperability.

## Installation

The development version can be installed from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("kohze/transientDTU")
```

After acceptance into Bioconductor, install the release version with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("transientDTU")
```

## Quick start

```r
library(transientDTU)
data("transientExampleSE")
data("transientExamplePairwise")
stage_order <- S4Vectors::metadata(transientExampleSE)$stage_order

auditTransientInput(
    transientExamplePairwise,
    stage_order,
    usage = transientExampleSE
)

result <- runTransientDTU(
    transientExamplePairwise,
    stage_order = stage_order,
    usage = transientExampleSE,
    gene_name_col = "gene_name",
    gene_q_col = "gene_q",
    gene_q_threshold = 0.05,
    panel_size = 6
)

result
candidateTable(result)
```

`transientExampleSE` follows the standard Bioconductor feature-by-sample
layout: usage values are in an assay, feature annotations are in `rowData`,
sample design variables are in `colData`, and experiment-level declarations
are in `metadata`. `transientExamplePairwise` is an `S4Vectors::DataFrame`
containing the matching generic contrast table. Both are seeded, synthetic,
and safe to use in examples and offline integration tests.

Upstream result columns can have arbitrary names; map them through the column
arguments of `makeStageDTU()` or `runTransientDTU()`. The signed effect must
always be `focal usage - comparator usage`. See the
[upstream input guide](inst/UPSTREAM_INPUT.md) for the required contract and
method-specific reshaping guidance.

Two-group studies use `min_comparators = 1`; the paper-compatible multi-group
default remains two comparators. See the
[20-researcher use-case analysis](USE_CASE_ANALYSIS.md) for supported,
conditional, and deliberately out-of-scope designs.

## Scope and statistical boundary

The package does not claim episode-level false-discovery-rate control. It
applies a transparent biological decision rule after upstream inference.
Simulation, full-family multiplicity adjustment, joint sensitivity analysis,
and independent validation remain necessary when the output supports a
scientific claim.

The default `flank_missing = "available"` option reproduces the motivating
paper rule. It records flank completeness for every episode. Use
`flank_missing = "complete"` when both finite flank summaries are required.

## Development status

Version 0.99.0 is prepared for Bioconductor review. The package includes unit
tests, executable vignettes, an optional full-paper regression test, GitHub
Actions, formal algorithm documentation, and AI-assistance provenance.

## Licence

Artistic License 2.0.
