# transientDTU 0.99.2

- Adopt the MIT License across package metadata and repository documentation.
- Add `simulateMethylationExpression()` and `plotMethylationExpression()` for
  reproducible example points and descriptive association, trajectory, and
  CpG-state panels modeled on the original thesis application.
- Link the README to the original thesis application motivating the new panel
  interface.
- Simplify README rendering by removing the zoomable workflow diagram.
- Revalidate the source package, examples, tests, and rebuilt vignettes.

# transientDTU 0.99.1

- Initial Bioconductor submission version.
- Generic pairwise-DTU input contract for ordered multi-group designs.
- Bounded diverge-reconverge episode detection with explicit edge-case rules.
- Optional replicate-separation and reciprocal-exchange annotations.
- Deterministic gene-level candidate ranking.
- Joint threshold-grid stability analysis.
- Known-truth simulation and benchmarking helpers.
- `SummarizedExperiment` interoperability, standard Bioconductor example data,
  and a formal result class.
- Reviewer-facing workflow documentation, complete related-method references,
  and an installed recipe for reproducing the seeded example data.
- Correct preservation of gene identifiers when simulated gene-level adjusted
  p-values are mapped back to pairwise rows.
- Structured input audits for schema, directed-contrast coverage, and
  replicate-cell completeness.
- First-class two-group operation through `min_comparators = 1`.
- Optional episode-length bounds, multi-stage flanks, and irregular-stage
  coordinate annotations.
- Complete, quantile, and median replicate summaries with explicit missing
  value and minimum-replicate policies.
- Faster stage aggregation; a 96,000-row synthetic table improved from about
  27 seconds to about 3 seconds in the local single-process benchmark.
- Formal linkage to the companion framework article, including a paper-output
  regression that reproduces all 1,348 archived episodes and the exact ordered
  six-gene panel.
- An explicit installed software citation matching the manuscript bibliography.
