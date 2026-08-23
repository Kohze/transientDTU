# transientDTU 0.99.0

- Initial Bioconductor submission version.
- Generic pairwise-DTU input contract for ordered multi-group designs.
- Bounded diverge-reconverge episode detection with explicit edge-case rules.
- Optional replicate-separation and reciprocal-exchange annotations.
- Deterministic gene-level candidate ranking.
- Joint threshold-grid stability analysis.
- Known-truth simulation and benchmarking helpers.
- `SummarizedExperiment` interoperability and a formal result class.
- Structured input audits for schema, directed-contrast coverage, and
  replicate-cell completeness.
- First-class two-group operation through `min_comparators = 1`.
- Optional episode-length bounds, multi-stage flanks, and irregular-stage
  coordinate annotations.
- Complete, quantile, and median replicate summaries with explicit missing
  value and minimum-replicate policies.
- Faster stage aggregation; a 96,000-row synthetic table improved from about
  27 seconds to about 3 seconds in the local single-process benchmark.
