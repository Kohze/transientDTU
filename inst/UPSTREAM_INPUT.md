# Preparing upstream differential-usage evidence

`transientDTU` is a post-inference decision layer. It does not require a
particular DTU model, but it does require the complete, directed contrast
family used for discovery.

## Atomic input contract

Each row describes one feature, one focal group, one different comparator, and
one stage. Required fields are:

| Meaning | Default column | Requirement |
|---|---|---|
| feature identifier | `feature_id` | one gene mapping only |
| gene identifier | `gene_id` | stable across stages |
| focal group | `focal_group` | at least two groups overall |
| comparator group | `comparator_group` | different from focal |
| ordered stage | `stage` | present in `stage_order` |
| signed usage effect | `effect` | focal minus comparator |
| adjusted component p-value | `q_value` | finite and in [0, 1] |

An optional `gene_q` column can impose an upstream gene-level screen. Adjust
component p-values across the full family that will be searched, not separately
within genes or only among selected candidates.

Supply both directions only when both groups may act as the focal group. Do not
manufacture a missing direction by negating an effect unless the associated
inferential statistic is valid for that directed contrast. Duplicate atomic
rows are rejected.

Run `auditTransientInput()` before discovery. For a two-group study, set
`min_comparators = 1`. For a targeted scan containing only prespecified focal
groups, provide `expected_groups` and optionally `expected_focal_groups`; reverse
focal directions are not required unless they are part of the intended search.

## Common upstream methods

- **satuRn:** extract transcript-level contrasts for every group pair at every
  stage, retain the model's adjusted transcript-level values, and reshape the
  contrast label into focal and comparator columns.
- **DEXSeq:** test the required stage-specific group contrasts, export feature
  identifiers, gene identifiers, signed usage or interaction effects, and an
  adjustment performed over the complete searched family.
- **DRIMSeq:** use its gene/transcript tests as upstream evidence, then construct
  directed stage-specific contrasts from the fitted model. A single omnibus
  p-value cannot supply direction or comparator coherence by itself.
- **limma:** `fitPairwiseUsage()` is a convenience adapter for continuous usage
  matrices and a saturated group-by-stage cell-means model. It is useful for
  prototyping and processed usage data, not a replacement for count-aware DTU
  inference when transcript counts are available.

Always record the upstream model, normalization, contrast matrix, multiplicity
family, feature filtering, and software versions. The output does not carry
episode-level FDR control merely because its input contains adjusted p-values.
For single-cell or repeated-measures studies, cells and repeated observations
must not be treated as independent biological replicates unless the upstream
design genuinely supports that assumption.
