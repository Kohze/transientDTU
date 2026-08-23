# Formal transientDTU decision rule

## Input unit

The atomic input is one feature, focal group, comparison group, and ordered
stage. It contains a signed focal-minus-comparator usage effect and an adjusted
p-value. Gene identifiers and optional gene-level adjusted p-values accompany
the feature.

## Stage state

For every feature, focal group, and stage, pairwise rows are collapsed into one
stage state. A non-zero state requires:

1. at least `min_comparators` distinct comparator groups;
2. finite adjusted p-values strictly below `q_threshold` for every required
   focal-versus-comparator test;
3. absolute effects at least `effect_threshold` for every comparator;
4. identical non-zero effect signs;
5. comparator spread strictly below `comparator_tolerance`; and
6. when requested, a finite gene-level adjusted p-value strictly below
   `gene_q_threshold`.

Comparator spread is `max(effect) - min(effect)`. Because every effect is
focal minus comparator, this equals the range of the comparison-group means.
The state is `1` for higher focal usage, `-1` for lower focal usage, and `0`
otherwise.

## Episode

Consecutive equal non-zero states form a candidate run. The run length must be
within `min_episode_stages` and `max_episode_stages`. A run is retained only
when the requested `flank_width` stages on each side are observed and the
maximum absolute focal-versus-comparator effect over those flanks is strictly
below `flank_tolerance`. Thus boundary runs are not episodes under the default
rule.
Missing stage rows break runs and cannot serve as reconverged flanks. The
default legacy-compatible `flank_missing = "available"` rule ignores non-finite
effect summaries at an otherwise observed flank, requires at least one finite
flank summary, and records completeness. The stricter `flank_missing =
"complete"` sensitivity requires finite summaries at every requested flank.
Defaults of one or more episode stages and one flank stage on each side exactly
retain the paper rule. Optional strictly increasing stage coordinates annotate
elapsed spans but do not change ordinal adjacency.

The episode key is feature, gene, focal group, start stage, end stage, and
direction. The weakest component evidence is the maximum adjusted p-value in
the run. Effects and stage counts are retained explicitly.

## Replicate separation

For every episode stage, a higher episode has separation
`min(focal replicates) - max(comparator replicates)`. A lower episode has
separation `min(comparator replicates) - max(focal replicates)`. The episode
separation is the minimum of these values across its stages. Strictly positive
separation means every focal replicate lies beyond every comparator replicate
in the required direction at every episode stage.

This is `replicate_method = "complete"`, the default. The optional `"quantile"`
method substitutes inner tail quantiles and `"median"` substitutes group
medians. They are descriptive robustness summaries, not inferential tests.
Missing values either make an episode incomplete or are omitted subject to a
minimum finite replicate count in the focal and every comparator group.

## Reciprocal exchange

An episode belongs to a reciprocal exchange when the same gene, focal group,
start stage, and end stage contains at least one higher and one lower feature
episode. This is a compositional pattern annotation, not an additional
hypothesis test.

## Candidate ranking

Episodes are grouped by gene, focal group, start stage, and end stage. Groups
are ordered by:

1. maximum episode weakest-component adjusted p-value, ascending;
2. maximum absolute usage effect, descending;
3. minimum replicate separation, descending when available;
4. number of episode stages, descending;
5. gene identifier, focal group, start-stage index, end-stage index, and a
   stable concatenation of feature identifiers, all ascending.

When one row per gene is requested, the first deterministically ordered row for
each gene is retained before applying the requested panel size.
