# Twenty-researcher use-case analysis

This is a design stress test, not a claim that twenty researchers used or
endorsed the software. Each researcher represents a plausible combination of
data and aim. Scores are: **3 strong**, **2 conditional**, **1 limited**, and
**0 inappropriate**. “Before” describes the initial pre-audit implementation
audit; “after” includes the robustness work recorded in the current commit.

| # | Researcher, data, and aim | Before | After | Neutral assessment after improvement |
|---:|---|---:|---:|---|
| 1 | Developmental biologist; bulk RNA-seq, three lineages, eight stages; find bounded lineage-specific DTU | 3 | 3 | This is the reference design. Generic upstream contrasts, strict paper defaults, replicate checking, reciprocal annotation, and ranking all apply directly. |
| 2 | Clinical researcher; treated versus control over six visits; find temporary treatment responses | 2 | 3 | Two groups are now explicit with `min_comparators = 1`, including simulation and audit support. Claims still depend on appropriate repeated-measures inference upstream. |
| 3 | Atlas consortium; six tissues across twelve stages; find tissue-restricted episodes | 2 | 3 | Arbitrary group counts work; the faster aggregation removes a major scaling penalty. Complete directed contrasts grow quadratically in group count and remain the user's responsibility. |
| 4 | Single-cell analyst; donor-pseudobulk transcript usage by cell type and condition | 2 | 3 | Strong when donors, not cells, are replicate units and upstream inference models donor/batch structure. `SummarizedExperiment` and replicate-cell auditing fit this workflow. |
| 5 | Single-cell analyst treating thousands of cells as independent biological replicates | 0 | 0 | Still inappropriate. The package cannot repair pseudoreplication; pseudobulk or a valid hierarchical upstream model is required. |
| 6 | Chronobiologist; irregular sampling at 0, 1, 2, 4, 8, and 16 hours | 2 | 3 | Strictly increasing `stage_coordinates` now annotate elapsed spans. Detection remains ordinal, so a four-hour and an eight-hour gap are adjacent by design rather than interpolated. |
| 7 | Splicing researcher; interested only in one-stage pulses | 2 | 3 | `min_episode_stages = 1` and `max_episode_stages = 1` express this aim without post-hoc filtering. |
| 8 | Regeneration researcher; wants sustained two-to-four-stage transient programs | 2 | 3 | Inclusive episode-length bounds directly encode the target duration. |
| 9 | Researcher requiring two stable stages before and after every event | 1 | 3 | `flank_width = 2` enforces four observed flanks. Wider flanks reduce eligibility near boundaries, as they should. |
| 10 | Large cohort with occasional outlying usage estimates; wants a robustness screen | 1 | 3 | Quantile or median replicate summaries avoid letting one observation determine the descriptive score. They are not substitutes for inferential modelling. |
| 11 | Small cohort with sporadic missing usage; wants to retain evaluable events | 1 | 2 | `missing_values = "omit"` plus `min_replicates` is explicit and auditable. Results remain conditional because missingness mechanisms can bias the screen. |
| 12 | Rare-disease study with one sample per group-stage cell | 2 | 2 | Discovery can run without replicate filtering, and the audit flags single-sample cells. Biological replicate consistency cannot be claimed. |
| 13 | Analyst receiving an incomplete satuRn/DEXSeq contrast export | 1 | 3 | `auditTransientInput()` locates incomplete focal-stage families before discovery and reports coverage by focal group and stage. It cannot reconstruct missing inference. |
| 14 | Investigator scanning only one prespecified focal lineage against all others | 2 | 3 | Targeted focal tables are now audited without demanding unused reverse focal directions; `expected_focal_groups` makes the intent explicit. |
| 15 | Gene with many isoforms; aim is reciprocal isoform exchange | 3 | 3 | Any feature count per gene is supported. Exact shared intervals with higher and lower features are annotated as reciprocal. |
| 16 | Researcher interested in non-reciprocal single-transcript episodes | 3 | 3 | Set `reciprocal_only = FALSE`; ranking does not demote the underlying thesis figures or force reciprocal biology. |
| 17 | Researcher interested in several episodes from the same gene | 2 | 3 | Detection already returns all episodes; set `unique_genes = FALSE` to retain multiple ranked gene-interval rows. |
| 18 | Researcher studying monotonic divergence or an event touching the first/last stage | 0 | 0 | Deliberately out of scope. Those patterns lack observed reconvergence on both sides and require a different trajectory model. |
| 19 | Proteomics or isoform-usage researcher with bounded proportions and ordered groups | 2 | 2 | The generic decision layer can operate on proportions, but terminology, upstream error models, and validation are transcript-focused. This is conditional reuse, not validated domain transfer. |
| 20 | Methods group comparing several upstream DTU engines on genome-scale tables | 2 | 3 | All engines can map to the same atomic contract, input audits expose coverage differences, and threshold stability supports common decision settings. Upstream model calibration and multiplicity remain outside the package. |

## What changed because of the panel

Five changes recurred across otherwise different researchers:

1. `auditTransientInput()` now reports schema errors, feature-to-gene conflicts,
   duplicate contrasts, incomplete directed families, sample-cell counts,
   missing usage, and the inference properties that require manual confirmation.
2. Two-group designs are supported explicitly with `min_comparators = 1`; the
   paper-compatible default remains two comparators.
3. Episode duration, flank depth, and irregular numeric stage coordinates are
   represented in the detector and result tables.
4. Replicate robustness can use complete separation, inner quantiles, or
   medians, with explicit missing-value and minimum-replicate policies.
5. Stage aggregation was vectorized. On the local Windows/R 4.4 single-process
   benchmark, a synthetic 1,000-gene table with 96,000 atomic rows and 48,000
   stage cells fell from 27.2 seconds to 2.8 seconds; episode detection took
   2.2 seconds. These timings are diagnostic, not cross-platform guarantees.

## Boundaries that were not “fixed”

The audit did not turn unsuitable questions into supported ones. The package
still does not infer DTU from counts, correct pseudoreplication, model subjects
or batches, verify that q-values use the right family, provide episode-level
FDR control, call boundary/monotonic trajectories, or establish biological
validation. Those boundaries are scientifically important and remain explicit.
