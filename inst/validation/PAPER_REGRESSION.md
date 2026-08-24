# Paper-output regression validation

Validation date: 23 August 2026
Package version: 0.99.2
Script: `inst/scripts/validate-paper-regression.R`

The optional local integration test passed against the retained manuscript
outputs. It reconstructed the standardized stage evidence, applied the default
legacy-compatible flank rule, recalculated replicate separation, annotated
reciprocal exchanges, and applied the new deterministic terminal tie-breaks.

Verified results:

- 1,350 statistical episodes before replicate filtering;
- 1,348 replicate-consistent episodes after filtering;
- exact feature, gene, region, interval, direction, replicate and reciprocal
  identities for all 1,348 archived rows;
- numeric agreement within `1e-12` for effect, component-q, flank and replicate
  separation fields; and
- the exact ordered six-gene panel: `Scg3`, `Gpm6a`, `Ntrk2`, `Tecr`, `Armc8`,
  and `Bin1`.

Input SHA-256 values:

| Input | SHA-256 |
|---|---|
| `transient_regional_stage_evaluations_all.csv` | `BDCD747C336FD7688D4DB1253191A3E7E7734A05C2DFEBD1B0F2C7E1713D1067` |
| `transient_regional_filtered_isoform_fractions.csv` | `179D59263A26122767B9432CAB6F25A786E8BEE5D544056D515E407B9F2B7CB1` |
| `transient_regional_isoform_episodes.csv` | `5B872ABCD8DA6926626005FCEFEC162AF3AB719E8E33A8BBC57E17E9E6BB3E69` |
| `transient_regional_top_candidates.csv` | `ED1E5D95C37A6FD11C30F20F1143A5E8720A6F96F82BA29769BFF1E99B497903` |

The source data are not bundled with the package. This report records the
result and hashes without expanding the software package's redistribution
scope.
