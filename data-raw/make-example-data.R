# Recreate the package's synthetic example dataset.
# Run from the package root after loading the development package.

transientExample <- simulateTransientDTU(
    n_genes = 12,
    n_stages = 8,
    n_replicates = 2,
    signal = 0.35,
    noise_sd = 0.01,
    seed = 20260823
)

save(
    transientExample,
    file = file.path("data", "transientExample.rda"),
    compress = "xz",
    version = 3
)
