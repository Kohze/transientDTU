# Contributing to transientDTU

Bug reports, reproducible edge cases, documentation corrections, and focused
pull requests are welcome. Please search existing issues first and use a small
synthetic example whenever real data cannot be shared.

Before opening a pull request:

1. Create a branch from `main` and keep the change focused.
2. Add or update a `testthat` test for changed behaviour.
3. Run `devtools::document()` and `devtools::test()`.
4. Run `R CMD build .` and `R CMD check <tarball>`.
5. Run `R CMD check --as-cran` on the tarball with current R-devel. Then run
   BiocCheck as a supplementary interoperability audit:

   ```r
   BiocCheck::BiocCheck(
       "transientDTU_1.0.0.tar.gz",
       `new-package` = TRUE
   )
   ```

6. Update `NEWS.md` when user-visible behaviour changes.

Do not commit identifying study data, credentials, generated check folders, or
large derived files. Contributions must be available under Artistic-2.0. State
the use of generative AI in the pull request when it materially assisted code,
tests, documentation, or review, and verify every contributed line.
