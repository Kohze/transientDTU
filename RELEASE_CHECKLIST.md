# Bioconductor submission checklist

The codebase is locally release-clean. Complete these maintainer- and
infrastructure-owned steps before opening a Bioconductor review:

- [ ] Confirm Robin Gounder's maintainer email accepts Bioconductor build and
      review mail and add an ORCID to `DESCRIPTION` if available.
- [ ] Confirm Russell Hamilton's contributor credit and agreement to distribute
      the credited material under the package licence.
- [ ] Confirm `transientDTU` remains unused, case-insensitively, in current or
      past CRAN and Bioconductor packages.
- [ ] Create the public `kohze/transientDTU` GitHub repository, push the
      package-only `main` branch, and then add verified `URL`, `BugReports`,
      and `repository-code` metadata.
- [ ] Run `R CMD build`, `R CMD check`, and `BiocCheck` in the current
      Bioconductor devel environment and address or justify every finding.
- [ ] Review the BiocCheck function-length NOTE. Prioritize extracting clear,
      testable helpers from the validation-heavy functions, or explain why a
      sequential implementation is safer where no clean boundary exists.
- [ ] Recreate the GitHub Actions, issue templates, and pkgdown configuration
      on a separate infrastructure branch, as required by the current
      Bioconductor submission instructions, and confirm its checks pass.
- [ ] Consider switching the vignettes to `BiocStyle::html_document` once
      `BiocStyle` is available in the development environment.
- [ ] Decide whether to add a licence-compatible public real-data vignette;
      the seeded synthetic example remains the offline test fixture.
- [ ] Submit the public repository through the Bioconductor contribution portal
      and respond promptly in the public review issue.
- [ ] After the paper receives a DOI, make its article citation the primary
      entry returned by `citation("transientDTU")`.

Do not commit the motivating study data. The optional regression script reads
those files only through `TRANSIENTDTU_PAPER_DIR`; the repository and source
package include only hashes and aggregate expected results.
