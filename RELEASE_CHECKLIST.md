# Release checklist

The codebase is locally release-clean. Complete these repository- and
author-owned steps before Bioconductor submission:

- [ ] Confirm Robin Gounder's maintainer email and add an ORCID to `DESCRIPTION`
      if available.
- [ ] Confirm Russell Hamilton's contributor credit and preferred publication
      metadata.
- [ ] Create `kohze/transientDTU` on GitHub and push the local `main` branch.
- [ ] Confirm all three GitHub Actions workflows pass; the Bioconductor devel
      workflow is the authoritative current-environment check.
- [ ] Enable GitHub Pages with GitHub Actions as the source before running the
      pkgdown workflow.
- [ ] Confirm `transientDTU` remains unused, case-insensitively, in current,
      archived, and deprecated CRAN/Bioconductor packages.
- [ ] Create a tagged GitHub release only after metadata approval.
- [ ] Submit the public GitHub repository through the Bioconductor contribution
      portal and respond to review in the package issue.

Do not commit the motivating study data. The optional regression script reads
those files only through `TRANSIENTDTU_PAPER_DIR` and the repository includes
only hashes and aggregate expected results.
