# CRAN release checklist

The codebase is locally release-clean. Complete these maintainer- and
infrastructure-owned steps before uploading the source tarball:

- [ ] Confirm Robin Gounder's maintainer email accepts unfiltered CRAN mail and
      add an ORCID to `DESCRIPTION` if available.
- [ ] Confirm Russell Hamilton's contributor credit and agreement to distribute
      the credited material under the package licence.
- [ ] Confirm `transientDTU` remains unused, case-insensitively, in current and
      archived CRAN packages and current or deprecated Bioconductor packages.
- [ ] Run `R CMD build` with current R-patched or R-release.
- [ ] Run `R CMD check --as-cran` on that tarball with current R-devel.
- [ ] Submit the same tarball to Win-builder R-devel and obtain no ERRORs,
      WARNINGs, or unexplained NOTEs.
- [ ] Confirm release, devel, and old-release checks pass on Windows, macOS,
      and Linux.
- [ ] If a public repository is created, restore verified `URL` and
      `BugReports` fields in `DESCRIPTION` before the final build.
- [ ] Update `cran-comments.md` with the actual final check environments and
      any unavoidable NOTE explanations.
- [ ] Upload through the CRAN submission form and confirm the validation email.

Do not commit the motivating study data. The optional regression script reads
those files only through `TRANSIENTDTU_PAPER_DIR`; the repository and source
package include only hashes and aggregate expected results.
