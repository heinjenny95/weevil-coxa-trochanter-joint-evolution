# Release Checklist

Use this checklist before making the repository public or archiving it with Zenodo.

## Required Metadata

- Confirm that `CITATION.cff` uses the final GitHub repository URL.
- Replace the placeholder author entry in `CITATION.cff` with the final author list.
- Add the final manuscript citation once available.
- Add the Zenodo DOI after the first archived release.
- Confirm whether the MIT license is correct for the code release.

## Reproducibility

- Decide whether scripts should remain as an analysis record or be made fully rerunnable.
- If full rerunnability is required, replace placeholder roots with project-relative paths or a single configurable root variable.
- Confirm that no raw data, supplementary tables or figure-caption documents are accidentally included.
- Confirm that no temporary backup folders or `.Rhistory` files are included.

## GitHub Release

- Create a GitHub repository, for example `weevil-coxa-trochanter-joint-evolution`.
- Push the local repository to GitHub.
- Create a tagged release, for example `v0.1.0`.
- Connect the repository to Zenodo and archive the release.
- Add the Zenodo DOI badge and citation to `README.md` after Zenodo creates the DOI.

## Journal Submission

- Use this repository only for code availability.
- Keep supplementary tables, captions and figure inventories outside this GitHub code repository.
