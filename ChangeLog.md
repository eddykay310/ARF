## v2.2 — 2026-03-31
### Added
  - `ARF-package.R` with `@importFrom` declarations for `grDevices`, `stats`, and `utils`
    functions, resolving `R CMD check` undefined-global warnings.
  - `.Rbuildignore` to exclude non-package directories (`docs`, `data-raw`,`rRNAs`, 
    `test_data`) from the build.
  - `docs/dricARF_use_cases_documentation.md`
  - ssGSEA-based correction of ES2 (NES_randZ) scores.

### Changed
  - `DESCRIPTION`: licence field corrected from `use_gpl3_license()` to `GPL (>= 3)`;
    `bedr` and `tidyverse` removed from `Imports` and replaced with specific packages
    (`dplyr`, `magrittr`, `readr`, `RColorBrewer`, `circlize`, `wesanderson`);
    `HelloRanges` moved to `Suggests`.
  - All roxygen2 documentation overhauled: many fixes across `ARF_platform.R`, `dripARF.R`,
    and `dricARF.R`, including corrected `@param` names, completed `@return` tags,
    fixed `@examples`, resolved typos, and corrected wrong titles/keywords.
  - `README.md` comprehensively updated.

## v2.1
  - All species support for both dricARF and dripARF.
  - New functions for PDB parsing and lift-overing; ARF_parse_PDB_ribosome(), ARF_convert_Ribo3D_pos(), 
    dripARF_get_RP_proximity_sets(), and dricARF_liftover_collision_sets().
  - dricARF output figure updated.

## v2.0
  - Integrating 3D ribosome structure analysis to ARF
  - All relevant functions now require the rRNA fasta file that is used for mapping.
  - First release of dricARF - Differential Ribosome Collision Prediction pipeline
  - Bug fixes (closing issues until now)

## v1.0
  - Publication release
  - A few bug fixes

## v0.9.1 - 2022-04-14
  - Paper revision release
  - rRNA-RP contact point distance threshold changed
  - New experimental function that allows altering the Contact Point threshold and possible directionality in predictions
  - A few technical updates
  - RP-rRNA distance matrices are added to the repo

## v0.9 - 2021-11-26
  - First release
  - Paper submission version
