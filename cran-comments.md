## Test environments
* local: macOS Sonoma 14.2.1, R 4.5.2, `R CMD check --as-cran`
* win-builder: R-release, R version 4.5.2 Patched (2026-02-13 r89426 ucrt), Windows x86_64-w64-mingw32
* win-builder: R-devel, R Under development (unstable) (2026-02-25 r89481 ucrt), Windows x86_64-w64-mingw32
* win-builder: R-oldrelease, R version 4.4.3 Patched (2026-02-12 r89426 ucrt), Windows x86_64-w64-mingw32

## R CMD check results
* local (`--as-cran`): 0 errors | 0 warnings | 0 notes
* win-builder R-release: 0 errors | 0 warnings | 1 note
* win-builder R-devel: 0 errors | 0 warnings | 1 note
* win-builder R-oldrelease: 0 errors | 0 warnings | 1 note

## Notes
* This is a new submission.
* The only NOTE is `Possibly misspelled words in DESCRIPTION`.
* The flagged terms are an intentional domain-specific acronym (DIF) and a proper name (author name), and are not typos.

## Downstream dependencies
No reverse dependency checks were run.
