## R CMD check results

0 errors | 0 warnings | 0 notes (local macOS `--as-cran`).

## Test environments

* local macOS (Darwin 25.3.0), R 4.5.2 — `R CMD check --as-cran`
* win-builder R-devel
* win-builder R-release
* macOS builder (mac.r-project.org)
* R-hub: linux, windows, macos, macos-arm64

## Changes since 0.4.3

This is a feature/maintenance release moving 0.4.3 → 0.4.9. Highlights:

* New higher-order infrastructure: `chain_structure()`, `markov_order_test()`,
  `markov_stability()`, `passage_time()`, `path_dependence()`, and
  `transition_entropy()`, each with full `print`/`summary`/`plot` methods and
  `netobject_group` dispatch.
* New `build_mlvar()` (multilevel VAR) and `build_mmm()` (mixed Markov models)
  with byte-equivalence to `mlVAR` and `seqHMM` respectively.
* `build_mcml()` correctness fixes: model-derived weights, symmetrize,
  labels-vs-data alignment; `cluster_summary()` matrix path now aggregation-only
  (`type =` removed — see DESCRIPTION).
* Bootstrap gains a `boundary` argument for inclusive vs strict consistency-range
  tests.
* New `plot_state_frequencies()` generic with marimekko / bars / mosaic styles.
* Documentation cleanup: ASCII-only R sources; all exported argument lists
  now match Rd `\usage`; `as_tna` example updated for new `cluster_summary`
  signature; Suggests trimmed (removed `BiasedUrn`, `HyperG`, `RSpectra`,
  `pkgdown` — none referenced in installed code, tests, or vignettes).

## Notes on URL checks

The DOI `https://doi.org/10.1145/3706468.3706513` (ACM) returns HTTP 403
when accessed by automated URL checkers. The DOI is valid and resolves
correctly in a browser. ACM restricts automated access to their landing
pages. The same DOI format (`<doi:10.1145/3706468.3706513>`) is used in
the DESCRIPTION field where it is validated by the CRAN submission
infrastructure.

## Note on spelling

The DESCRIPTION references include author surnames (Saqr, López-Pernas)
which are flagged by the spell checker as misspelled words. These are
proper names and correct as written.

## Downstream dependencies

The package has no reverse dependencies on CRAN.

The companion package `cograph` (visualization) uses Nestimate's
`cograph_network` class objects but does not formally depend on Nestimate.
