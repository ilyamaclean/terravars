# terravars 0.0.1

## Bug fixes

* `tests/testthat.R` referenced a nonexistent package (`topoclim`), so the
  test suite never ran under `R CMD check`; fixed to reference `terravars`.
  A stray duplicate file at the package root has been removed.
* `test-skyview.R`, `test-solar.R`, and `test-wind.R` called an old,
  nonexistent API (`sky_view_factor()`, `solar_index()`,
  `wind_shelter(wind_direction =, max_dist =)`); rewritten against the
  current `skyview()`/`solarindex()`/`wind_shelter()` functions.
* `twi(method = "standard")` no longer silently produces a large finite
  value on perfectly flat ground; it now produces `Inf`, as documented.
* `.check_projected()` referenced `terra::is.projected()`, which does not
  exist, and was never called by anything; it now uses
  `terra::is.lonlat()` and is wired into every function documented as
  requiring a projected CRS (`solarindex()`, `horizon()`, `skyview()`,
  `twostream()`, `canopy3d()`, `wind_shelter()`, `wind_altitude()`,
  `coastalexposure()`).
* `twostream()`, `skyview()`, and `canopy3d()` now validate the documented
  `lref + ltra < 1` constraint instead of silently returning `NaN` when
  it's violated.
* `basindelin()`, `basinmerge()`, and `flowacc()` now validate their
  `SpatRaster` arguments (previously the only exported functions in the
  package with no input validation at all); `basinmerge()`'s `bsn` and
  `flowacc()`'s `weight`/`bsn` are now checked against `dtm` for matching
  extent/resolution/CRS via `.check_compatible()`.
* `dtm100m`, the example DTM used throughout the package's own
  documentation, was itself entirely undocumented; it now has a proper
  help page.

## Internal / cleanup

* Removed `renumberbasin()` and `two_stream_cpp()` (unused C++ entry
  points with no R caller).
* Removed `basinCpp()`'s unused `dun` parameter (and the corresponding
  unused array built for it in `basindelin()`).
* `DESCRIPTION` no longer claims a `RcppArmadillo` dependency that isn't
  actually used or linked; `stats` and `utils` (used via `::`) are now
  declared in `Imports`.
* Minor formatting cleanup in `R/wetness.R` for consistency with the rest
  of the package; documented the `"saga"` method's internal contributing-
  area cap in `twi()`'s Details.
* Added dedicated tests for `basindelin()`, `basinmerge()`, and
  `fillsinks()` (previously untested), then factored the "pad with a 9999
  border" boilerplate repeated across `basindelin()`, `basinmerge()`,
  `fillsinks()`, and `flowacc()` into a single shared `.pad1()` helper
  (`R/utils.R`), with the new tests as a regression safety net.
* Removed `inst/extdata/dtm100m.tif`, an unused bundled file not
  referenced anywhere in `R/`, `tests/`, or `vignettes/`.
