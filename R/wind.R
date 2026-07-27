# ============================================================
# wind.R
#
# Wind and coastal-exposure functions -- grouped in one file since
# they're all small, complementary inputs to wind/climate downscaling:
#   - wind_shelter()      topographic wind shelter coefficient
#   - wind_altitude()     orographic speed-up (elevation) coefficient
#   - coastalexposure()   proportion sea upwind (coastal exposure)
# ============================================================

#' Wind shelter coefficient
#'
#' Estimates how much a point's wind speed is reduced by the terrain
#' upwind of it, for one or more wind directions. For each point, this
#' looks back the way the wind is coming from and finds the steepest
#' angle up to any upwind terrain within a search distance - the same
#' terrain horizon calculation used by [horizon()] elsewhere in this
#' package - then converts that angle into a multiplier you can apply
#' directly to a wind speed: close to 1 for an exposed spot with nothing
#' upwind to block the wind (e.g. a ridge top or open field), falling
#' towards 0 for a spot tucked behind higher ground in that direction
#' (e.g. downwind of a ridge or in the lee of a hill).
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must use a projected CRS with metre units. Required unless
#'   `hor` is supplied.
#' @param wdir Numeric vector. The direction(s) the wind is blowing *from*
#'   (degrees, N = 0, clockwise) - e.g. 270 for a wind coming from the
#'   west. Supplying more than one value returns one layer per direction.
#' @param obs_hgt Numeric. Height above ground (metres) at which the
#'   coefficient should apply - e.g. `10` for a wind speed measured or
#'   modelled at a standard 10 m mast, rather than right at the surface.
#'   Passed through to [horizon()]: a higher observer sees over more of
#'   the upwind terrain, so shelter falls off (and the coefficient rises
#'   towards 1) as `obs_hgt` increases. Default `0` (ground level). Only
#'   used when `hor` is not supplied.
#' @param hor Optional multi-layer `SpatRaster` of horizon angles already
#'   computed by [horizon()], covering (at least) the azimuths implied by
#'   `wdir` (i.e. `wdir \%\% 360`) - supply
#'   `hor = horizon(dtm, azi = wdir \%\% 360, ...)`. Useful when
#'   you've already computed a horizon stack for another purpose, or want
#'   a faster (if less exhaustive) search via [horizon()]'s own `nstep`
#'   rather than the full search `wind_shelter()` otherwise always does -
#'   `dtm` and `obs_hgt` are both ignored when this is supplied.
#'
#' @return A `SpatRaster` with one layer per wind direction, named
#'   `wind_shelter_<wdir>`.
#' @export
#'
#' @seealso [horizon()]
#'
#' @details
#' The upwind search itself is exactly [horizon()]'s own default: every
#' intervening cell out to the edge of the raster is checked, in whichever
#' direction(s) `wdir` implies. There is no shorter-range option exposed
#' here directly - if a faster, sparser search matters more than
#' exhaustiveness (e.g. on a very large raster), compute it yourself via
#' [horizon()]'s own `nstep` and pass the result in through `hor` instead.
#' The resulting angle is always \eqn{\ge} 0 - a point with
#' upwind terrain *below* it is treated as fully open (angle 0), not as
#' having a negative angle. This is a deliberate consequence of reusing
#' [horizon()]'s own convention (needed so its output can also feed
#' [skyview()] and [solarindex()]), and differs from the original
#' Winstral et al. (2002) "Sx" index, which allows negative angles to
#' flag terrain that sits *above* its surroundings (e.g. a ridge top) as
#' extra-exposed rather than merely unsheltered. If you need that
#' distinction preserved, this function is not currently a substitute for
#' it.
#'
#' The angle-to-coefficient transform is:
#' \deqn{\text{coef} = 1 - \arctan(17 \tan(\text{angle})) / 1.65}
#' following Ryan (1977). It maps an angle of 0 (no upwind obstruction) to
#' a coefficient of 1 (full wind speed), falling off steeply even for
#' shallow horizon angles of a few degrees, and asymptoting towards (but
#' never quite reaching) 0 as the angle approaches 90.
#'
#' @references
#' Ryan, B.C. (1977). A mathematical model for diagnosis and prediction of
#' surface winds in mountainous terrain. *Journal of Applied Meteorology*,
#' 16, 571-584.
#'
#' Winstral, A., Elder, K. & Davis, R.E. (2002). Spatial snow modeling of
#' wind-redistributed snow using terrain-based parameters. *Journal of
#' Hydrometeorology*, 3, 524-538.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#'
#' # Shelter coefficient from the prevailing south-westerly
#' ws  <- wind_shelter(dtm, wdir = c(225, 270))
#' plot(ws)
#'
#' # Same, but for wind speed at a 10 m mast height
#' ws10 <- wind_shelter(dtm, wdir = c(225, 270), obs_hgt = 10)
#'
#' # Reusing a horizon stack computed once for several directions
#' ha  <- horizon(dtm, azi = c(225, 270) \%\% 360)
#' ws2 <- wind_shelter(wdir = c(225, 270), hor = ha)
#' }
wind_shelter <- function(dtm = NULL, wdir = seq(0, 315, by = 45),
                          obs_hgt = 0, hor = NULL) {
  if (!is.null(dtm)) {
    stopifnot(inherits(dtm, "SpatRaster"))
    .check_projected(dtm)
  }

  # The terrain that can block/slow wind arriving FROM `wdir` sits in that
  # same compass direction -- e.g. wind from the west (wdir = 270) is
  # blocked by a ridge to the west, so we search directly toward `wdir`.
  # This matches horizon()'s own "search toward this azimuth" convention
  # (the same one solarindex() relies on when searching toward the sun's
  # azimuth to find terrain shading it). Earlier versions of this function
  # instead searched wdir + 180 -- the downwind side -- inherited from the
  # original (now-retired) wind_shelter_cpp(); found to be backwards while
  # comparing against coastalexposure()'s C++, which searches azimuth =
  # wdir directly with no flip. See the handoff addendum for the worked
  # example that caught this.
  upwind_az <- wdir %% 360

  if (is.null(hor)) {
    if (is.null(dtm))
      stop("'dtm' must be supplied when 'hor' is not.", call. = FALSE)

    # nstep left at horizon()'s own default (0): always search every
    # intervening cell out to the edge of the raster, rather than trading
    # accuracy for speed via a shorter/sparser search. Callers who want
    # that trade-off can compute horizon() themselves with nstep set and
    # pass it in via 'hor' instead.
    hor <- horizon(dtm, azi = upwind_az, obs_hgt = obs_hgt)
  } else {
    stopifnot(inherits(hor, "SpatRaster"))
  }

  want    <- paste0("horizon_", upwind_az)
  missing <- setdiff(want, names(hor))
  if (length(missing) > 0)
    stop("'hor' is missing layer(s) for the azimuth(s) implied by 'wdir': ",
         paste(missing, collapse = ", "),
         ". Supply hor = horizon(dtm, azi = wdir %% 360, ...).",
         call. = FALSE)
  hor <- hor[[want]]

  # Ryan (1977) shelter transform, ported from the reference
  # implementation: coef = 1 - atan(0.17 * 100 * tan(angle)) / 1.65
  result <- 1 - atan(17 * tan(hor * pi / 180)) / 1.65
  names(result) <- paste0("wind_shelter_", wdir)
  result
}

#' Wind altitude (elevation) coefficient
#'
#' Estimates how much a point's wind speed is increased by sitting on
#' locally elevated, exposed ground - the orographic "speed-up" effect
#' seen on hilltops and ridges, where wind accelerates as it's forced up
#' and over higher terrain. Complements [wind_shelter()], which handles
#' the opposite case (wind slowed by terrain blocking it from a specific
#' direction): together the two give a coefficient that can distinguish a
#' high point exposed to the wind from an equally high point tucked in the
#' lee of something taller (see Details).
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must use a projected CRS with metre units.
#' @param ref_elev Optional numeric scalar. Elevation above sea level
#'   (metres) that the wind speed being downscaled is actually valid at -
#'   e.g. a weather station's elevation. Every cell's coefficient is then
#'   based on how much higher (if at all) it sits above this single fixed
#'   reference. Mutually exclusive with `coarse_dtm`. If neither is
#'   supplied, defaults to the mean elevation of `dtm` itself.
#' @param coarse_dtm Optional `SpatRaster`, a coarser-resolution elevation
#'   model - e.g. matching the grid of the (coarse) wind data being
#'   downscaled. Resampled to `dtm`'s grid with bilinear interpolation and
#'   used as a spatially varying reference elevation instead of one fixed
#'   number: each fine cell's coefficient reflects how much it sticks up
#'   above what the coarse data effectively "saw" at that location.
#'   Reprojected first if its CRS differs from `dtm`'s. Mutually exclusive
#'   with `ref_elev`.
#' @param meas_hgt Numeric. Height above ground (metres) that the wind
#'   speed being downscaled applies at (and that the output coefficient
#'   should therefore also be understood to apply at) - e.g. a standard
#'   2 m or 10 m reference height. Not the same thing as elevation above
#'   sea level (that's `ref_elev`/`coarse_dtm`) or [wind_shelter()]'s
#'   `obs_hgt`, though it plays an analogous role. Default `2`. Values
#'   below 20 cm are raised to 20 cm with a warning (the underlying
#'   log-height-profile relationship performs poorly below that); values
#'   at or below `5.42/67.8` (~8 cm) are rejected outright, since the
#'   formula is undefined there.
#'
#' @return A `SpatRaster` of wind altitude coefficients (\eqn{\ge} 1,
#'   dimensionless), named `wind_altitude`.
#' @export
#'
#' @seealso [wind_shelter()], [horizon()]
#'
#' @details
#' Only elevation *above* the reference has any effect: a cell at or below
#' the reference elevation gets a coefficient of exactly 1 (no speed-up),
#' never less. This is deliberate - any wind *reduction* associated with
#' being at lower elevation is really a consequence of being sheltered by
#' higher terrain in a specific direction, which is what [wind_shelter()]
#' already accounts for. Folding a reduction into this (direction-blind)
#' function as well would double up with, or contradict, whatever
#' [wind_shelter()] says for that same cell and wind direction - so this
#' function only ever pushes the coefficient up, or leaves it at 1, and
#' the two are meant to be multiplied together (e.g. in a downscaling
#' pipeline) rather than used as alternatives.
#'
#' For a cell `ed` metres above the reference elevation, the coefficient
#' follows the same log-height-profile relationship used for a simple
#' measurement-height conversion, treating `ed` as if it were extra
#' measurement height added on top of `meas_hgt`:
#' \deqn{\text{coef} = \max\left(1,\ \frac{\log(67.8(\text{ed} + \text{meas\_hgt}) - 5.42)}{\log(67.8\,\text{meas\_hgt} - 5.42)}\right)}
#' following Ryan (1977).
#'
#' `coarse_dtm` is resampled to `dtm`'s grid with bilinear interpolation,
#' so the reference elevation surface varies smoothly and directly reflects
#' `coarse_dtm`'s own elevation values at every fine cell, rather than
#' duplicating each coarse cell's value uniformly across the whole block it
#' covers. That block-duplication (nearest-neighbour) approach was tried
#' first and discarded: it bakes a hard elevation step into the reference
#' surface at every coarse-cell boundary that has nothing to do with actual
#' terrain, and the log-height-profile transform below amplifies those fake
#' steps into a visibly blocky, checkerboard-like coefficient surface.
#'
#' @references
#' Ryan, B.C. (1977). A mathematical model for diagnosis and prediction of
#' surface winds in mountainous terrain. *Journal of Applied Meteorology*,
#' 16, 571-584.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#'
#' # Default: reference elevation = mean elevation of dtm
#' wa <- wind_altitude(dtm)
#' plot(wa)
#'
#' # Reference elevation from a known station height
#' wa_station <- wind_altitude(dtm, ref_elev = 145)
#'
#' # Reference elevation from a coarser DTM (e.g. matching a coarse wind
#' # model's own grid)
#' dtm_coarse <- aggregate(dtm, fact = 10)
#' wa_coarse  <- wind_altitude(dtm, coarse_dtm = dtm_coarse)
#'
#' # Combine with wind_shelter() for the full picture: exposed high ground
#' # sped up, but only where it isn't also sheltered from this direction
#' ws <- wind_shelter(dtm, wdir = 270)
#' downscale_coef <- wa * ws
#' }
wind_altitude <- function(dtm, ref_elev = NULL, coarse_dtm = NULL, meas_hgt = 2) {
  stopifnot(inherits(dtm, "SpatRaster"))
  .check_projected(dtm)

  if (!is.null(ref_elev) && !is.null(coarse_dtm))
    stop("Supply either 'ref_elev' or 'coarse_dtm', not both.", call. = FALSE)

  # True domain boundary of the log-height-profile formula (67.8*h - 5.42
  # must stay positive); 0.2 m is the looser, "performs poorly below this"
  # floor already documented for the same formula elsewhere in this package.
  floor_hgt <- 5.42 / 67.8
  if (meas_hgt <= floor_hgt)
    stop("'meas_hgt' must be greater than ", round(floor_hgt, 3),
         " m for the log-height-profile formula to be defined.",
         call. = FALSE)
  if (meas_hgt < 0.2) {
    warning("Wind-height profile function performs poorly below 20 cm so ",
            "'meas_hgt' converted to 20 cm.", call. = FALSE)
    meas_hgt <- 0.2
  }

  if (is.null(coarse_dtm)) {
    if (is.null(ref_elev))
      ref_elev <- terra::global(dtm, "mean", na.rm = TRUE)[1, 1]
    ref <- ref_elev
  } else {
    stopifnot(inherits(coarse_dtm, "SpatRaster"))
    if (!terra::same.crs(coarse_dtm, dtm))
      coarse_dtm <- terra::project(coarse_dtm, dtm)
    # Bilinear so the reference elevation surface varies smoothly and uses
    # coarse_dtm's own elevation values directly at every fine cell, rather
    # than duplicating each coarse cell's value uniformly across the whole
    # block it covers (nearest-neighbour bakes a fake hard step into the
    # reference at every coarse-cell boundary, which the log-height-profile
    # transform below then amplifies into a blocky coefficient surface).
    ref <- terra::resample(coarse_dtm, dtm, method = "bilinear")
  }

  ed  <- dtm - ref
  h   <- ed + meas_hgt
  raw <- suppressWarnings(log(67.8 * h - 5.42) / log(67.8 * meas_hgt - 5.42))

  # Only cells genuinely above the reference get a speed-up; everything
  # else (including any NaN the log would produce below its domain, which
  # only ever arises on the ed <= 0 side) is exactly 1, not merely clamped.
  coef <- terra::ifel(ed > 0, raw, 1)
  coef <- terra::ifel(coef < 1, 1, coef)   # safety net for edge-of-domain rounding

  names(coef) <- "wind_altitude"
  coef
}

#' Coastal exposure (proportion sea upwind)
#'
#' Estimates how "coastal" each point is in a specified upwind direction:
#' the proportion of sea versus land sampled along that direction (higher
#' = more exposed to open water), weighted so that nearby coastline
#' matters more than distant coastline. Wind that's
#' travelled a long way over open sea behaves differently (usually
#' milder, moister, and less variable) than wind that's mostly crossed
#' land, so this is intended as an input to downstream wind/climate
#' downscaling rather than a wind speed multiplier in its own right,
#' unlike [wind_shelter()]/[wind_altitude()].
#'
#' @param landsea A `SpatRaster` with `NA` cells representing sea and any
#'   non-`NA` value representing land - the finest-resolution land/sea
#'   mask available. Values are computed for the whole of `landsea`'s own
#'   extent (see Details for how cells near its edge are handled). Must
#'   use a projected CRS with metre units.
#' @param wdir Either a numeric scalar giving the direction (degrees,
#'   N = 0, clockwise) the wind is blowing *from* - the direction
#'   searched for land/sea coverage - or the string `"all"`, which
#'   instead averages the result over 8 directions (`seq(0, 315, by =
#'   45)`, the same set used as [wind_shelter()]'s default) to give an
#'   omnidirectional coastal exposure index. See Details.
#' @param coarse Optional `SpatRaster`, or list of `SpatRaster` objects,
#'   of additional land/sea masks at coarser resolutions than `landsea`
#'   (`NA` = sea, matching `landsea`'s convention). Lets the search reach
#'   much further out than `landsea`'s own extent would allow, without
#'   needing an enormous fine-resolution raster - see Details for how
#'   multiple resolutions are combined. All masks (including `landsea`)
#'   must share the same CRS. Any element that's a wrapped
#'   `terra::PackedSpatRaster` instead of a `SpatRaster` (as returned by
#'   `terra::wrap()` - the form bundled package datasets such as [dtmc]
#'   are stored in, since a plain `SpatRaster` can't be saved to disk) is
#'   unwrapped automatically, so a wrapped dataset can be passed straight
#'   in without unwrapping it by hand first. `NULL` (default): behaves
#'   exactly as if only `landsea` existed.
#' @param n Numeric, roughly 1--3, greater than 0. Controls how strongly
#'   nearby land/sea coverage is weighted relative to distant coverage -
#'   see Details. Default `2` (the original calibration). `1` weights the
#'   near and far field more evenly; higher values weight the near field
#'   more heavily still. The search always reaches all the way out to the
#'   edge of the available data regardless of `n` - lower values of `n`
#'   just need more sample points to cover that same distance evenly,
#'   rather than stopping any sooner.
#' @param jitter Logical. If `TRUE`, each sampled point's search azimuth
#'   is randomly perturbed by up to +/-10 degrees, instead of every sample
#'   following an exact straight line at azimuth `wdir` - see Details.
#'   Default `FALSE`.
#'
#' @return A `SpatRaster` of coastal exposure values (0 = every sampled
#'   point upwind was land, i.e. fully sheltered from the coast; 1 = every
#'   sampled point upwind was sea, i.e. fully exposed to open water),
#'   covering the whole of `landsea`'s own extent.
#' @export
#'
#' @seealso [wind_shelter()], [wind_altitude()]
#'
#' @details
#' For each land cell, this samples land/sea coverage at a set of
#' distances along `wdir`, then takes the mean proportion of those samples
#' that were sea. There is
#' no explicit distance-based weight applied to that mean - the "nearby
#' coverage matters more" effect instead comes entirely from how the
#' sample *distances themselves* are spaced: they follow
#' \eqn{d_k = k^n \times \text{resolution}}, so raising `n` packs
#' relatively more samples into the near field (and fewer into the far
#' field) for the same total sample count, which pulls a plain average of
#' those samples towards whatever's nearby. `n = 2` reproduces the
#' original "inverse distance squared" calibration; `n = 1` gives evenly
#' (linearly) spaced samples instead, weighting the near and far field
#' more equally. This is a deliberate port of how the reference
#' implementation actually works, not a literal \eqn{1/\text{distance}^n}
#' term applied during averaging.
#'
#' The largest `k` used is always whatever's needed to reach `maxdist`
#' (see the `coarse` paragraph below for how `maxdist` itself is set) -
#' solved directly from \eqn{k^n \ge \text{maxdist} / \text{resolution}} -
#' so the search reaches all the way out to the edge of the available data
#' for every value of `n`, not just `n = 2`. Lower `n` needs more sample
#' points to cover the same distance evenly rather than in fast-widening
#' jumps, so a low `n` on a large search radius costs more to compute -
#' e.g. covering a 50 km radius at 100 m resolution takes on the order of
#' 170 samples at `n = 2` but around 4000 at `n = 1`. Because of this, when
#' `n < 2` the search distance is capped at twice `landsea`'s own diagonal
#' extent, regardless of how much further `coarse` would otherwise let it
#' reach - see the `coarse` paragraph below for why this matters in
#' practice and when you'll see a warning about it.
#'
#' When `coarse` is supplied, `landsea` and every mask in `coarse` are
#' first sorted by resolution (finest first) regardless of the order you
#' supply them in. For each sampled point, the finest level whose extent
#' actually contains that point is used; coarser levels are only
#' consulted when a point falls outside every finer level's extent. This
#' means a coarse mask only ever gets used for the far field, exactly
#' where a fine mask wouldn't reach anyway - you don't need to manually
#' work out a distance cutoff for each resolution. The maximum search
#' distance is extended to cover whichever supplied level has the
#' largest extent, using that level's true diagonal
#' (`sqrt(width^2 + height^2)`) - the distance needed to guarantee the
#' search reaches that level's edge from anywhere inside it, in any
#' direction. A very large `coarse` level (this package's bundled
#' [dtmc] dataset's second element, for instance, spans hundreds of
#' kilometres) can push that search distance far beyond what anyone
#' searching from a single `landsea`-sized region actually needs - and,
#' per the `n` paragraph above, the number of samples needed to reach it
#' evenly can get very large at low `n`. So at `n < 2` only, the search
#' distance is capped to twice `landsea`'s own diagonal extent, and a
#' `warning()` is issued if that cap actually reduced anything (i.e. if
#' `coarse` would otherwise have reached further). `n >= 2` is never
#' capped, since the sample count there grows mildly enough (roughly
#' `sqrt(maxdist / resolution)` at `n = 2`) that even a very large
#' `coarse` level stays cheap. If you genuinely need the search to reach
#' all the way out to a large `coarse` level's edge at a low `n`, raise
#' `n` toward 2, or supply a smaller `coarse` level instead.
#'
#' `wdir` is searched directly (azimuth = `wdir`, no adjustment) - the
#' same convention used by [horizon()] and [wind_shelter()].
#'
#' `jitter = TRUE` perturbs each sample point's azimuth independently (a
#' different random offset for every combination of output cell and
#' sample distance), rather than shifting a whole transect by one fixed
#' amount - this fans each transect out into a narrow +/-10-degree cone,
#' which smooths out grid-alignment/stair-step artifacts that can
#' otherwise appear at certain wind directions. The perturbation is drawn
#' from a fast deterministic hash of the cell's row/column and sample
#' index (not R's random number generator, which can't safely be called
#' from the parallel C++ loop this runs in) - so it is not affected by
#' `set.seed()`, but is fully reproducible: the same inputs always
#' produce the same output raster.
#'
#' `wdir = "all"` runs the whole search once per direction (with `jitter`
#' applied independently for each, if requested) and averages the 8
#' resulting rasters together, giving a single direction-agnostic
#' exposure value per cell.
#'
#' A sample point that falls outside every supplied level's extent
#' (`landsea`'s own, and any `coarse` level's) is simply skipped rather
#' than counted as sea or land, so a cell's value always comes from the
#' mean of whichever sample points *did* land somewhere with data - never
#' an error, and never a value dragged toward sea or land by missing data.
#' In practice this only matters near the edges of the data you supply:
#' a cell close to `landsea`'s own boundary (or, if `coarse` is supplied,
#' close to the boundary of whichever level covers that area) has fewer
#' far-field samples to draw on than one further in, so its value leans
#' more heavily on its near-field neighbourhood. If you want every cell in
#' some smaller region of interest to have the full complement of
#' far-field samples, supply a `landsea` (and `coarse`, if used) that
#' extends comfortably beyond that region, and crop the result to the
#' region you actually want afterwards.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' landsea <- rast(landsea100m)   # NA = sea, elevation/any value = land
#'
#' # Coastal exposure for a northeasterly and a westerly wind
#' ce_ne <- coastalexposure(landsea, wdir = 45)
#' ce_w  <- coastalexposure(landsea, wdir = 270)
#' plot(c(ce_ne, ce_w), main = c("Northeasterly", "Westerly"))
#'
#' # Extending the search with a coarser regional mask, so a small fine
#' # raster can still "see" coastline much further away
#' landsea_coarse <- aggregate(landsea, fact = 20, fun = "modal", na.rm = TRUE)
#' ce_far <- coastalexposure(landsea, wdir = 270, coarse = landsea_coarse)
#'
#' # A bundled package dataset works the same way, passed straight in --
#' # its elements are wrapped (terra::PackedSpatRaster) and get unwrapped
#' # automatically, no need to call rast()/unwrap() on them first
#' ce_regional <- coastalexposure(landsea, wdir = 270, coarse = dtmc)
#'
#' # A gentler near/far balance than the default n = 2
#' ce_linear <- coastalexposure(landsea, wdir = 270, n = 1)
#'
#' # Smooth out grid-alignment artifacts with a bit of angular jitter
#' ce_jit <- coastalexposure(landsea, wdir = 270, jitter = TRUE)
#'
#' # Omnidirectional exposure, averaged over 8 directions
#' ce_all <- coastalexposure(landsea, wdir = "all")
#' }
coastalexposure <- function(landsea, wdir, coarse = NULL, n = 2, jitter = FALSE) {
  stopifnot(inherits(landsea, "SpatRaster"))
  .check_projected(landsea, arg = "landsea")
  stopifnot(is.logical(jitter), length(jitter) == 1L, !is.na(jitter))
  stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n > 0)
  if (is.character(wdir)) {
    wdir <- match.arg(wdir, "all")
    dirs <- seq(0, 315, by = 45)
  } else {
    stopifnot(is.numeric(wdir), length(wdir) == 1L)
    dirs <- wdir
  }

  # A wrapped SpatRaster (terra::PackedSpatRaster) can't be used directly --
  # it's the on-disk/serialisable form terra::wrap() produces (and the form
  # bundled package datasets like dtmc are stored in, since a live
  # SpatRaster holds an external C++ pointer that can't survive save()/
  # load()). Unwrap any such element automatically so callers can pass a
  # wrapped dataset straight through without unwrapping it by hand first.
  unwrap_if_packed <- function(r) {
    if (inherits(r, "PackedSpatRaster")) terra::unwrap(r) else r
  }

  if (is.null(coarse)) {
    coarse <- list()
  } else if (inherits(coarse, "PackedSpatRaster")) {
    coarse <- list(terra::unwrap(coarse))
  } else if (inherits(coarse, "SpatRaster")) {
    coarse <- list(coarse)
  } else {
    stopifnot(is.list(coarse))
    coarse <- lapply(coarse, unwrap_if_packed)
    if (!all(vapply(coarse, inherits, logical(1), "SpatRaster")))
      stop("'coarse' must be a SpatRaster, a (possibly wrapped) list of ",
           "SpatRaster objects, or a wrapped SpatRaster.", call. = FALSE)
  }

  levels <- c(list(landsea), coarse)
  if (!all(vapply(levels, function(r) terra::same.crs(r, landsea), logical(1))))
    stop("'landsea' and every mask in 'coarse' must share the same CRS.",
         call. = FALSE)

  # Sort finest-to-coarsest so the C++ side can just try levels in order
  # and stop at the first one that covers a given sample point -- not
  # assuming the caller listed 'coarse' in any particular order.
  resos  <- vapply(levels, function(r) terra::res(r)[1], numeric(1))
  ord    <- order(resos)
  levels <- levels[ord]
  resos  <- resos[ord]

  to_binary <- function(r) {
    b <- r * 0 + 1
    b[is.na(b)] <- 0
    b
  }
  mats <- lapply(levels, function(r) as.matrix(to_binary(r), wide = TRUE))
  exts <- lapply(levels, terra::ext)

  reso <- resos[1]   # finest level's resolution -- sets near-field sample spacing

  # maxdist is the largest supplied level's true diagonal
  # (sqrt(width^2 + height^2)) -- the distance needed to guarantee a ray
  # starting anywhere inside that level's extent, travelling in any
  # direction, reaches all the way to that level's edge before running
  # out of samples. An earlier version of this formula
  # (sqrt(xres*width + yres*height)) grew only with the *square root* of
  # the raster's size rather than linearly with it, so it undershot the
  # true reach badly for anything but a tiny raster -- e.g. a 22 km x
  # 18 km, 100 m-resolution level got a maxdist of only 2000 m from that
  # formula, versus this level's actual ~28 km diagonal, silently
  # capping every search (regardless of `n`) at a small fraction of the
  # raster's real extent. Using the true diagonal fixes that for every
  # raster size and every `n`.
  maxdist_of <- function(r) {
    ex <- terra::ext(r)
    sqrt((ex$xmax - ex$xmin)^2 + (ex$ymax - ex$ymin)^2)
  }
  maxdist <- max(vapply(levels, maxdist_of, numeric(1)))

  # At n < 2, the sample count needed to reach 'maxdist' evenly grows fast
  # as 'maxdist' grows (kmax scales as (maxdist/reso)^(1/n), i.e. linearly
  # in maxdist itself at n = 1) -- fine when 'maxdist' is set by 'landsea'
  # alone, but a large 'coarse' level (e.g. this package's own bundled
  # 'dtmc' dataset, whose second element spans hundreds of km) can push
  # 'maxdist' -- and so the sample count -- far higher than anyone
  # searching from a single 'landsea'-sized region actually needs. Cap the
  # search distance at low n to twice landsea's own diagonal extent
  # (comfortably further than landsea itself, but nowhere near a whole
  # continental-scale coarse level) rather than letting it grow unchecked,
  # and warn when the cap actually reduces anything so it's never silent.
  # n >= 2 is left uncapped -- kmax's growth there is mild enough (e.g.
  # sqrt(maxdist/reso) at n = 2) that even a very large 'coarse' level
  # stays cheap.
  if (n < 2) {
    landsea_cap <- 2 * maxdist_of(landsea)
    if (maxdist > landsea_cap) {
      warning(sprintf(
        "'n' = %.3g is below 2, and the search distance implied by 'coarse' (%.0f m) would need a very large number of samples to reach evenly. Capping the search distance to %.0f m (2x 'landsea's own diagonal extent) instead -- increase 'n' toward 2, or supply a smaller 'coarse' level, to search further without this cap.",
        n, maxdist, landsea_cap), call. = FALSE)
      maxdist <- landsea_cap
    }
  }

  # Sample distances follow s = k^n * reso for a sequence of k running from
  # 0 up to whatever value reaches at least 'maxdist' -- solved directly
  # from k^n * reso >= maxdist, i.e. k >= (maxdist / reso)^(1/n) -- rather
  # than a fixed k range. A fixed range (this package's own original
  # calibration used k up to 125, in steps of 1/8) only reaches 'maxdist'
  # for n around 2 or higher; for n = 1 it silently capped the search at
  # 125 * reso regardless of how much further 'maxdist' actually allowed,
  # so a westward (or any direction) search could run out of samples well
  # short of the raster's edge and wrongly read cells as fully sheltered
  # (0) that do have sea further along the same line, just beyond that
  # fixed cap. Deriving the range from 'maxdist' itself removes that cap
  # entirely: the search always marches out to the true edge of the
  # available data, for any n. The eighth-unit step size is unchanged from
  # the original calibration; only the upper bound is now dynamic.
  kmax <- max(8, ceiling((maxdist / reso)^(1 / n) * 8))
  s <- (c(0, (8:kmax) / 8))^n * reso
  s <- s[s <= maxdist]

  # Output always covers the whole of 'landsea' itself -- no separate
  # target-region cropping. A sample point that falls outside every level's
  # extent is simply skipped by coastal_exposure_cpp() (see Details), so
  # cells near landsea's own edge still get a value, just from fewer
  # far-field samples than cells further in.
  lss <- to_binary(landsea)
  lsm <- as.matrix(lss, wide = TRUE)
  e   <- terra::ext(landsea)

  mask_xmin <- vapply(exts, function(x) x$xmin, numeric(1))
  mask_xmax <- vapply(exts, function(x) x$xmax, numeric(1))
  mask_ymin <- vapply(exts, function(x) x$ymin, numeric(1))
  mask_ymax <- vapply(exts, function(x) x$ymax, numeric(1))

  jitter_deg <- if (jitter) 10 else 0

  rasters <- lapply(dirs, function(d) {
    lsr <- coastal_exposure_cpp(lsm, reso, e$xmin, e$ymax, s, d,
                                 mats, resos,
                                 mask_xmin, mask_xmax, mask_ymin, mask_ymax,
                                 jitter_deg)
    .rast(lsr, lss)
  })

  if (length(dirs) > 1L) {
    out <- terra::mean(terra::rast(rasters))
    nm  <- "coastalexposure_all"
  } else {
    out <- rasters[[1]]
    nm  <- paste0("coastalexposure_", dirs)
  }

  # coastal_exposure_cpp()/rasters above carry the raw sampled land
  # proportion (0 = sea, 1 = land); flipped here to sea proportion (0 =
  # land/fully sheltered, 1 = sea/fully exposed) so higher values mean
  # more coastally exposed, matching the function's own name -- flipping
  # after averaging rather than before is equivalent (mean(1-x) ==
  # 1-mean(x)) but cheaper, only ever done once regardless of how many
  # directions were averaged. Renamed explicitly afterwards since
  # arithmetic on a SpatRaster isn't guaranteed to preserve layer names.
  out <- 1 - out
  names(out) <- nm
  out
}
