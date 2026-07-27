# Internal R-side utilities for terravars
# Not exported.

# Check that a SpatRaster has a projected (metre-unit) CRS. terra has no
# is.projected(); a raster is projected iff it is not geographic (lon/lat).
.check_projected <- function(r, arg = "dtm") {
  if (terra::is.lonlat(r)) {
    stop(sprintf(
      "'%s' must have a projected CRS with metre units (e.g. UTM, OSGB36).",
      arg
    ), call. = FALSE)
  }
}

# Coerce a NULL/NA/scalar/SpatRaster argument to an nr x nc matrix, filling
# with `default` when NULL or a length-1 NA. Shared by twostream() and
# skyview() for optional per-cell parameters (x, lref, ltra, gref, ...).
.to_mat <- function(v, default, nr, nc) {
  # SpatRaster is checked first: length() on a SpatRaster returns its layer
  # count (terra's method), not "is this a scalar" -- so a single-layer
  # raster would otherwise trip length(v) == 1L below and then hit is.na(v),
  # which returns a SpatRaster, not a scalar logical, and crashes the `&&`
  # with "invalid 'y' type in 'x && y'".
  if (inherits(v, "SpatRaster"))
    return(as.matrix(v, wide = TRUE))
  if (is.null(v) || (length(v) == 1L && is.na(v)))
    return(matrix(default, nr, nc))
  matrix(as.numeric(v)[1L], nr, nc)
}

# Check the documented `lref + ltra < 1` constraint (a leaf can't reflect
# and transmit more light than hits it) shared by twostream()/skyview()/
# canopy3d(). Violating it silently drives the two-stream algebra's `a`
# term negative, producing NaN (not an error) deep inside twostream_tme_cpp()
# / canopy_diffuse_frac_cpp() -- checked here, once resolved to matrices, so
# every caller gets the same clear error instead.
.check_leaf_optics <- function(lref_m, ltra_m) {
  bad <- !is.na(lref_m) & !is.na(ltra_m) & (lref_m + ltra_m >= 1)
  if (any(bad)) {
    stop(sum(bad), " cell(s) have 'lref' + 'ltra' >= 1 (a leaf can't ",
         "reflect and transmit more light than hits it); reduce 'lref' ",
         "and/or 'ltra' so their sum is < 1 everywhere.", call. = FALSE)
  }
}

# Check that two SpatRasters have the same extent and resolution
.check_compatible <- function(r1, r2, arg1 = "x", arg2 = "y") {
  if (!terra::compareGeom(r1, r2, stopOnError = FALSE)) {
    stop(sprintf(
      "'%s' and '%s' do not have the same extent, resolution, or CRS.",
      arg1, arg2
    ), call. = FALSE)
  }
}

# Extract per-cell latitude and longitude matrices from a SpatRaster.
# For projected CRS (e.g. OSGB36, UTM) the native x/y coordinates are
# reprojected to WGS84 geographic (EPSG:4326) via terra::project().
# For rasters already in geographic CRS the native coordinates are
# returned directly.
# Returns a list with elements $lats and $lons, each an nrow x ncol matrix.
.latslonsfromr <- function(r) {
  e  <- terra::ext(r)
  rs <- terra::res(r)
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)

  # Cell-centre y coordinates (top row first, matching R matrix row order)
  lats_vec <- rep(seq(e$ymax - rs[2] / 2,
                      e$ymin + rs[2] / 2,
                      length.out = nr), nc)

  # Cell-centre x coordinates (sorted so they group by column, column-major)
  lons_vec <- rep(seq(e$xmin + rs[1] / 2,
                      e$xmax - rs[1] / 2,
                      length.out = nc), nr)
  lons_vec <- lons_vec[order(lons_vec)]

  if (terra::is.lonlat(r)) {
    # Already geographic: native coords are lon/lat
    lons <- array(lons_vec, dim = c(nr, nc))
    lats <- array(lats_vec, dim = c(nr, nc))
  } else {
    # Projected: reproject cell centres to WGS84
    xy <- cbind(x = lons_vec, y = lats_vec)
    ll <- terra::project(xy,
                         from = terra::crs(r),
                         to   = "EPSG:4326")
    lons <- array(ll[, 1], dim = c(nr, nc))
    lats <- array(ll[, 2], dim = c(nr, nc))
  }

  list(lats = lats, lons = lons)
}

# Plain (non-slope-corrected) direct-beam canopy extinction coefficient,
# K_be(theta, x) -- the ellipsoidal leaf-angle-distribution formula of
# Campbell & Norman (1998), a function of solar zenith angle and the leaf
# angle distribution parameter `x` only (no dependence on surface slope or
# on the solar index). Mirrors kBeCpp() in src/utils.h exactly.
#   zen_deg - solar zenith angle (degrees), scalar
#   xv_m    - leaf angle distribution parameter, matrix: the ratio of mean
#             leaf-area projection onto a horizontal plane to that onto a
#             vertical plane, following the Campbell & Norman (1998)
#             ellipsoidal convention (1 = spherical/isotropic, 0 = vertical
#             leaves, Inf = horizontal leaves)
# Returns a matrix of K_be values, same dims as xv_m (NA where xv_m is NA).
.k_be_r <- function(zen_deg, xv_m) {
  zen_r <- zen_deg * pi / 180
  cos_z <- cos(zen_r)
  sin_z <- sin(zen_r)
  tan_z <- if (cos_z > 1e-10) sin_z / cos_z else 1e10

  k <- matrix(NA_real_, nrow(xv_m), ncol(xv_m))

  # `valid` excludes NA/NaN cells up front so the branch masks below are
  # always plain TRUE/FALSE (never NA) -- without this, `!is.finite(NA)`
  # is TRUE, which would misclassify NA cells as "horizontal" (x = Inf)
  # instead of leaving them as NA in the output.
  valid    <- !is.na(xv_m)
  is_sph   <- valid & xv_m == 1
  is_vertl <- valid & xv_m == 0            # x = 0: vertical leaves
  is_horzl <- valid & !is.finite(xv_m)     # x = Inf: horizontal leaves (NA excluded by `valid`)
  is_gen   <- valid & !is_sph & !is_vertl & !is_horzl

  k[is_sph] <- if (cos_z > 1e-10) 0.5 / cos_z else 6000
  # Vertical leaves, azimuthally random: G(zenith) = (2/pi) * sin(zenith)
  # (mean of |cos(phi)| over a uniform leaf azimuth is 2/pi), so
  # K_be = G / cos(zenith) = (2/pi) * tan(zenith) -- NOT tan(zenith) alone,
  # which is the projection for a single leaf whose azimuth is locked to
  # face the sun exactly, not the azimuthally-averaged canopy value assumed
  # everywhere else here (G as a function of zenith only). Mirrors
  # kBeCpp()'s xv == 0.0 branch in src/utils.h exactly.
  k[is_vertl] <- (2 / pi) * tan_z
  k[is_horzl] <- 1.0
  xg <- xv_m[is_gen]
  k[is_gen]   <- sqrt(xg^2 + tan_z^2) / (xg + 1.774 * (xg + 1.182)^(-0.733))

  k
}

# Slope-corrected direct-beam canopy extinction coefficient (vectorised,
# R-side). Mirrors canopy_kd() in src/twostream.cpp exactly, applied here
# as a Beer-Lambert attenuation layer on top of an already-computed
# geometric solar index (see solarindex()'s `pai`/`x` arguments).
#   zen_deg - solar zenith angle (degrees), scalar
#   xv_m    - leaf angle distribution parameter, matrix (see .k_be_r())
#   si_m    - solar index on the sloped surface (cos of incidence angle,
#             clamped to [0, 1]), matrix, same dims as xv_m
# Returns a matrix of extinction coefficients, capped at 6000 (near-zero
# transmission) to avoid Inf when si_m is ~0.
#
# Note: no explicit si_m-NA masking is needed here beyond what `ifelse()`
# already does below -- when si_m is NA, `ifelse(si_m > 1e-10, ...)` itself
# evaluates to NA regardless of k's value, so NA propagates correctly either
# way.
.canopy_extinction <- function(zen_deg, xv_m, si_m) {
  zen_r <- zen_deg * pi / 180
  cos_z <- cos(zen_r)

  k  <- .k_be_r(zen_deg, xv_m)
  kd <- ifelse(si_m > 1e-10, k * cos_z / si_m, 6000)
  pmin(kd, 6000)
}

# Canopy interception rate G(theta, x) = K_be(theta, x) * cos(theta) -- the
# Ross-Nilson G-function: the mean projection of unit leaf area onto the
# plane perpendicular to the beam, as a fraction of that leaf area.
# Independent of surface slope (unlike .canopy_extinction()'s slope-
# corrected kd) -- it characterises the canopy/beam geometry alone. Used as
# the limiting (PAI -> Inf) canopy-interception term in solarindex()'s
# `ground = FALSE` option (see solarindex() Details). Clamped to [0, 1].
#   zen_deg - solar zenith angle (degrees), scalar
#   xv_m    - leaf angle distribution parameter, matrix (see .k_be_r())
# Returns a matrix, same dims as xv_m.
.canopy_intercept <- function(zen_deg, xv_m) {
  zen_r <- zen_deg * pi / 180
  cos_z <- cos(zen_r)
  g <- .k_be_r(zen_deg, xv_m) * cos_z
  pmin(pmax(g, 0), 1)
}

# Coerce a SpatRaster/PackedSpatRaster (or already-plain matrix/array) to a
# plain matrix (single layer) or array (multi-layer). Used by basindelin()
# ahead of its C++ call, which expects a plain matrix, not a SpatRaster.
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r <- terra::rast(r)
  if (class(r)[1] == "SpatRaster") {
    if (dim(r)[3] == 1) {
      m <- as.matrix(r, wide = TRUE)
    } else m <- as.array(r)
  } else {
    m <- r
  }
  return(m)
}

# Wrap a plain matrix/array back into a SpatRaster, taking its extent and
# CRS from a template raster `tem`. Used by basindelin() to turn the
# basin-ID matrix returned by basinCpp() back into a SpatRaster.
.rast <- function(m, tem) {
  r <- terra::rast(m)
  terra::ext(r) <- terra::ext(tem)
  terra::crs(r) <- terra::crs(tem)
  r
}

# Pad a plain nr x nc matrix by 1 cell in every direction, filling the new
# border with `fill`. Shared by basindelin()/basinmerge()/fillsinks()/
# flowacc(), which each build a padded matrix before their C++ worker call
# -- the 1-cell border doubles as "every real cell on the raster's own
# edge already touches an outlet" for those workers, with no special-
# casing needed for the raster boundary itself.
#   m           - the matrix to pad
#   fill        - value for the new border (e.g. 9999 for an elevation
#                 matrix's NA/no-data sentinel, NA for a basin-ID/weight
#                 companion array)
#   na_as_fill  - if TRUE, interior NA cells in `m` are also replaced with
#                 `fill` before padding (used for the elevation matrix
#                 itself, whose NA/no-data cells share the same 9999
#                 sentinel as the border); FALSE (default) leaves `m`'s own
#                 NA pattern untouched and only pads the new border (used
#                 for basin-ID/weight companion arrays, whose own NA cells
#                 mean something different and must not be touched here).
.pad1 <- function(m, fill, na_as_fill = FALSE) {
  nr <- dim(m)[1]
  nc <- dim(m)[2]
  if (na_as_fill) m[is.na(m)] <- fill
  out <- array(fill, dim = c(nr + 2, nc + 2))
  out[2:(nr + 1), 2:(nc + 1)] <- m
  out
}

# Single-pass nearest-neighbour infill of NA cells in a plain matrix `m`
# (row/col convention matching .is()/.rast()). Any NA cell touching at
# least one non-NA cell among its 8 neighbours is replaced by that
# neighbour's value; cardinal neighbours (distance 1) are checked before
# diagonal ones (distance sqrt(2)), and within a tie the first match in
# scan order wins, so the result is deterministic. An NA cell with no
# non-NA neighbour at all (i.e. two or more cells from any real data) is
# left NA -- deliberately not an iterative flood-fill, since a focal
# calculation with a 3x3 window (terra::terrain()'s slope/aspect) can
# never look further than one cell out anyway, so filling further would
# never change any real cell's result.
.nn_fill <- function(m) {
  nr <- nrow(m); nc <- ncol(m)
  pad <- matrix(NA_real_, nr + 2, nc + 2)
  pad[2:(nr + 1), 2:(nc + 1)] <- m

  filled <- m
  need   <- is.na(m)

  offsets <- list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1),   # cardinal (dist 1)
                  c(-1, -1), c(-1, 1), c(1, -1), c(1, 1)) # diagonal (dist sqrt(2))

  for (o in offsets) {
    if (!any(need)) break
    nb <- pad[(2:(nr + 1)) + o[1], (2:(nc + 1)) + o[2]]
    ok <- need & !is.na(nb)
    filled[ok] <- nb[ok]
    need[ok]   <- FALSE
  }
  filled
}

# Nearest-neighbour-filled version of `dtm`, used by solarindex()/
# twostream()/twi()'s `fillna` argument in place of a flat elevation-0
# substitution. terra::terrain()'s slope/aspect calculation needs a full
# 3x3 neighbourhood, so it comes out NA not just over genuine NA cells
# (sea, data gaps) but at every cell touching one (an NA "halo"), AND at
# the outermost ring of `dtm` itself, since cells there have neighbours
# that fall outside the raster entirely -- both are really the same
# problem (a focal window missing a real value), so both are fixed the
# same way here. Returns a list:
#   $extended -- `dtm` extended by 1 cell in every direction (using only
#     dtm's own data, no larger dataset required) and nearest-neighbour
#     filled, for terra::terrain() -- crop back to `dtm`'s own extent
#     afterwards.
#   $interior -- the same fill, without the extra border ring, on
#     `dtm`'s own extent -- for anywhere else `dtm`'s elevation values
#     feed a calculation directly (e.g. flowacc()/ray-marching), so
#     those see the same fill values as the slope/aspect calculation at
#     every interior cell.
.fillna_dtm <- function(dtm) {
  ext <- terra::extend(dtm, 1)
  filled_m <- .nn_fill(.is(ext))
  ext_filled <- .rast(filled_m, ext)
  list(
    extended = ext_filled,
    interior = terra::crop(ext_filled, dtm)
  )
}

# ---------------------------------------------------------------------------
# Cost estimate + warning for full raster-wide ray-marched sky-direction
# scans (skyview()'s hgt/Goudriaan-3D diffuse path, canopy3d()'s diffuse
# fraction). Both re-examine the *entire* raster once per requested sky
# direction (ndir * n_elev directions in total), so total cost scales with
# both the direction count AND the raster's cell count -- a warning keyed
# on direction count alone can't tell a small DTM (fast regardless of how
# many directions) from a large one (slow even at modest direction counts),
# and was firing on the package's own shipped defaults (ndir=32, n_elev=6)
# on an ordinary-sized DTM before this fix.
#
# .RAY_MARCH_CELLS_PER_SEC is calibrated from one real benchmark reported by
# the package author: 192 directions (the defaults above) across a 187 x 214
# DTM took ~2 seconds in practice, i.e. 187*214*192 / 2 ~= 3.8e6
# cell-directions/second -- halved here as a safety margin, since per-ray
# cost also depends on canopy height/structure (a ray stops early once it
# climbs above the tallest canopy in the scene -- see veg_path_diffuse_cpp()
# -- so patchier or taller canopy means more steps before that early exit)
# and on how many CPU cores are actually available to the OpenMP-parallel
# C++ loop, neither of which one benchmark on one machine can capture. This
# is a heuristic estimate, not a guarantee -- recalibrate if it proves
# consistently mis-calibrated in practice (e.g. on much larger DTMs, or
# machines without OpenMP support).
.RAY_MARCH_CELLS_PER_SEC <- 1.9e6
.RAY_MARCH_WARN_SEC      <- 15

# Format a seconds estimate as a short, human-readable string.
.format_est_time <- function(sec) {
  if (sec < 60)         sprintf("%.0f seconds", sec)
  else if (sec < 3600)  sprintf("%.1f minutes", sec / 60)
  else                  sprintf("%.1f hours", sec / 3600)
}

# n_dir_total - total sky directions requested (ndir * n_elev)
# nr, nc      - raster dimensions the ray-march runs across
# label       - short description of which calculation this is, used in the
#               warning message (e.g. "skyview()'s hgt diffuse attenuation")
.warn_ray_march_cost <- function(n_dir_total, nr, nc, label) {
  est_sec <- (n_dir_total * nr * nc) / .RAY_MARCH_CELLS_PER_SEC
  if (est_sec > .RAY_MARCH_WARN_SEC) {
    warning(sprintf(
      "%s: %d total ray-marched sky directions (ndir x n_elev) across a %d x %d raster -- estimated around %s, based on raster size as well as direction count (a rough estimate, not a guarantee). Consider reducing 'ndir'/'n_elev', or working with a smaller extent.",
      label, n_dir_total, nr, nc, .format_est_time(est_sec)
    ), call. = FALSE, immediate. = TRUE)
  }
}
