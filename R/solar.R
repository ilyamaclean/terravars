# ============================================================
# solar.R
#
# Solar / radiation functions:
#   - solarindex()   direct beam solar index on a sloped surface
#   - horizon()      terrain horizon elevation angle
#   - skyview()      fraction of sky hemisphere visible from each cell
#   - twostream()    below-canopy radiation via two-stream approximation
# ============================================================

#' Solar index
#'
#' Computes how directly the sun's rays strike a sloped surface at each point
#' in the landscape - a slope facing the sun scores higher than one facing
#' away from it, and flat ground facing straight up scores highest of all
#' when the sun is overhead. The result is a plain multiplier between 0 (the
#' sun's rays are grazing along the surface, contributing no direct light)
#' and 1 (the rays hit the surface square-on): multiply it by the direct
#' sunlight intensity you'd measure on a flat surface to get the intensity on
#' the sloped surface instead. (Technically, this is the cosine of the angle
#' between the sun's rays and the surface's normal, clamped to \[0, 1\].)
#' Optional add-ons let you also account for the sun being blocked by nearby
#' hills (`shade`) or by an overlying plant canopy (`pai`/`x`/`hgt`).
#'
#' @details
#' Two usage modes are determined by `tme`:
#'
#' **Annual mode** (`tme = "annual"`): gives a typical year-round solar
#' index without you having to supply a specific list of dates and times.
#' Internally, sun positions are sampled across the full year 2023 and
#' grouped into a manageable number of representative positions, so the
#' result is accurate without the cost of checking all 8760 hours one by
#' one. `wgt` still defaults to `"mean"` here, the same as elsewhere,
#' weighting every daylight position equally - but `wgt = "clearsky"` is
#' usually the better choice for a representative annual estimate, since it
#' weights sunnier-looking sun positions more heavily rather than treating
#' persistently overcast ones as just as representative. It isn't applied
#' automatically and needs to be requested explicitly.
#'
#' **Time-series mode** (`tme` is a set of specific dates/times, in UTC):
#' the solar index is worked out separately for each timestep you supply,
#' then combined into a single result according to `wgt` (or kept separate
#' per timestep if `wgt = "all"`). By default (`cenloc = TRUE`) the sun's
#' position is calculated once for the middle of the map and used
#' everywhere - see the `cenloc` argument below for when you'd want
#' per-cell positions instead.
#'
#' **Weighting** (`wgt`): controls how multiple timesteps get combined into
#' one result. `"mean"` (default) is a simple average across all the
#' timesteps you supplied, counting any night-time ones as zero - so a set
#' of daytime-only times gives a daytime average, while a full 24-hour
#' sequence gives a lower day-and-night average. `"clearsky"` instead
#' leans on timesteps with strong, clear-sky sunlight, often a better match
#' for how much solar energy a spot actually receives over typical weather.
#' `"all"` skips averaging and returns one map per timestep (not available
#' in annual mode).
#'
#' **Terrain shading** (`shade = TRUE`): also switches the solar index to
#' zero at any point currently hidden from the sun behind surrounding hills
#' or ridges, on top of the plain slope/aspect effect described above. This
#' means checking the horizon in a number of directions, pre-computed once
#' and then reused for every sun position - see [horizon()] for how the
#' speed/accuracy tradeoff (`nstep`) works. `ndir` is an upper limit on how
#' many directions get checked, not a fixed cost: with `cenloc = TRUE` and
#' a specific set of times in `tme` (i.e. not annual mode), the exact sun
#' azimuths needed are already known in advance, so only that many
#' directions are actually searched - one direction for a single timestep,
#' not a full `ndir`-direction sweep - falling back to the full `ndir`
#' sweep only once there are more distinct daytime azimuths than `ndir`
#' itself. Annual mode and `cenloc = FALSE` always use the full sweep,
#' since neither has a fixed, known-in-advance set of exact azimuths to
#' target instead.
#'
#' **Vegetation attenuation** (`pai`): supplying a plant area index raster
#' - a measure of how much leaf and branch material sits above each point,
#' see the `pai` argument below - reduces the solar index to account for
#' some of that direct sunlight being blocked by the canopy, similarly to
#' how light dims passing through a hedge or forest edge. The higher the
#' sun sits above the horizon, the shorter its path through the canopy and
#' the less it's blocked; that angle-dependence, together with how the
#' canopy's leaves are typically oriented (`x`), uses the same underlying
#' calculation as [twostream()], so the two functions should give closely
#' matching results for the same inputs.
#'
#' **Diagonal vegetation paths** (`hgt`): the `pai`-only option above
#' assumes sunlight travels straight down through whatever canopy sits
#' directly above each point - a reasonable shortcut when the canopy is
#' fairly uniform. Supplying canopy height (`hgt`) instead traces the sun's
#' actual path sideways across neighbouring points' canopies too, so a beam
#' slipping in over a gap or a shorter patch of vegetation next door is
#' blocked less than the straight-down shortcut would suggest, while one
#' arriving through a taller, denser neighbour is blocked more. This gives
#' a more realistic result but is noticeably slower to compute (it can't be
#' pre-computed once and reused the way terrain shading can), currently
#' requires `cenloc = TRUE`, and prints a warning if you request it for
#' more than 100 sun positions at once.
#'
#' **Combined canopy and ground interception** (`ground = FALSE`): by
#' default, the result describes sunlight that has passed *through* the
#' canopy and reached the ground beneath it. Setting `ground = FALSE`
#' instead describes the sunlight caught by the canopy and ground
#' *together*, treating them as one combined surface - useful when what
#' you actually care about is the vegetation surface itself (e.g. a grass
#' sward), not what's happening on the soil underneath it. Only applies
#' when `pai` is also supplied. When `shade = TRUE`, a cell hidden from the
#' sun by surrounding terrain returns 0 here just as it does in the
#' `ground = TRUE` case: if the direct beam is blocked by terrain it reaches
#' neither the ground nor the canopy above it, so both parts of the combined
#' surface are dark.
#'
#' **Plant area index convention** (`pai_horiz`): a supplied `pai` can be
#' quoted per unit horizontal ground area (the usual reporting convention,
#' `pai_horiz = TRUE`, default) or per unit sloped ground-surface area
#' (`pai_horiz = FALSE`). The two are identical on flat ground and differ by
#' `cos(slope)` on a slope; the conversion is applied internally so the same
#' physical canopy attenuates the beam identically under either convention.
#' This is purely an area-units/optical-path-length adjustment - leaf angles
#' are assumed fixed relative to the horizontal (gravity) and do not rotate
#' with the ground, regardless of `pai_horiz`. It affects only the simpler
#' pai-only path; the `hgt` ray-marching path always uses `pai` as supplied,
#' since it derives the through-canopy path length from the terrain surface
#' directly.
#'
#' **NA cells in `dtm`** (`fillna = TRUE`): a DTM with `NA` cells - most
#' commonly the sea, at a coastline - causes a knock-on problem beyond just
#' those cells themselves: slope and aspect are worked out from a small
#' neighbourhood around each cell, so any land cell touching the sea also
#' comes out `NA`, even though it has a perfectly good elevation of its
#' own. The same thing happens, for the same reason, at the outermost ring
#' of `dtm` itself, since those cells' neighbourhoods reach past the edge
#' of the raster entirely. `fillna = TRUE` works around both by filling
#' each affected cell from its nearest real neighbour before the
#' calculation (not a flat substitute value - a coastal cell picks up
#' whatever elevation the nearest land actually has, and an edge cell
#' picks up whatever its nearest interior neighbour has), then masking the
#' result back to `dtm`'s real `NA` pattern afterwards - so the sea itself
#' stays `NA` in the output, but the coastline right next to it, and the
#' raster's own edge, no longer do. Extending `dtm` with a larger buffer
#' raster isn't needed for this - the edge fill only ever uses `dtm`'s own
#' data.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must have a projected CRS with metre units.
#' @param tme A set of dates/times (`POSIXct`/`POSIXlt`, **in UTC**) to
#'   compute the solar index for, or the character string `"annual"` for a
#'   typical year-round estimate covering many sun positions automatically
#'   (see Details).
#' @param wgt Character string controlling how multiple timesteps are
#'   combined into one result: `"mean"` (default) - simple average,
#'   counting night-time as zero; `"clearsky"` - weighted toward
#'   clear-sky-bright timesteps; `"all"` - no averaging, one map per
#'   timestep (not available with `tme = "annual"`). See Details.
#' @param slp Optional `SpatRaster` of slope (degrees). Worked out from
#'   `dtm` automatically if left as `NULL`.
#' @param asp Optional `SpatRaster` of aspect/orientation (degrees, N = 0,
#'   clockwise). Worked out from `dtm` automatically if left as `NULL`.
#' @param shade Logical. Set to `TRUE` to also zero out points currently
#'   hidden from the sun by surrounding hills or ridges (see Details).
#'   Default `FALSE` (only the point's own slope/aspect is considered).
#' @param ndir Integer. Upper limit on the number of directions used to
#'   check the surrounding horizon when `shade = TRUE`. Default 32 - higher
#'   gives a more precise horizon at extra computation cost. Only actually
#'   used in full when there are that many (or more) distinct sun positions
#'   to place directions at; a single timestep, for instance, only ever
#'   needs its own direction, not a full sweep (see Details).
#' @param nstep Integer. Passed to [horizon()] to control how thoroughly the
#'   horizon is searched. `0` (default) checks carefully all the way to the
#'   edge of the map; a positive integer trades some accuracy for speed by
#'   checking more sparsely further away. See [horizon()] for details.
#' @param cenloc Logical. Controls whether the sun's position is worked out
#'   once for the whole map (fast, and fine for most maps) or separately for
#'   every single grid cell (slower, but more accurate for maps so large
#'   that the sun's position genuinely differs noticeably from one side to
#'   the other). `TRUE` (default) uses one shared position; a warning
#'   appears if the map spans more than roughly 1 degree of latitude or
#'   longitude, as a hint that `FALSE` may be worth the extra time. Ignored
#'   in annual mode, which always uses one shared position.
#' @param pai Optional `SpatRaster` of plant area index (m^2 m^-2) - a
#'   measure of how much leaf, stem and branch material sits above each
#'   point (higher values mean a denser canopy that blocks more light).
#'   Supplying this reduces the solar index to account for that shading -
#'   see Details. `NULL` (default) applies no vegetation shading at all.
#' @param x Leaf angle distribution parameter, used together with `pai`.
#'   This controls how the canopy's leaves are typically angled, which
#'   affects how much light passing through gets blocked versus slipping
#'   through gaps: leave it at the default (1) for a "typical" mixed
#'   canopy if you're not sure, use a value close to 0 for canopies made up
#'   mostly of upright/vertical leaves (e.g. many grasses and some
#'   conifers), or a large value (or `Inf`) for canopies of mostly flat,
#'   horizontal leaves (e.g. many broadleaf trees). (Technically, this is
#'   the Campbell & Norman (1998) ellipsoidal leaf-angle-distribution
#'   parameter - the ratio of leaf area projected onto a horizontal plane
#'   versus a vertical one.) May be a single value or a `SpatRaster` if it
#'   varies across the map. Ignored when `pai` is `NULL`.
#' @param hgt Optional `SpatRaster` of canopy height above the local ground
#'   (metres). Supplying this alongside `pai` switches on the more
#'   realistic diagonal calculation described in Details, instead of
#'   assuming sunlight travels straight down through each point's own
#'   canopy column. Requires `cenloc = TRUE`. Cells with `pai > 0` but
#'   `hgt <= 0` (or missing) are treated as having no canopy there, with a
#'   warning. `NULL` (default) uses the simpler, faster straight-down
#'   approximation.
#' @param ground Logical. Only relevant when `pai` is supplied. `TRUE`
#'   (default) returns sunlight that reaches the ground beneath the
#'   canopy. `FALSE` instead returns the sunlight caught by the canopy and
#'   ground combined, treating them as one surface - useful for short
#'   vegetation like grass where the canopy itself is what you're
#'   interested in (see Details).
#' @param pai_horiz Logical. Only relevant when `pai` is supplied (and only
#'   affects the simpler pai-only path, not the `hgt` ray-marching path).
#'   Declares which area the supplied `pai` is measured against. `TRUE`
#'   (default) - `pai` is per unit **horizontal** ground area, the usual
#'   way leaf/plant area index is reported. `FALSE` - `pai` is per unit
#'   **sloped** ground-surface area. On flat ground the two are identical;
#'   on a slope they differ by a factor of `cos(slope)`, which is applied
#'   internally so that a fixed physical amount of canopy attenuates the
#'   beam by the same amount whichever convention you use. Leaf angles are
#'   unaffected either way (see Details). `pai_horiz` changes the numerical
#'   result on sloping ground, so set it to match how your `pai` data was
#'   actually derived.
#' @param fillna Logical. If `TRUE`, `NA` cells in `dtm` (e.g. sea, or any
#'   other no-data gap), and the outer ring of `dtm` itself, are
#'   temporarily filled from their nearest real neighbour for the purposes
#'   of the calculation, and the result is masked back to `dtm`'s own `NA`
#'   pattern afterwards. Fixes a specific problem at the edge of `NA`
#'   regions and at the raster's own edge: slope/aspect worked out from
#'   `dtm` would otherwise come out `NA` not just over the `NA` cells
#'   themselves but at every valid land cell that touches one - e.g. an
#'   entire coastline's worth of land cells picking up `NA` purely for
#'   being next to the sea - and at every cell along `dtm`'s own boundary,
#'   which has no neighbouring cell at all just past the edge. Default
#'   `FALSE` (previous behaviour, `NA` cells left as `NA` and their
#'   neighbours affected as described). See Details.
#'
#' @return A single-layer `SpatRaster` of solar index values
#'   (dimensionless, 0-1), or a multi-layer `SpatRaster` when `wgt = "all"`.
#'   Layer name is `"solar_index"` in time-series mode and
#'   `"solar_index_annual"` in annual mode.
#' @export
#'
#' @seealso [horizon()], [skyview()], [twostream()], [canopy3d()]
#'
#' @references
#' Bird, R.E. & Hulstrom, R.L. (1981). A simplified clear sky model for
#' direct and diffuse insolation on horizontal surfaces. \emph{SERI/TR-642-761}.
#'
#' Kasten, F. (1965). A new table and approximate formula for relative
#' optical air mass. \emph{Archives for Meteorology, Geophysics and
#' Bioclimatology, Series B}, 14, 206--223.
#'
#' Campbell, G.S. & Norman, J.M. (1998). \emph{An Introduction to
#' Environmental Biophysics}, 2nd ed. Springer, New York. (canopy
#' extinction coefficient / leaf angle distribution used by `pai`/`x`).
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#'
#' # --- Mean over a sequence of daytime hours (summer solstice) ---
#' tme <- seq(as.POSIXct("2023-06-21 05:00", tz = "UTC"),
#'            as.POSIXct("2023-06-21 20:00", tz = "UTC"),
#'            by = "1 hour")
#' si <- solarindex(dtm, tme)
#' plot(si)
#'
#' # --- Clear-sky-weighted annual mean ---
#' si_ann <- solarindex(dtm, "annual", wgt = "clearsky")
#' plot(si_ann)
#'
#' # --- Annual mean with terrain shading ---
#' si_ann_s <- solarindex(dtm, "annual", shade = TRUE)
#' plot(si_ann_s)
#'
#' # --- Per-timestep stack ---
#' si_stack <- solarindex(dtm, tme, wgt = "all")
#' plot(si_stack)
#'
#' # --- Per-cell solar positions (large DTM spanning > 1 degree) ---
#' si_percell <- solarindex(dtm, tme, cenloc = FALSE)
#'
#' # --- Attenuated by a canopy (plant area index and leaf angle, Lizard, Cornwall) ---
#' pai  <- rast(pai)
#' vegx <- rast(vegx)
#' si_veg <- solarindex(dtm, tme, pai = pai, x = vegx)
#' plot(si_veg)
#'
#' # --- Diagonal ray-marched attenuation (canopy height varies spatially) ---
#' hgt <- rast(hgt)
#' si_veg3d <- solarindex(dtm, tme, pai = pai, x = vegx, hgt = hgt)
#' plot(si_veg3d)
#'
#' # --- Combined canopy+ground interception (e.g. short-sward energy balance) ---
#' si_canopy <- solarindex(dtm, tme, pai = pai, x = vegx, ground = FALSE)
#' plot(si_canopy)
#'
#' # --- DTM with NA cells at a coastline: avoid an NA halo on land ---
#' si_coast <- solarindex(dtm, tme, fillna = TRUE)
#' plot(si_coast)
#' }
solarindex <- function(dtm,
                       tme,
                       wgt    = "mean",
                       slp    = NULL,
                       asp    = NULL,
                       shade  = FALSE,
                       ndir   = 32L,
                       nstep  = 0L,
                       cenloc = TRUE,
                       pai    = NULL,
                       x      = NULL,
                       hgt    = NULL,
                       ground = TRUE,
                       pai_horiz = TRUE,
                       fillna = FALSE) {

  stopifnot(inherits(dtm, "SpatRaster"))
  .check_projected(dtm)
  stopifnot(is.logical(ground), length(ground) == 1L, !is.na(ground))
  stopifnot(is.logical(pai_horiz), length(pai_horiz) == 1L, !is.na(pai_horiz))
  wgt <- match.arg(wgt, c("mean", "clearsky", "all"))

  # ---- Optional NA fill (see `fillna` below) ----
  # NA cells in `dtm` (sea, or any other no-data gap) make terra::terrain()
  # return NA slope/aspect not just at the NA cells themselves but at every
  # valid land cell that touches one, since its focal window includes an NA
  # neighbour -- e.g. every coastal land cell next to the sea -- and, for
  # the same reason, at the outermost ring of `dtm` itself, since those
  # cells' focal windows reach outside the raster entirely. When
  # fillna = TRUE both are fixed via .fillna_dtm(): NA cells (and, for the
  # slope/aspect calculation only, a 1-cell border extending past `dtm`'s
  # own edge) are temporarily filled from their nearest real neighbour, and
  # the result is masked back to `dtm`'s own NA pattern at every return
  # point below, so genuinely-NA cells stay NA in the output while land
  # cells next to them (or at the raster's edge) get a real value instead
  # of NA.
  dtm_na <- dtm
  fill_ext <- NULL
  if (isTRUE(fillna)) {
    fill_ext <- .fillna_dtm(dtm)
    dtm <- fill_ext$interior
  }
  .mask_fillna <- function(r) if (isTRUE(fillna)) terra::mask(r, dtm_na) else r

  # ---- Resolve annual mode ----
  annual_mode <- is.character(tme) && length(tme) == 1L &&
                 tolower(tme) == "annual"
  if (annual_mode) {
    if (wgt == "all")
      stop("wgt = \"all\" is not supported in annual mode - this would produce ~8760 layers.",
           call. = FALSE)
    cenloc <- TRUE   # annual always uses centre position
  } else if (!inherits(tme, c("POSIXct", "POSIXlt"))) {
    stop("'tme' must be a POSIXct/POSIXlt vector or the character string \"annual\".",
         call. = FALSE)
  }

  # ---- Slope / aspect ----
  # When fillna = TRUE, computed from fill_ext$extended (1 cell wider than
  # `dtm` in every direction, see .fillna_dtm()) so the outer ring gets a
  # real value too, then cropped back to `dtm`'s own extent.
  if (is.null(slp)) {
    slp <- if (isTRUE(fillna)) {
      terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm_na)
    } else terra::terrain(dtm, "slope", unit = "degrees")
  }
  if (is.null(asp)) {
    asp <- if (isTRUE(fillna)) {
      terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm_na)
    } else terra::terrain(dtm, "aspect", unit = "degrees")
  }
  sl_m <- as.matrix(slp, wide = TRUE)
  as_m <- as.matrix(asp, wide = TRUE)
  nr   <- terra::nrow(dtm)
  nc   <- terra::ncol(dtm)

  # ---- Vegetation attenuation (optional) ----
  use_veg <- !is.null(pai)
  if (use_veg) {
    stopifnot(inherits(pai, "SpatRaster"))
    pai_m <- as.matrix(pai, wide = TRUE)
    x_m   <- if (inherits(x, "SpatRaster")) {
      # Checked before the is.null()/is.na() branch below: length() on a
      # SpatRaster returns its layer count (terra's method), not 1 meaning
      # "scalar" -- so a single-layer raster like this would otherwise trip
      # length(x) == 1L and then hit is.na(x), which returns a SpatRaster,
      # not a scalar logical, and crashes `&&` ("invalid 'y' type").
      as.matrix(x, wide = TRUE)
    } else if (is.null(x) || (length(x) == 1L && is.na(x))) {
      matrix(1.0, nr, nc)
    } else {
      matrix(as.numeric(x)[1L], nr, nc)
    }

    # ---- Plant area index area convention (pai_horiz) ----
    # The per-pixel-column canopy extinction below (.canopy_extinction /
    # .veg_transmission's non-hgt path) is written for `pai` expressed per
    # unit *sloped* ground-surface area: its cos(zen)/si factor is the
    # optical path-length elongation of a beam through a slope-parallel
    # canopy layer, which is naturally quoted per sloped area. When the user
    # supplies pai per unit *horizontal* area instead (pai_horiz = TRUE, the
    # usual LAI/PAI convention), convert to per-sloped area first: a fixed
    # amount of leaf spread over a tilted patch covers more sloped surface
    # than its horizontal footprint by 1/cos(slope), so its per-sloped
    # density is lower -- hence multiply by cos(slope). Leaf *angles* are
    # unaffected (they stay gravity-referenced via G(zen, x)); this is purely
    # the area-units / path-length term. The hgt (diagonal ray-marching) path
    # derives its path length from the DTM surface geometry directly, so it
    # takes pai as supplied and this conversion is not applied there.
    pai_col <- pai_m
    if (isTRUE(pai_horiz)) pai_col <- pai_m * cos(sl_m * pi / 180)
  }

  # ---- 3D vegetation ray-marching (optional; requires pai + cenloc) ----
  # A beam travelling diagonally across neighbouring pixels' canopies,
  # not just straight down through the target pixel's own column --
  # see veg_path_atten_cpp() for the ray-marching/vertically-uniform-
  # column approximation this relies on.
  use_hgt <- !is.null(hgt)
  if (use_hgt) {
    if (!use_veg)
      stop("'hgt' requires 'pai' to also be supplied.", call. = FALSE)
    if (!cenloc)
      stop("'hgt' (diagonal-path vegetation attenuation) requires cenloc = TRUE: ",
           "it ray-marches a single shared solar position across the whole raster ",
           "per timestep, which is incompatible with cenloc = FALSE's per-cell ",
           "solar positions.", call. = FALSE)
    if (isTRUE(pai_horiz))
      warning("pai_horiz = TRUE has no effect on the hgt (diagonal ray-marching) ",
              "path: that path derives the through-canopy optical path length from ",
              "the DTM surface geometry directly, so pai is used as supplied. The ",
              "horizontal-vs-sloped-area distinction only affects the simpler ",
              "pai-only (no hgt) path.", call. = FALSE)
    stopifnot(inherits(hgt, "SpatRaster"))
    hgt_m <- as.matrix(hgt, wide = TRUE)

    bad_hgt <- !is.na(pai_m) & pai_m > 0 & (is.na(hgt_m) | hgt_m <= 0)
    if (any(bad_hgt))
      warning(sum(bad_hgt), " cell(s) have pai > 0 but hgt <= 0 (or NA); vegetation ",
              "attenuation is treated as zero (no canopy) at those cells.",
              call. = FALSE)

    elev_m <- as.matrix(dtm, wide = TRUE)
    res    <- terra::res(dtm)[1]
  }

  # Local helper: vegetation transmission fraction (tau = fraction of the
  # geometric solar index `mat` that reaches the ground), via diagonal
  # ray-marching (hgt) or the simpler per-pixel-column approximation (pai
  # alone). Requires use_veg = TRUE.
  .veg_transmission <- function(mat, zen, azi) {
    if (use_hgt) {
      veg_path_atten_cpp(elev_m, hgt_m, pai_m, x_m, res, zen, azi, nstep)
    } else {
      kd <- .canopy_extinction(zen, x_m, mat)
      exp(-kd * pai_col)
    }
  }

  # Local helper: apply vegetation attenuation to a per-solar-position
  # matrix. When `ground = TRUE` (default), returns the fraction of `mat`
  # transmitted to the ground, mat * tau -- unchanged from the pre-`ground`
  # behaviour. When `ground = FALSE`, instead returns the canopy+ground
  # combined interception, (1 - tau) * G(zen, x) + tau * mat, where
  # G(zen, x) is the canopy-interception limit as PAI -> Inf (see
  # .canopy_intercept()) -- relevant for short-sward/grass canopies where
  # canopy and ground jointly intercept the beam (see Details). Used
  # wherever vegetation attenuation can apply, including the cenloc = FALSE
  # branch below (which can never have use_hgt = TRUE, guarded above).
  .apply_veg <- function(mat, zen, azi) {
    if (!use_veg) return(mat)
    tau <- .veg_transmission(mat, zen, azi)
    if (ground) {
      mat * tau
    } else {
      g <- .canopy_intercept(zen, x_m)
      # Terrain shading must zero the canopy term as well as the ground term.
      # `mat` (from solar_index_cpp) is already 0 wherever the cell is hidden
      # from the sun by surrounding terrain, but `g` -- the canopy
      # interception limit -- is a function of solar zenith and leaf angle
      # only and knows nothing about the local horizon, so without this it
      # would keep its full sunlit value in terrain shade and (1 - tau) * g
      # would wrongly light up shaded valley floors under any canopy. When
      # the beam is blocked by terrain it reaches neither the ground nor the
      # canopy above it, so both terms must vanish. (The sun-below-flat-
      # horizon / night case is already handled: there cos(zen) <= 0, so
      # .canopy_intercept()'s own clamp drives g to 0.) The mask is rebuilt
      # from the same horizon slice passed to solar_index_cpp, so it matches
      # exactly which cells that call zeroed for shade.
      if (use_shade) {
        hz <- hor_for_az(azi)
        if (length(hz) > 0L) g[(90 - zen) <= hz] <- 0
      }
      (1 - tau) * g + tau * mat
    }
  }

  # ---- Geographic coordinates ----
  ll   <- .latslonsfromr(dtm)
  cr   <- (nr + 1L) %/% 2L
  cc   <- (nc + 1L) %/% 2L
  lat0 <- ll$lats[cr, cc]
  lon0 <- ll$lons[cr, cc]

  if (cenloc) {
    lat_span <- diff(range(ll$lats, na.rm = TRUE))
    lon_span <- diff(range(ll$lons, na.rm = TRUE))
    if (lat_span > 1 || lon_span > 1)
      warning("DTM spans more than 1 degree; solar position is computed only at the ",
              "raster centre (cenloc = TRUE), which may introduce errors near the ",
              "grid edges. Set cenloc = FALSE to compute per-cell positions.",
              call. = FALSE)
  }

  # ---- Time-series sun positions, resolved early ----
  # Needed here rather than down in TIME-SERIES MODE below, so the shade
  # direction-count optimisation in the horizon stack (next) can see the
  # actual sun azimuths before deciding how many directions to search.
  # No-op in annual mode, which works out its own representative positions
  # later and never uses `tme`/`years`/etc. directly.
  if (!annual_mode) {
    if (!inherits(tme, "POSIXlt")) tme <- as.POSIXlt(tme, tz = "UTC")
    years  <- as.integer(tme$year + 1900L)
    months <- as.integer(tme$mon  + 1L)
    days   <- as.integer(tme$mday)
    hours  <- tme$hour + tme$min / 60.0 + tme$sec / 3600.0
    n_t    <- length(years)

    if (cenloc) {
      pos_mat <- solar_pos_cpp(lat0, lon0, years, months, days, hours)
      zen_v   <- pos_mat[, 1L]   # solar zenith (degrees)
      azi_v   <- pos_mat[, 2L]   # solar azimuth (degrees)
    }
  }

  # ---- Horizon stack (built once if shade = TRUE) ----
  # Searches `ndir` evenly-spaced directions by default. But when the
  # exact sun azimuths needed are already known (cenloc = TRUE time-series
  # mode, above) and there are no more of them than `ndir`, search exactly
  # those instead -- a single timestep only ever needs its own azimuth,
  # not a full ndir-direction sweep. Falls back to the original
  # evenly-spaced sweep whenever there are more distinct daytime azimuths
  # than `ndir` (no worse than before), and always for annual mode (whose
  # representative-position count isn't known until after this point) and
  # cenloc = FALSE (whose per-cell azimuths vary continuously across the
  # raster and are interpolated in C++ from an evenly-spaced cube, which
  # needs a regular grid of directions to interpolate between).
  use_shade <- isTRUE(shade)
  exact_az  <- NULL
  if (use_shade && !annual_mode && cenloc) {
    az <- sort(unique(round(azi_v[zen_v < 90.0], 4)))
    if (length(az) > 0L && length(az) <= ndir) exact_az <- az
  }

  if (use_shade) {
    hor_azimuths <- if (!is.null(exact_az)) exact_az
                     else seq(0, 360 - 360 / ndir, length.out = ndir)
    hor_r    <- horizon(dtm, azi = hor_azimuths, nstep = nstep)
    n_az_hor <- terra::nlyr(hor_r)
    hor_cube <- array(NA_real_, dim = c(nr, nc, n_az_hor))
    for (k in seq_len(n_az_hor))
      hor_cube[, , k] <- as.matrix(hor_r[[k]], wide = TRUE)
  } else {
    hor_azimuths <- numeric(0)
    hor_cube     <- array(0.0, dim = c(0L, 0L, 0L))
  }
  az_step <- 360.0 / ndir   # used by the cenloc = FALSE C++ interpolation path

  # Local helper: return the horizon slice nearest to a given azimuth --
  # circular nearest-match against whatever directions were actually
  # searched above, correct whether that's the original evenly-spaced
  # ndir sweep or the reduced exact-match set. Returns an empty matrix
  # when shade = FALSE (signals no shading to C++).
  hor_for_az <- function(azi_deg) {
    if (!use_shade) return(matrix(nrow = 0, ncol = 0))
    d <- abs(((azi_deg - hor_azimuths + 180) %% 360) - 180)
    hor_cube[, , which.min(d)]
  }

  # =================================================================
  # ANNUAL MODE: binning over representative sun positions
  # =================================================================
  if (annual_mode) {

    raw <- annual_solar_weights_cpp(lat0, lon0, 2023L, 15.0, 60.0, 1013.25, 1.0)
    if (nrow(raw) == 0L)
      stop("No daytime sun positions found; check DTM location.", call. = FALSE)

    # Bin weights: clearsky irradiance or uniform count
    w_raw <- if (wgt == "clearsky") raw[, 3L] else rep(1.0, nrow(raw))

    # Bin into 5-degree zenith x 5-degree azimuth cells
    bsz     <- 5.0
    bin_key <- floor(raw[, 1L] / bsz) * 1000L + floor(raw[, 2L] / bsz)
    keys    <- unique(bin_key)
    n_b     <- length(keys)
    bw      <- bzen <- bazi <- numeric(n_b)
    for (i in seq_len(n_b)) {
      sel     <- bin_key == keys[i]
      ws      <- w_raw[sel]
      bw[i]   <- sum(ws)
      bzen[i] <- stats::weighted.mean(raw[sel, 1L], ws)
      bazi[i] <- stats::weighted.mean(raw[sel, 2L], ws)
    }

    # Retain only the most influential bins (fixed internal thresholds)
    keep <- bw >= 0.005 * max(bw)
    bw   <- bw[keep]; bzen <- bzen[keep]; bazi <- bazi[keep]
    if (length(bw) > 40L) {
      ord  <- order(bw, decreasing = TRUE)[seq_len(40L)]
      bw   <- bw[ord]; bzen <- bzen[ord]; bazi <- bazi[ord]
    }

    if (use_hgt && length(bw) > 100L)
      warning(length(bw), " solar positions requested with vegetation ray-marching ",
              "(pai + hgt); this is considerably more expensive per position than ",
              "the pai-only or hgt-free paths.", call. = FALSE, immediate. = TRUE)

    acc  <- matrix(0.0, nr, nc)
    wtot <- 0.0
    for (i in seq_along(bw)) {
      mat <- solar_index_cpp(sl_m, as_m, bzen[i], bazi[i], hor_for_az(bazi[i]))
      mat <- .apply_veg(mat, bzen[i], bazi[i])
      acc  <- acc + mat * bw[i]
      wtot <- wtot + bw[i]
    }

    result <- if (wtot > 0) acc / wtot else matrix(0.0, nr, nc)
    out <- terra::rast(dtm)
    terra::values(out) <- result
    names(out) <- "solar_index_annual"
    return(.mask_fillna(out))
  }

  # =================================================================
  # TIME-SERIES MODE
  # =================================================================
  # tme/years/months/days/hours/n_t already resolved above, ahead of the
  # horizon stack.

  if (use_hgt && n_t > 100L)
    warning(n_t, " solar positions requested with vegetation ray-marching (pai + ",
            "hgt); this is considerably more expensive per position than the ",
            "pai-only or hgt-free paths. Consider reducing the number of ",
            "timesteps, or omit `hgt` (falls back to the faster per-pixel-column ",
            "approximation) for large runs.", call. = FALSE, immediate. = TRUE)

  # -----------------------------------------------------------------
  # cenloc = TRUE: single solar position at DTM centre for all cells
  # -----------------------------------------------------------------
  if (cenloc) {

    # pos_mat/zen_v/azi_v already computed above, ahead of the horizon
    # stack.
    wt_v <- if (wgt == "clearsky") clearsky_cpp(zen_v) else rep(1.0, n_t)

    if (wgt == "all") {
      layers <- vector("list", n_t)
      for (i in seq_len(n_t)) {
        mat <- if (zen_v[i] >= 90.0) {
          matrix(0.0, nr, nc)
        } else {
          mat_i <- solar_index_cpp(sl_m, as_m, zen_v[i], azi_v[i], hor_for_az(azi_v[i]))
          .apply_veg(mat_i, zen_v[i], azi_v[i])
        }
        r <- terra::rast(dtm)
        terra::values(r) <- mat
        names(r) <- format(tme[i], "%Y%m%d_%H%M")
        layers[[i]] <- r
      }
      return(.mask_fillna(terra::rast(layers)))
    }

    acc  <- matrix(0.0, nr, nc)
    wtot <- 0.0
    for (i in seq_len(n_t)) {
      if (zen_v[i] >= 90.0 || wt_v[i] <= 0.0) next
      mat  <- solar_index_cpp(sl_m, as_m, zen_v[i], azi_v[i],
                               hor_for_az(azi_v[i]))
      mat  <- .apply_veg(mat, zen_v[i], azi_v[i])
      acc  <- acc + mat * wt_v[i]
      wtot <- wtot + wt_v[i]
    }

    # "mean": divide by n_t so night-time steps count as zero
    # (consistent with solar_index_tme_cpp convention)
    result <- if (wgt == "mean") acc / n_t
              else if (wtot > 0)  acc / wtot
              else                 matrix(0.0, nr, nc)

  # -----------------------------------------------------------------
  # cenloc = FALSE: per-cell solar position from lat/lon matrices
  # -----------------------------------------------------------------
  } else {

    lat_m <- ll$lats
    lon_m <- ll$lons

    if (wgt == "mean" && !use_veg) {
      # Delegate entirely to the optimised C++ function
      result <- solar_index_tme_cpp(sl_m, as_m, lat_m, lon_m,
                                     years, months, days, hours,
                                     hor_cube, az_step = az_step)
    } else {
      # Centre-location solar positions - used for clearsky weights, and
      # (when pai is supplied) as an approximate per-timestep zenith angle
      # for vegetation attenuation. The solar zenith barely varies across a
      # typical DTM at any given instant, so this is a good approximation
      # even when per-cell positions are used for the geometric index.
      pos_cen <- solar_pos_cpp(lat0, lon0, years, months, days, hours)
      wt_v    <- if (wgt == "clearsky") clearsky_cpp(pos_cen[, 1L])
                 else rep(1.0, n_t)

      if (wgt == "all") {
        layers <- vector("list", n_t)
        for (i in seq_len(n_t)) {
          mat <- solar_index_tme_cpp(sl_m, as_m, lat_m, lon_m,
                                      years[i], months[i], days[i], hours[i],
                                      hor_cube, az_step = az_step)
          mat <- .apply_veg(mat, pos_cen[i, 1L], pos_cen[i, 2L])
          r <- terra::rast(dtm)
          terra::values(r) <- mat
          names(r) <- format(tme[i], "%Y%m%d_%H%M")
          layers[[i]] <- r
        }
        return(.mask_fillna(terra::rast(layers)))
      }

      acc  <- matrix(0.0, nr, nc)
      wtot <- 0.0
      for (i in seq_len(n_t)) {
        if (wt_v[i] <= 0.0) next
        mat  <- solar_index_tme_cpp(sl_m, as_m, lat_m, lon_m,
                                     years[i], months[i], days[i], hours[i],
                                     hor_cube, az_step = az_step)
        mat  <- .apply_veg(mat, pos_cen[i, 1L], pos_cen[i, 2L])
        acc  <- acc + mat * wt_v[i]
        wtot <- wtot + wt_v[i]
      }
      # "mean": divide by n_t, consistent with the delegated fast path
      # above (night-time / zero-weight steps count as zero); "clearsky":
      # divide by the summed clear-sky weight.
      result <- if (wgt == "mean") acc / n_t
                else if (wtot > 0)  acc / wtot
                else                 matrix(0.0, nr, nc)
    }
  }

  out <- terra::rast(dtm)
  terra::values(out) <- result
  names(out) <- "solar_index"
  .mask_fillna(out)
}


#' Terrain horizon angle
#'
#' For each point in the landscape, works out how high up the sky the
#' surrounding terrain reaches in one or more compass directions - in
#' effect, "if I stood here and looked that way, how far up would I have to
#' tilt my head before I stopped seeing hills and started seeing open sky?"
#' A high horizon angle means nearby terrain blocks a lot of the view in
#' that direction (e.g. deep in a valley); an angle near 0 means the view
#' in that direction is basically flat and open. This underpins several
#' other functions in this package - terrain shading of sunlight
#' ([solarindex()]), sky-view factor ([skyview()]), and wind shelter
#' ([wind_shelter()]) - and can also be used on its own.
#'
#' @details
#' Two search modes are available, controlled by `nstep`:
#'
#' **Careful mode** (`nstep = 0`, default): checks every point along the
#' way, all the way out to the edge of the map. This is the most accurate
#' option and is recommended for complex terrain - urban areas, narrow
#' valleys, or anywhere the blocking terrain might sit at an irregular
#' distance away.
#'
#' **Fast mode** (`nstep > 0`): checks progressively fewer, more widely
#' spaced points further away (at 1^2, 2^2, ..., `nstep`^2 cell widths - e.g.
#' with `nstep = 15` and 100 m resolution, the checked distances are 100,
#' 400, 900, ..., 22,500 m). This is substantially faster and is a
#' reasonable shortcut because closer terrain generally blocks more of the
#' view than distant terrain of the same height, so it's more likely to be
#' what actually sets the horizon - but it can miss blocking terrain that
#' sits at one of the skipped intermediate distances. Best suited to open
#' landscapes where the horizon is typically set by distant topography
#' rather than something close by.
#'
#' This function is used internally by [solarindex()], [skyview()], and
#' [wind_shelter()]. Call it directly when you want to work the horizon out
#' once and reuse it across several calls to those functions, rather than
#' having each one recompute it separately.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must have a projected CRS with metre units.
#' @param azi Numeric vector. The compass direction(s) to check (degrees,
#'   N = 0, clockwise). One raster layer is returned per direction.
#'   Defaults to the eight cardinal and intercardinal directions (N, NE, E,
#'   ...).
#' @param nstep Integer. `0` (default) checks carefully all the way to the
#'   edge of the map. A positive integer switches to the faster, sparser
#'   search described in Details, checking that many points (maximum
#'   search distance = `nstep`^2 x resolution, in metres).
#' @param obs_hgt Numeric. Height above ground (metres) of the observer
#'   doing the looking, added to each point's own elevation before
#'   searching. `0` (default) matches the usual "standing on the ground"
#'   horizon. Raising this lets an elevated viewpoint - e.g. a wind sensor
#'   or turbine hub mounted a few metres up - "see over" more of the
#'   surrounding terrain than a ground-level point would, lowering the
#'   horizon angle it finds. Distinct from `hgt` elsewhere in this package
#'   (canopy height): this is the height of the *observer*, not of any
#'   vegetation.
#'
#' @return A `SpatRaster` with one layer per azimuth, containing terrain
#'   horizon elevation angles in degrees. Layer names are `horizon_<azimuth>`.
#' @export
#'
#' @seealso [solarindex()], [skyview()], [wind_shelter()]
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#'
#' # Eight directions (default)
#' ha <- horizon(dtm)
#' plot(ha)
#'
#' # 36 directions for sky-view factor - linear mode (accurate)
#' ha36 <- horizon(dtm, azi = seq(0, 350, by = 10))
#'
#' # Faster quadratic mode for open landscape
#' ha_q <- horizon(dtm, azi = seq(0, 315, by = 45), nstep = 15)
#'
#' # Pre-compute once and reuse across functions
#' ha32 <- horizon(dtm, azi = seq(0, 348.75, length.out = 32))
#' svf  <- skyview(hor = ha32)
#'
#' # Horizon as seen from 10 m up (e.g. a mast-mounted sensor), which sees
#' # over more of the surrounding terrain than a ground-level point
#' ha_10m <- horizon(dtm, obs_hgt = 10)
#' }
horizon <- function(dtm,
                    azi     = seq(0, 315, by = 45),
                    nstep   = 0L,
                    obs_hgt = 0) {
  stopifnot(inherits(dtm, "SpatRaster"))
  .check_projected(dtm)

  res  <- terra::res(dtm)[1]
  elev <- as.matrix(dtm, wide = TRUE)

  layers <- lapply(azi, function(az) {
    mat <- horizon_angle_cpp(elev, res, az, as.integer(nstep), obs_hgt)
    r   <- terra::rast(dtm)
    terra::values(r) <- mat
    names(r) <- paste0("horizon_", az)
    r
  })

  terra::rast(layers)
}


#' Sky-view factor
#'
#' Computes the sky-view factor (SVF) for each point in the landscape -
#' roughly, how much of the overall sky is visible from that point, as a
#' fraction between 0 (completely enclosed, e.g. the bottom of a deep gorge)
#' and 1 (a full, unobstructed view of the sky, e.g. the middle of a flat
#' open field). This matters because it controls how much diffuse sky light
#' - light scattered around by the atmosphere and clouds, as opposed to the
#' sun's direct beam - reaches that point: a valley floor surrounded by
#' hills receives less diffuse light than an exposed hilltop, even under
#' identical overcast skies. With no vegetation, this is purely about the
#' shape of the surrounding terrain - how much sky is blocked by nearby
#' hills and ridges. Optionally supplying `pai` also accounts for an
#' overlying plant canopy blocking part of the sky as well (see Details).
#'
#' @details
#' **Terrain-only** (`pai = NULL`, the default): purely a function of the
#' surrounding terrain shape, checked in `ndir` directions around each
#' point.
#'
#' **Vegetation attenuation** (`pai` supplied, no `hgt`): reduces the
#' terrain-only sky-view factor to also account for an overlying canopy
#' blocking part of the sky. This uses one shared "how much sky light gets
#' through the canopy" fraction for the whole point, which doesn't change
#' with the sun's position and - for both calculation options below -
#' doesn't vary by direction either, since both assume the canopy stretches
#' out uniformly in every direction (see the `hgt` option below to instead
#' account for real local variation in canopy shape). The combined result
#' is simply the terrain-only sky-view factor multiplied by this canopy
#' fraction. Two calculation options are available via `model`:
#' \describe{
#'   \item{`"twostream"`}{(default) the more detailed and accurate option,
#'     using the same underlying physics as [twostream()]'s `fdifdown`
#'     output - which also accounts for light bouncing between leaves and
#'     the ground before it eventually escapes or reaches the ground
#'     (rather than assuming it's simply absorbed on first contact). Needs
#'     a few extra inputs describing how reflective the leaves and ground
#'     are (`x`, `lref`, `ltra`, `gref`).}
#'   \item{`"goudriaan"`}{a simpler, faster approximation (Jan Goudriaan's
#'     simplified version of the two-stream model) that treats the canopy
#'     more like a plain filter: light passing through gets progressively
#'     dimmer the more leaf material it passes, with no bouncing-around
#'     effects included. Less accurate than `"twostream"`, and only needs
#'     `lref`/`ltra` (`x` and `gref` are ignored) - its simpler maths is
#'     also what makes the full 3D option below possible.}
#' }
#'
#' **3D vegetation attenuation** (`pai` and `hgt` both supplied, with
#' `model = "goudriaan"`): rather than assuming the canopy is a uniform
#' blanket in every direction, this instead looks at the actual shape of
#' the terrain and canopy in the surrounding neighbourhood - so a gap in
#' the canopy to the north, or a taller stand of trees to the south,
#' genuinely changes the result in that direction rather than being
#' averaged away. It works by checking many individual directions across
#' the sky (`ndir` compass directions x `n_elev` bands of sky height) and
#' combining what's actually visible in each into one overall sky-view
#' factor that already folds in both terrain and canopy blocking together.
#' This is noticeably slower than the options above, since it re-examines
#' the whole map once for every direction checked (`ndir * n_elev` checks
#' in total) - a warning appears if the estimated runtime (based on both
#' the direction count and the raster's size, since a big map is slower to
#' re-examine than a small one) looks like it's heading past around 15
#' seconds; the defaults (`ndir = 32`, `n_elev = 6`) run quickly on an
#' ordinary-sized DTM. It's only
#' available with `model = "goudriaan"`: the more detailed two-stream
#' option treats the canopy as an infinite uniform layer with no concept of
#' "this direction" versus "that direction", so it can't be evaluated this
#' way - `hgt` is ignored, with a warning, whenever
#' `model = "twostream"`.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must have a projected CRS with metre units. Required unless
#'   `hor` is supplied (and always required for the `hgt` 3D option, since
#'   the actual elevation surface is needed, not just pre-computed horizon
#'   angles).
#' @param ndir Integer. Number of compass directions checked around each
#'   point when computing the terrain horizon. Default 32 - higher gives a
#'   more precise result at greater computational cost. Also used as the
#'   number of directions checked for the `hgt` 3D option (see Details).
#' @param nstep Integer. Passed to [horizon()], and to the `hgt` 3D option
#'   if used, to control how thoroughly the surroundings are searched. `0`
#'   (default) checks carefully all the way to the edge of the map; a
#'   positive integer trades some accuracy for speed. See [horizon()] for
#'   details.
#' @param hor Optional multi-layer `SpatRaster` of horizon angles already
#'   computed by [horizon()]. Supplying this skips recomputing them from
#'   `dtm`/`ndir`/`nstep` - useful when you've already computed a horizon
#'   stack for another purpose and want to reuse it. `dtm` may then be left
#'   out entirely, unless the `hgt` 3D option is also being used.
#' @param pai Optional `SpatRaster` of plant area index (m^2 m^-2) - a
#'   measure of how much leaf and branch material sits above each point
#'   (higher values mean a denser canopy that blocks more sky). Supplying
#'   this reduces the sky-view factor to account for that shading - see
#'   Details. `NULL` (default) applies no vegetation shading.
#' @param x Leaf angle distribution parameter, used with `pai` only when
#'   `model = "twostream"`. Controls how the canopy's leaves are typically
#'   angled: leave it at the default (1) for a "typical" mixed canopy if
#'   unsure, use a value close to 0 for canopies of mostly upright/vertical
#'   leaves (e.g. many grasses and some conifers), or a large value (or
#'   `Inf`) for canopies of mostly flat, horizontal leaves (e.g. many
#'   broadleaf trees). Follows the same Campbell & Norman (1998)
#'   ellipsoidal convention as [solarindex()] and [twostream()]. May be a
#'   scalar or `SpatRaster`. Ignored when `pai` is `NULL` or
#'   `model = "goudriaan"`.
#' @param lref Leaf reflectance (0--1) - the fraction of light hitting a
#'   leaf that bounces straight back off it. `NULL` or `NA` uses 0.4
#'   (matching [twostream()]'s default). May be a scalar or `SpatRaster`.
#'   Ignored when `pai` is `NULL`.
#' @param ltra Leaf transmittance (0--1) - the fraction of light hitting a
#'   leaf that passes straight through it. `NULL` or `NA` uses 0.2
#'   (matching [twostream()]'s default). Note `lref + ltra` must be < 1
#'   (a leaf can't reflect and transmit more light than hits it). May be a
#'   scalar or `SpatRaster`. Ignored when `pai` is `NULL`.
#' @param gref Ground reflectance (0--1) - the fraction of light hitting
#'   the bare ground beneath the canopy that bounces back upward (higher
#'   for pale, dry, or snow-covered ground; lower for dark, wet soil).
#'   Default 0.15 (matching [twostream()]'s default). Used only when
#'   `model = "twostream"`. May be a scalar or `SpatRaster`. Ignored when
#'   `pai` is `NULL` or `model = "goudriaan"`.
#' @param hgt Optional `SpatRaster` of vegetation height (metres), matching
#'   `dtm`'s grid. Supplying this together with `pai` and
#'   `model = "goudriaan"` switches on the more realistic 3D calculation
#'   described in Details, which looks at the actual local shape of the
#'   canopy instead of assuming a uniform blanket. Ignored (with a warning)
#'   otherwise - including whenever `model = "twostream"`, which has no way
#'   to use it.
#' @param model Either `"twostream"` (default, more accurate) or
#'   `"goudriaan"` (simpler and faster, and the only option that supports
#'   `hgt`'s 3D calculation); selects which canopy calculation is used when
#'   `pai` is supplied. See Details.
#' @param n_elev Integer. Number of bands of sky height checked (from the
#'   horizon up to directly overhead) for the `hgt` 3D option only. Default
#'   6 - higher gives a more precise result at greater computational cost
#'   (total checks scale as `ndir * n_elev`).
#'
#' @return A `SpatRaster` of sky-view factor values (dimensionless, 0--1).
#' @export
#'
#' @seealso [horizon()], [solarindex()], [twostream()], [canopy3d()]
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' svf <- skyview(dtm)
#' plot(svf)
#'
#' # Pre-compute horizon angles and pass in directly
#' ha  <- horizon(dtm, azi = seq(0, 348.75, length.out = 32))
#' svf <- skyview(hor = ha)
#'
#' # Attenuated by a canopy, full two-stream model (Lizard, Cornwall)
#' pai_r  <- rast(pai)
#' vegx_r <- rast(vegx)
#' svf_veg <- skyview(dtm, pai = pai_r, x = vegx_r)
#' plot(svf_veg)
#'
#' # Same, using Goudriaan's simplified Beer-Lambert approximation instead
#' svf_goud <- skyview(dtm, pai = pai_r, model = "goudriaan")
#'
#' # Full 3D ray-marched diffuse attenuation using canopy height
#' hgt_r <- rast(hgt)
#' svf_3d <- skyview(dtm, pai = pai_r, hgt = hgt_r, model = "goudriaan",
#'                    n_elev = 6)
#' plot(svf_3d)
#' }
skyview <- function(dtm    = NULL,
                    ndir   = 32L,
                    nstep  = 0L,
                    hor    = NULL,
                    pai    = NULL,
                    x      = NULL,
                    lref   = NULL,
                    ltra   = NULL,
                    gref   = 0.15,
                    hgt    = NULL,
                    model  = "twostream",
                    n_elev = 6L) {

  model <- match.arg(model, c("twostream", "goudriaan"))

  if (!is.null(dtm)) {
    stopifnot(inherits(dtm, "SpatRaster"))
    .check_projected(dtm)
  }

  if (is.null(hor)) {
    if (is.null(dtm))
      stop("'dtm' must be supplied when 'hor' is not.", call. = FALSE)
    azimuths <- seq(0, 360 - 360 / ndir, length.out = ndir)
    hor      <- horizon(dtm, azi = azimuths, nstep = nstep)
  }

  n_az  <- terra::nlyr(hor)
  nrows <- terra::nrow(hor)
  ncols <- terra::ncol(hor)
  # Azimuths matching hor's own layers, regardless of whether hor was built
  # above or supplied pre-computed -- assumes the standard evenly-spaced
  # convention used by horizon()/skyview() (documented requirement for the
  # hgt ray-marched path below).
  azimuths <- seq(0, 360 - 360 / n_az, length.out = n_az)

  cube <- array(NA_real_, dim = c(nrows, ncols, n_az))
  for (k in seq_len(n_az)) {
    cube[, , k] <- as.matrix(hor[[k]], wide = TRUE)
  }

  result <- sky_view_factor_cpp(cube)

  # ---- Vegetation attenuation of diffuse radiation (optional) ----
  use_veg      <- !is.null(pai)
  use_hgt_diff <- use_veg && !is.null(hgt) && model == "goudriaan"

  if (!is.null(hgt) && !use_hgt_diff) {
    warning("'hgt' is only used when 'pai' is supplied together with ",
            "model = \"goudriaan\" (the two-stream model has no 3D ",
            "generalisation to diagonal ray-marched paths); 'hgt' is ",
            "ignored here.", call. = FALSE)
  }

  if (use_veg) {
    stopifnot(inherits(pai, "SpatRaster"))
    pai_m  <- as.matrix(pai, wide = TRUE)
    lref_m <- .to_mat(lref, 0.4, nrows, ncols)
    ltra_m <- .to_mat(ltra, 0.2, nrows, ncols)
    .check_leaf_optics(lref_m, ltra_m)

    if (use_hgt_diff) {

      # ---- 3D diffuse ray-marching (Goudriaan, hgt supplied) ----
      if (is.null(dtm))
        stop("'dtm' must be supplied for the diagonal (hgt) diffuse ",
             "calculation (model = \"goudriaan\" with hgt): elevation data ",
             "is needed for ray-marching, not just pre-computed horizon ",
             "angles.", call. = FALSE)
      stopifnot(inherits(hgt, "SpatRaster"))
      hgt_m  <- as.matrix(hgt, wide = TRUE)
      elev_m <- as.matrix(dtm, wide = TRUE)
      res    <- terra::res(dtm)[1]

      bad_hgt <- !is.na(pai_m) & pai_m > 0 & (is.na(hgt_m) | hgt_m <= 0)
      if (any(bad_hgt))
        warning(sum(bad_hgt), " cell(s) have pai > 0 but hgt <= 0 (or NA); ",
                "vegetation attenuation is treated as zero (no canopy) at ",
                "those cells.", call. = FALSE)

      # Cost estimate accounts for raster size as well as direction count --
      # see .warn_ray_march_cost() in utils.R for the calibration/rationale.
      n_dir_total <- n_az * n_elev
      .warn_ray_march_cost(n_dir_total, nrows, ncols,
                            "skyview()'s hgt diffuse attenuation")

      # Elevation-angle bin edges/midpoints spanning the hemisphere (0-90
      # degrees); each bin's exact (no-canopy) contribution is the closed-
      # form sin(e)cos(e) antiderivative sin^2(e), so summing bins telescopes
      # to exactly cos^2(h) when canopy transmission is 1 everywhere --
      # matching sky_view_factor_cpp()'s formula exactly in the no-canopy
      # limit, regardless of n_elev.
      edges <- seq(0, 90, length.out = n_elev + 1L)
      mids  <- (edges[-1L] + edges[-(n_elev + 1L)]) / 2

      acc <- matrix(0.0, nrows, ncols)
      for (j in seq_len(n_az)) {
        h_j <- cube[, , j]   # this azimuth's terrain horizon angle (degrees)
        for (i in seq_len(n_elev)) {
          tau_ij <- veg_path_diffuse_cpp(elev_m, hgt_m, pai_m, lref_m, ltra_m,
                                          res, 90 - mids[i], azimuths[j], nstep)
          # Bin's visible fraction: 0 if entirely below this cell's horizon,
          # the full closed-form weight if entirely above, a partial
          # closed-form weight if the horizon falls inside the bin.
          lo_eff <- pmax(edges[i], h_j)
          w_ij   <- pmax(0, sin(edges[i + 1L] * pi / 180)^2 -
                            sin(lo_eff        * pi / 180)^2)
          acc    <- acc + tau_ij * w_ij
        }
      }
      # This fully replaces (not multiplies) the terrain-only `result`
      # computed above -- terrain horizon blocking is already folded into
      # `w_ij` per direction, so `acc / n_az` is already the combined
      # terrain+canopy sky-view factor, not a separate canopy-only fraction.
      result <- acc / n_az

    } else if (model == "twostream") {

      # ---- Two-stream diffuse fraction (scalar, no hgt) ----
      x_m    <- .to_mat(x,    1.0,  nrows, ncols)
      gref_m <- .to_mat(gref, 0.15, nrows, ncols)
      dfrac  <- canopy_diffuse_frac_cpp(pai_m, lref_m, ltra_m, x_m, gref_m)
      result <- result * dfrac

    } else {

      # ---- Goudriaan simplified diffuse fraction (scalar, no hgt) ----
      # exp(-sqrt(alpha) * K_diffuse * pai), K_diffuse = 1, alpha = leaf
      # absorptivity = 1 - lref - ltra. Independent of `x`/`gref` (ignored
      # for this model -- unlike the two-stream fraction, which uses both).
      alpha_m <- pmax(1 - lref_m - ltra_m, 0)
      dfrac   <- exp(-sqrt(alpha_m) * pai_m)
      result  <- result * dfrac
    }
  }

  template <- if (!is.null(dtm)) dtm else hor[[1]]
  out <- terra::rast(template)
  terra::values(out) <- result
  names(out) <- "sky_view_factor"
  if (!is.null(dtm)) out <- terra::mask(out, dtm)
  out
}


#' Two-stream canopy radiation model
#'
#' Computes what fraction of the sunlight arriving above a canopy actually
#' makes it down to the ground beneath, worked out separately for the sun's
#' direct beam and for diffuse sky light (light scattered across the rest of
#' the sky), and accounting for the terrain's slope and how reflective the
#' leaves and ground are. This is a more physically detailed calculation
#' than simply assuming a fixed proportion of light is blocked: it also
#' tracks light that bounces between leaves and the ground one or more times
#' before it eventually reaches the ground or escapes back out, rather than
#' treating any light that touches a leaf as lost. (Technically, this is the
#' Sellers (1985) two-stream approximation, with the sun's angle relative to
#' the slope handled via the Campbell & Norman (1998) extinction-coefficient
#' formulation.)
#'
#' @details
#' Each output is a fraction between 0 and 1: multiply it by the
#' corresponding above-canopy sunlight amount (`swdir` for the sun's direct
#' beam, `swdif` for diffuse sky light) to get the actual amount reaching
#' the ground. Three fractions are returned per cell, describing three
#' different fates for the light:
#'
#' \describe{
#'   \item{`fdirdown`}{How much of the sun's direct beam reaches the ground
#'     still as a direct beam, having passed through the canopy without
#'     being intercepted. Depends on the sun's angle, the slope, and how
#'     dense the canopy is (`pai`).}
#'   \item{`fdifdown`}{How much of the diffuse sky light reaches the
#'     ground. Doesn't depend on the sun's position at all - only on how
#'     dense the canopy is and how reflective/absorbent the leaves are.}
#'   \item{`fdirscat`}{How much of the sun's direct beam gets caught by
#'     the canopy, bounced around, and eventually reaches the ground as
#'     scattered, diffuse-like light rather than passing straight through.
#'     This is the extra light that a simpler "canopy just dims light
#'     passing through it" model would miss entirely.}
#' }
#'
#' To get the total below-canopy sunlight, add all three together: total
#' downward = `swdir * fdirdown + swdif * fdifdown + swdir * fdirscat`,
#' where `swdir` and `swdif` are the above-canopy direct and diffuse
#' sunlight amounts you're starting from.
#'
#' **Leaf reflectance and transmittance** (`lref`, `ltra`) default to
#' typical values for broadband sunlight (0.4 and 0.2). If you're working
#' with a specific part of the light spectrum instead (e.g.
#' photosynthetically active radiation, PAR), supply your own values as
#' rasters.
#'
#' **Plant area index convention** (`pai_horiz`): `pai` can be supplied per
#' unit horizontal ground area (`TRUE`, default - the usual convention) or
#' per unit sloped ground-surface area (`FALSE`), the two differing by
#' `cos(slope)` on sloping ground. The direct-beam fraction (`fdirdown`)
#' treats the canopy as a slope-parallel layer and so needs `pai` per sloped
#' area, whereas the diffuse fraction (`fdifdown`) solves the two-stream
#' equations for a horizontally-homogeneous layer and so needs `pai` per
#' horizontal area; the scattered fraction (`fdirscat`) draws on both.
#' Whichever convention you supply, each fraction is computed from `pai`
#' converted to its own convention, so a fixed physical canopy attenuates
#' the beam identically either way and only the reported area basis differs.
#' On flat ground `pai_horiz` has no effect. Leaf angles are always
#' referenced to the horizontal and never rotate with the ground.
#'
#' When you supply more than one time in `tme`, results are averaged across
#' all of them by default (`output = "mean"`). Night-time timesteps count
#' as zero toward the direct and scattered fractions, so supply only
#' daytime times if you want a daytime-only average. `fdifdown` doesn't
#' depend on the sun's position at all, so it comes out identical
#' regardless of which timesteps you include.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must use a projected CRS with metre units. Slope and aspect are
#'   worked out from it automatically.
#' @param pai A `SpatRaster` of plant area index (m^2 m^-2) - a measure of
#'   how much leaf, stem and branch material sits above each point (higher
#'   values mean a denser canopy that blocks more light).
#' @param tme A `POSIXct` or `POSIXlt` vector of times **in UTC**. Unlike
#'   some other functions in this package, the sun's position here is
#'   always worked out separately for every single grid cell from its own
#'   geographic coordinates (more accurate over large maps, but slower).
#' @param x Leaf angle distribution parameter. Controls how the canopy's
#'   leaves are typically angled, which affects how much light passing
#'   through gets blocked versus slipping through gaps: leave it at the
#'   default (1) for a "typical" mixed canopy if unsure, use a value close
#'   to 0 for canopies of mostly upright/vertical leaves (e.g. many
#'   grasses and some conifers), or a large value (or `Inf`) for canopies
#'   of mostly flat, horizontal leaves (e.g. many broadleaf trees).
#'   (Technically, the ratio of leaf area projected onto a horizontal plane
#'   versus a vertical one, following the Campbell & Norman (1998)
#'   ellipsoidal convention.) May be a scalar, `SpatRaster`, or `NA`/`NULL`.
#' @param lref Leaf reflectance (0--1) - the fraction of light hitting a
#'   leaf that bounces straight back off it. `NULL` or `NA` uses 0.4.
#'   May be a scalar, `SpatRaster`, or `NA`/`NULL`.
#' @param ltra Leaf transmittance (0--1) - the fraction of light hitting a
#'   leaf that passes straight through it. `NULL` or `NA` uses 0.2.
#'   Note `lref + ltra` must be < 1 (a leaf can't reflect and transmit more
#'   light than hits it). May be a scalar, `SpatRaster`, or `NA`/`NULL`.
#' @param gref Ground reflectance (0--1) - the fraction of light hitting
#'   the bare ground beneath the canopy that bounces back upward (higher
#'   for pale, dry, or snow-covered ground; lower for dark, wet soil).
#'   Default 0.15. May be a scalar or `SpatRaster`.
#' @param output Character. `"mean"` (default) returns a single 3-layer
#'   `SpatRaster` averaged across all the timesteps in `tme`. `"stack"`
#'   instead returns a named list of three `SpatRaster` objects
#'   (`fdirdown`, `fdifdown`, `fdirscat`), each with one layer per
#'   timestep, so you can inspect individual timesteps rather than only
#'   the average.
#' @param pai_horiz Logical. Declares which area the supplied `pai` is
#'   measured against: `TRUE` (default) - per unit **horizontal** ground
#'   area (the usual leaf/plant area index convention); `FALSE` - per unit
#'   **sloped** ground-surface area. Identical on flat ground; on a slope
#'   they differ by `cos(slope)`. Unlike [solarindex()] (which is purely a
#'   direct-beam calculation), the two streams here need `pai` in different
#'   conventions - the direct beam is attenuated as a slope-parallel layer
#'   (per sloped area) while the diffuse two-stream assumes a
#'   horizontally-homogeneous layer (per horizontal area) - so each is fed
#'   `pai` converted to its own convention internally. Set `pai_horiz` to
#'   match how your `pai` data was derived; it changes the result on
#'   sloping ground. See Details.
#' @param fillna Logical. If `TRUE`, `NA` cells in `dtm` (e.g. sea, or any
#'   other no-data gap), and the outer ring of `dtm` itself, are
#'   temporarily filled from their nearest real neighbour before
#'   slope/aspect are worked out, and every output layer is masked back to
#'   `dtm`'s own `NA` pattern afterwards - avoids slope/aspect (and
#'   everything downstream of them) coming out `NA` at land cells that
#'   merely touch an `NA` cell, e.g. an entire coastline's worth of cells
#'   picking up `NA` purely for being next to the sea, and at the raster's
#'   own edge, where the same problem occurs because there's no
#'   neighbouring cell at all just past the boundary. Default `FALSE`.
#'   See [solarindex()]'s own `fillna` for the full rationale.
#'
#' @return A 3-layer `SpatRaster` (layers: `fdirdown`, `fdifdown`,
#'   `fdirscat`) when `output = "mean"`, or a named list of three
#'   `SpatRaster` objects when `output = "stack"`.
#' @export
#'
#' @seealso [solarindex()], [horizon()], [canopy3d()]
#'
#' @references
#' Sellers, P.J. (1985). Canopy reflectance, photosynthesis and
#' transpiration. \emph{International Journal of Remote Sensing}, 6,
#' 1335--1372.
#'
#' Campbell, G.S. & Norman, J.M. (1998). \emph{An Introduction to
#' Environmental Biophysics}, 2nd ed. Springer, New York.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm  <- rast(dtm100m)
#' pai  <- rast(pai)          # plant area index (Lizard, Cornwall)
#' vegx <- rast(vegx)         # leaf angle distribution parameter
#'
#' tme <- as.POSIXct("2023-06-21 10:00", tz = "UTC")
#'
#' # Single instant, default leaf angle distribution (spherical)
#' ts <- twostream(dtm, pai, tme)
#' plot(ts)
#'
#' # Single instant, with the real (spatially variable) leaf angle distribution
#' ts_x <- twostream(dtm, pai, tme, x = vegx)
#'
#' # Mean over a day
#' tme_day <- seq(as.POSIXct("2023-06-21 05:00", tz = "UTC"),
#'                as.POSIXct("2023-06-21 20:00", tz = "UTC"),
#'                by = "1 hour")
#' ts_day <- twostream(dtm, pai, tme_day, x = vegx)
#'
#' # Hourly stack with spatially variable leaf properties
#' lref <- rast(dtm); terra::values(lref) <- 0.35
#' ts_stack <- twostream(dtm, pai, tme_day, x = vegx, lref = lref, output = "stack")
#'
#' # DTM with NA cells at a coastline: avoid an NA halo on land
#' ts_coast <- twostream(dtm, pai, tme, fillna = TRUE)
#' }
twostream <- function(dtm, pai, tme,
                      x      = NULL,
                      lref   = NULL,
                      ltra   = NULL,
                      gref   = 0.15,
                      output = "mean",
                      pai_horiz = TRUE,
                      fillna = FALSE) {

  stopifnot(inherits(dtm, "SpatRaster"))
  .check_projected(dtm)
  stopifnot(inherits(pai, "SpatRaster"))
  stopifnot(is.logical(pai_horiz), length(pai_horiz) == 1L, !is.na(pai_horiz))
  output <- match.arg(output, c("mean", "stack"))

  # ---- Optional NA fill (see `fillna` below / solarindex()'s own version
  # of this for the full rationale) -- NA cells in `dtm`, and the outer
  # ring of `dtm` itself (whose focal window reaches past the raster's own
  # edge), are temporarily filled from their nearest real neighbour via
  # .fillna_dtm() so terra::terrain()'s slope/aspect calculation doesn't
  # come out NA at every valid land cell touching one (e.g. a coastline) or
  # at the raster's edge, then every output layer is masked back to `dtm`'s
  # real NA pattern below.
  dtm_na <- dtm
  fill_ext <- NULL
  if (isTRUE(fillna)) {
    fill_ext <- .fillna_dtm(dtm)
    dtm <- fill_ext$interior
  }
  .mask_fillna <- function(r) if (isTRUE(fillna)) terra::mask(r, dtm_na) else r

  # Slope and aspect from dtm -- when fillna = TRUE, computed on
  # fill_ext$extended (1 cell wider than `dtm` in every direction) so the
  # outer ring gets a real value too, then cropped back to `dtm`'s extent.
  slp <- if (isTRUE(fillna)) {
    terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm_na)
  } else terra::terrain(dtm, "slope", unit = "degrees")
  asp <- if (isTRUE(fillna)) {
    terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm_na)
  } else terra::terrain(dtm, "aspect", unit = "degrees")
  sl_m <- as.matrix(slp, wide = TRUE)
  as_m <- as.matrix(asp, wide = TRUE)
  pa_m <- as.matrix(pai, wide = TRUE)

  nr <- terra::nrow(dtm)
  nc <- terra::ncol(dtm)

  x_m    <- .to_mat(x,    1.0,  nr, nc)
  lref_m <- .to_mat(lref, 0.4,  nr, nc)
  ltra_m <- .to_mat(ltra, 0.2,  nr, nc)
  gref_m <- .to_mat(gref, 0.15, nr, nc)
  .check_leaf_optics(lref_m, ltra_m)

  # Lat/lon from DTM
  ll    <- .latslonsfromr(dtm)
  lat_m <- ll$lats
  lon_m <- ll$lons

  # Time components
  if (!inherits(tme, "POSIXlt")) tme <- as.POSIXlt(tme, tz = "UTC")
  years  <- as.integer(tme$year + 1900L)
  months <- as.integer(tme$mon  + 1L)
  days   <- as.integer(tme$mday)
  hours  <- tme$hour + tme$min / 60.0 + tme$sec / 3600.0
  n_t    <- length(years)

  # Helper: matrix -> named SpatRaster layer
  .make_layer <- function(mat, nm) {
    r <- terra::rast(dtm)
    terra::values(r) <- mat
    names(r) <- nm
    r
  }

  # Mean mode (or single timestep)
  if (output == "mean" || n_t == 1L) {
    res <- twostream_tme_cpp(sl_m, as_m, pa_m, lref_m, ltra_m, x_m, gref_m,
                              lat_m, lon_m, years, months, days, hours,
                              pai_horiz)
    return(.mask_fillna(terra::rast(list(
      .make_layer(res$fdirdown, "fdirdown"),
      .make_layer(res$fdifdown, "fdifdown"),
      .make_layer(res$fdirscat, "fdirscat")
    ))))
  }

  # Stack mode: one layer per timestep for each output
  dir_layers  <- vector("list", n_t)
  dif_layers  <- vector("list", n_t)
  scat_layers <- vector("list", n_t)

  for (i in seq_len(n_t)) {
    res <- twostream_tme_cpp(sl_m, as_m, pa_m, lref_m, ltra_m, x_m, gref_m,
                              lat_m, lon_m,
                              years[i], months[i], days[i], hours[i],
                              pai_horiz)
    tstr <- format(tme[i], "%Y%m%d_%H%M")
    dir_layers [[i]] <- .make_layer(res$fdirdown, paste0("fdirdown_", tstr))
    dif_layers [[i]] <- .make_layer(res$fdifdown, paste0("fdifdown_", tstr))
    scat_layers[[i]] <- .make_layer(res$fdirscat, paste0("fdirscat_", tstr))
  }

  list(
    fdirdown = .mask_fillna(terra::rast(dir_layers)),
    fdifdown = .mask_fillna(terra::rast(dif_layers)),
    fdirscat = .mask_fillna(terra::rast(scat_layers))
  )
}


#' 3D canopy radiation model
#'
#' Computes what fraction of the sunlight arriving above a canopy actually
#' makes it down to the ground beneath - the same question [twostream()]
#' answers - but by tracing the sun's and sky's actual paths through the
#' real, uneven shape of the canopy and terrain nearby, rather than
#' assuming the canopy stretches out as a perfectly uniform blanket in
#' every direction the way [twostream()] does. That makes this function
#' better suited to patchy or uneven vegetation - gaps, clearings, a taller
#' stand of trees next door - at the cost of being considerably slower to
#' compute. It's the 3D counterpart to [twostream()], built from the same
#' underlying technique already used inside [solarindex()] (for direct
#' sunlight) and [skyview()] (for diffuse sky light) when their `hgt`
#' option is switched on.
#'
#' @details
#' Two fractions are returned per cell, each between 0 and 1:
#'
#' \describe{
#'   \item{`fdirdown`}{How much of the sun's direct beam reaches the
#'     ground, worked out by tracing the sun's actual path sideways across
#'     neighbouring points' canopies - rather than assuming it travels
#'     straight down through only the target point's own canopy. Depends
#'     on the sun's position and the shape of the canopy and terrain
#'     nearby.}
#'   \item{`fdifdown`}{How much of the diffuse sky light reaches the
#'     ground, worked out the same way but checked across many directions
#'     of the sky (`ndir` compass directions x `n_elev` bands of sky
#'     height) rather than just toward the sun. Doesn't depend on the
#'     sun's position - only on the canopy's density, its leaf properties,
#'     and its local shape.}
#' }
#'
#' Unlike [twostream()], there's no third `fdirscat` output here: this
#' simpler calculation doesn't track light bouncing between leaves and the
#' ground, which is what makes it practical to trace individual paths
#' through the canopy at all (see [skyview()]'s documentation for more on
#' this tradeoff between accuracy and the ability to do a 3D calculation).
#'
#' **On slope**: you won't find any separate handling of slope here, unlike
#' [twostream()] - because the calculation traces the sun's actual path
#' through the real shape of the terrain and canopy, a slope's effect on
#' how much canopy the light has to pass through comes out naturally, with
#' no separate correction needed. What this function does *not* do is
#' account for how squarely the light then lands on a sloped receiving
#' surface once it arrives there - for that, combine the result with a
#' (vegetation-free) result from [solarindex()] separately, the same way
#' [solarindex()]'s own `hgt` option combines the two internally.
#'
#' **On terrain shading**: this function only accounts for light being
#' blocked by vegetation, not by bare terrain (e.g. a distant ridge with no
#' canopy on it) - tracing a path through open, canopy-free ground doesn't
#' block anything here, unlike [horizon()]'s dedicated terrain-shading
#' calculation. If you also want terrain shading, combine the result with
#' [horizon()] or [skyview()] separately.
#'
#' **On the sun's position**: because tracing the sun's path is done once
#' for the whole map rather than separately at every point, this function
#' - like [solarindex()]'s equivalent `hgt` option - assumes one shared sun
#' position across the whole map at each timestep (worked out at the
#' middle of the map), rather than [twostream()]'s more precise per-point
#' calculation. A warning appears if your map is large enough (roughly more
#' than 1 degree of latitude/longitude) that this assumption may start to
#' matter.
#'
#' When you supply more than one time in `tme`, `fdirdown` is averaged
#' across all of them by default (`output = "mean"`), counting night-time
#' timesteps as zero - matching the convention used by [twostream()] and
#' [solarindex()]. `fdifdown` never changes with the sun's position, so it
#' comes out identical regardless of which timesteps you include.
#'
#' This function is considerably slower than [twostream()], since it
#' re-examines the whole map once for every sun position (for `fdirdown`)
#' and once for every sky direction checked (for `fdifdown` - `ndir *
#' n_elev` checks in total, done once regardless of how many timesteps you
#' supply). A warning appears if the number of sun positions exceeds 100, or
#' if the diffuse fraction's estimated runtime (based on both the direction
#' count and the raster's size, since a big map is slower to re-examine
#' than a small one) looks like it's heading past around 15 seconds; the
#' defaults (`ndir = 32`, `n_elev = 6`) run quickly on an ordinary-sized
#' DTM.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model. Must use a projected CRS with metre units.
#' @param pai A `SpatRaster` of plant area index (m^2 m^-2) - a measure of
#'   how much leaf, stem and branch material sits above each point (higher
#'   values mean a denser canopy that blocks more light).
#' @param hgt A `SpatRaster` of canopy height above the local ground
#'   (metres), matching `dtm`/`pai`'s grid. Cells with `pai > 0` but
#'   `hgt <= 0` (or missing) are treated as having no canopy there, with a
#'   warning.
#' @param tme A set of dates/times (`POSIXct`/`POSIXlt`, **in UTC**). The
#'   sun's position for tracing direct sunlight is worked out once per
#'   timestep, using the middle of the map (see Details).
#' @param x Leaf angle distribution parameter, used for the direct-beam
#'   fraction. Controls how the canopy's leaves are typically angled,
#'   which affects how much light passing through gets blocked versus
#'   slipping through gaps: leave it at the default (1) for a "typical"
#'   mixed canopy if unsure, use a value close to 0 for canopies of mostly
#'   upright/vertical leaves (e.g. many grasses and some conifers), or a
#'   large value (or `Inf`) for canopies of mostly flat, horizontal leaves
#'   (e.g. many broadleaf trees). Follows the same Campbell & Norman (1998)
#'   ellipsoidal convention as [solarindex()] and [twostream()]. May be a
#'   scalar or `SpatRaster`.
#' @param lref Leaf reflectance (0--1), used for the diffuse fraction - the
#'   fraction of light hitting a leaf that bounces straight back off it.
#'   `NULL` or `NA` uses 0.4 (matching [twostream()]'s default). May be a
#'   scalar or `SpatRaster`.
#' @param ltra Leaf transmittance (0--1), used for the diffuse fraction -
#'   the fraction of light hitting a leaf that passes straight through it.
#'   `NULL` or `NA` uses 0.2 (matching [twostream()]'s default). Note
#'   `lref + ltra` must be < 1 (a leaf can't reflect and transmit more
#'   light than hits it). May be a scalar or `SpatRaster`.
#' @param ndir Integer. Number of compass directions checked around the sky
#'   when tracing diffuse sky light. Default 32 (matching [skyview()]'s
#'   default) - higher gives a more precise result at greater computational
#'   cost.
#' @param n_elev Integer. Number of bands of sky height checked (from the
#'   horizon up to directly overhead), also for diffuse sky light. Default
#'   6 (matching [skyview()]'s default).
#' @param nstep Integer. Controls how thoroughly neighbouring points are
#'   searched when tracing a path, passed through to the underlying
#'   calculation. `0` (default) checks carefully all the way to the edge of
#'   the map; a positive integer trades some accuracy for speed. See
#'   [horizon()] for details.
#' @param output Character. `"mean"` (default) returns a single 2-layer
#'   `SpatRaster` (`fdirdown` averaged across timesteps, `fdifdown`
#'   unchanged since it doesn't vary by timestep). `"stack"` instead
#'   returns a named list of two `SpatRaster` objects: `fdirdown` (one
#'   layer per timestep) and `fdifdown` (a single layer).
#'
#' @return A 2-layer `SpatRaster` (layers: `fdirdown`, `fdifdown`) when
#'   `output = "mean"`, or a named list of two `SpatRaster` objects when
#'   `output = "stack"`.
#' @export
#'
#' @seealso [twostream()], [solarindex()], [skyview()], [horizon()]
#'
#' @references
#' Goudriaan, J. (1977). \emph{Crop micrometeorology: a simulation study}.
#' Centre for Agricultural Publishing and Documentation, Wageningen.
#'
#' Campbell, G.S. & Norman, J.M. (1998). \emph{An Introduction to
#' Environmental Biophysics}, 2nd ed. Springer, New York.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm  <- rast(dtm100m)
#' pai  <- rast(pai)          # plant area index (Lizard, Cornwall)
#' hgt  <- rast(hgt)          # canopy height (Lizard, Cornwall)
#' vegx <- rast(vegx)         # leaf angle distribution parameter
#'
#' tme <- as.POSIXct("2023-06-21 10:00", tz = "UTC")
#'
#' # Single instant, default leaf angle distribution / optical properties
#' c3 <- canopy3d(dtm, pai, hgt, tme)
#' plot(c3)
#'
#' # With the real leaf angle distribution and a coarser (faster) diffuse
#' # quadrature
#' c3_x <- canopy3d(dtm, pai, hgt, tme, x = vegx, ndir = 16L, n_elev = 4L)
#'
#' # Mean over a day, hourly stack
#' tme_day <- seq(as.POSIXct("2023-06-21 05:00", tz = "UTC"),
#'                as.POSIXct("2023-06-21 20:00", tz = "UTC"),
#'                by = "1 hour")
#' c3_stack <- canopy3d(dtm, pai, hgt, tme_day, x = vegx, output = "stack")
#' }
canopy3d <- function(dtm, pai, hgt, tme,
                     x      = NULL,
                     lref   = NULL,
                     ltra   = NULL,
                     ndir   = 32L,
                     n_elev = 6L,
                     nstep  = 0L,
                     output = "mean") {

  stopifnot(inherits(dtm, "SpatRaster"))
  .check_projected(dtm)
  stopifnot(inherits(pai, "SpatRaster"))
  stopifnot(inherits(hgt, "SpatRaster"))
  output <- match.arg(output, c("mean", "stack"))

  nr <- terra::nrow(dtm)
  nc <- terra::ncol(dtm)
  res <- terra::res(dtm)[1]

  elev_m <- as.matrix(dtm, wide = TRUE)
  pai_m  <- as.matrix(pai, wide = TRUE)
  hgt_m  <- as.matrix(hgt, wide = TRUE)

  x_m    <- .to_mat(x,    1.0, nr, nc)
  lref_m <- .to_mat(lref, 0.4, nr, nc)
  ltra_m <- .to_mat(ltra, 0.2, nr, nc)
  .check_leaf_optics(lref_m, ltra_m)

  bad_hgt <- !is.na(pai_m) & pai_m > 0 & (is.na(hgt_m) | hgt_m <= 0)
  if (any(bad_hgt))
    warning(sum(bad_hgt), " cell(s) have pai > 0 but hgt <= 0 (or NA); ",
            "vegetation attenuation is treated as zero (no canopy) at ",
            "those cells.", call. = FALSE)

  # ---- Diffuse fraction: solar-angle independent, computed once ----
  # Ray-marched azimuth x elevation-bin quadrature using Goudriaan's
  # diffuse approximation (see veg_path_diffuse_cpp()). Deliberately no
  # terrain-horizon truncation here (unlike skyview()'s hgt path): this is
  # meant as a pure canopy-transmission quantity, matching twostream()'s
  # fdifdown, which is also terrain/slope independent -- not a combined
  # terrain+canopy sky-view factor. Quadrature weights use the exact
  # closed-form sin(e)cos(e) antiderivative (sin^2(e)) per elevation bin,
  # so that if pai were 0 everywhere, summing all bins would telescope to
  # exactly 1 (the full unobstructed hemisphere), regardless of n_elev.
  azimuths <- seq(0, 360 - 360 / ndir, length.out = ndir)
  edges    <- seq(0, 90, length.out = n_elev + 1L)
  mids     <- (edges[-1L] + edges[-(n_elev + 1L)]) / 2

  # Cost estimate accounts for raster size as well as direction count -- see
  # .warn_ray_march_cost() in utils.R for the calibration/rationale.
  n_dir_total <- ndir * n_elev
  .warn_ray_march_cost(n_dir_total, nr, nc, "canopy3d()'s diffuse fraction")

  acc_dif <- matrix(0.0, nr, nc)
  for (j in seq_len(ndir)) {
    for (i in seq_len(n_elev)) {
      tau_ij  <- veg_path_diffuse_cpp(elev_m, hgt_m, pai_m, lref_m, ltra_m,
                                       res, 90 - mids[i], azimuths[j], nstep)
      w_i     <- sin(edges[i + 1L] * pi / 180)^2 - sin(edges[i] * pi / 180)^2
      acc_dif <- acc_dif + tau_ij * w_i
    }
  }
  fdif <- acc_dif / ndir

  # ---- Direct fraction: per solar position, single shared position ----
  # (ray-marching needs one shared zenith/azimuth per call -- same
  # constraint as solarindex()'s hgt path; see Details above)
  ll   <- .latslonsfromr(dtm)
  cr   <- (nr + 1L) %/% 2L
  cc   <- (nc + 1L) %/% 2L
  lat0 <- ll$lats[cr, cc]
  lon0 <- ll$lons[cr, cc]

  lat_span <- diff(range(ll$lats, na.rm = TRUE))
  lon_span <- diff(range(ll$lons, na.rm = TRUE))
  if (lat_span > 1 || lon_span > 1)
    warning("DTM spans more than 1 degree; solar position for the direct-",
            "beam ray-marching is computed only at the raster centre, ",
            "which may introduce errors near the grid edges.",
            call. = FALSE)

  if (!inherits(tme, "POSIXlt")) tme <- as.POSIXlt(tme, tz = "UTC")
  years  <- as.integer(tme$year + 1900L)
  months <- as.integer(tme$mon  + 1L)
  days   <- as.integer(tme$mday)
  hours  <- tme$hour + tme$min / 60.0 + tme$sec / 3600.0
  n_t    <- length(years)

  if (n_t > 100L)
    warning(n_t, " solar position(s) requested for direct-beam ",
            "ray-marching; this is considerably more expensive per ",
            "position than a closed-form model. Consider reducing the ",
            "number of timesteps.", call. = FALSE, immediate. = TRUE)

  pos_mat <- solar_pos_cpp(lat0, lon0, years, months, days, hours)
  zen_v   <- pos_mat[, 1L]
  azi_v   <- pos_mat[, 2L]

  .make_layer <- function(mat, nm) {
    r <- terra::rast(dtm)
    terra::values(r) <- mat
    names(r) <- nm
    r
  }

  if (output == "mean" || n_t == 1L) {
    acc_dir <- matrix(0.0, nr, nc)
    for (i in seq_len(n_t)) {
      if (zen_v[i] >= 90.0) next
      acc_dir <- acc_dir + veg_path_atten_cpp(elev_m, hgt_m, pai_m, x_m,
                                               res, zen_v[i], azi_v[i], nstep)
    }
    fdir <- acc_dir / n_t
    return(terra::rast(list(
      .make_layer(fdir, "fdirdown"),
      .make_layer(fdif, "fdifdown")
    )))
  }

  # Stack mode: one fdirdown layer per timestep; fdifdown is a single
  # layer (solar-angle independent, so identical across all timesteps).
  dir_layers <- vector("list", n_t)
  for (i in seq_len(n_t)) {
    mat <- if (zen_v[i] >= 90.0) {
      matrix(0.0, nr, nc)
    } else {
      veg_path_atten_cpp(elev_m, hgt_m, pai_m, x_m, res,
                          zen_v[i], azi_v[i], nstep)
    }
    tstr <- format(tme[i], "%Y%m%d_%H%M")
    dir_layers[[i]] <- .make_layer(mat, paste0("fdirdown_", tstr))
  }

  list(
    fdirdown = terra::rast(dir_layers),
    fdifdown = .make_layer(fdif, "fdifdown")
  )
}
