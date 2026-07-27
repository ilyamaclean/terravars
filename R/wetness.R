# ============================================================
# wetness.R
#
# Soil moisture / hydrological functions:
#   - twi()   topographic wetness index (standard / modified / saga)
# ============================================================
# ============================================================================ #
# ~~~~~~~~~ Basin delineation worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Delineate hydrological basins
#' @description Assigns every cell of a DTM to a basin. Works on the
#'   principle that any cell adjoining and uphill of another cell belongs to
#'   the same basin as that cell, so each basin grows outward and upward
#'   from a local pit until every reachable cell has been claimed. Cells at
#'   exactly the same elevation as a claimed neighbour always merge into the
#'   same basin, regardless of `method`. `NA` cells (e.g. sea) are left `NA`
#'   -- never assigned a basin.
#' @param dtm A `SpatRaster` of elevation values (metres).
#' @param method Character. `"steepest"` (default): a strictly higher
#'   neighbour only merges into the basin that owns its own single steepest
#'   downhill neighbour (ties for steepest are broken arbitrarily but
#'   deterministically), so each sloped cell's fate is decided purely by its
#'   local gradient rather than by the order basins happen to be discovered
#'   in. `"any"`: a strictly higher neighbour merges into a basin if it
#'   adjoins *any* already-claimed cell of that basin (this package's
#'   original behaviour, before `method` was introduced). Flats/ties behave
#'   identically under both methods.
#' @return A `SpatRaster` of basin IDs (integer), on the same grid as `dtm`.
#'   Basin ID *numbers* are assigned in whatever order basins happen to be
#'   discovered internally and carry no meaning beyond grouping.
#' @export
#' @importFrom Rcpp sourceCpp
#' @useDynLib terravars, .registration = TRUE
#'
#' @seealso [basinmerge()], [fillsinks()]
basindelin <- function(dtm, method = "steepest") {
  stopifnot(inherits(dtm, "SpatRaster"))
  method <- match.arg(method, c("any", "steepest"))
  dm  <- .is(dtm)
  dm2 <- .pad1(dm, 9999, na_as_fill = TRUE)
  # (2) create blank basin file
  bsn <- dm2 * NA
  bsn <- basinCpp(dm2, bsn, method)
  dd  <- dim(bsn)
  bsn <- bsn[2:(dd[1] - 1), 2:(dd[2] - 1)]
  if (class(bsn)[1] != 'matrix') bsn <- matrix(bsn, ncol = ncol(dtm), nrow = nrow(dtm))
  r <- .rast(bsn, dtm)
  return(r)
}

# ============================================================================ #
# ~~~~~~~~~ Basin merging worker functions here  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' @title Renumber a vector/matrix/array of IDs sequentially, gap-free
#' @description Remaps whatever distinct non-NA values are present in `v`
#'   to consecutive integers 1, 2, 3, ... in ascending order of their
#'   original value, preserving `v`'s shape and its `NA` pattern. Used as
#'   an explicit final step wherever an ID field could otherwise end up
#'   with gaps (e.g. after merging groups together, where surviving group
#'   labels are a sparse subset of the original ones).
#' @noRd
.renumber_seq <- function(v) {
  present <- is.na(v) == FALSE
  u   <- sort(unique(v[present]))
  out <- v
  out[present] <- match(v[present], u)
  out
}

#' Merge basins separated by a shallow barrier
#'
#' Cold air (and other dense-air drainage) can spill over a shallow ridge
#' separating two adjacent basins much more readily than it can flow uphill
#' over a tall one. Starting from a set of basins (e.g. from
#' [basindelin()]), this merges together any basins whose lowest crossing
#' point (their pour point, or col) sits close to basin-floor elevation,
#' treating the two as one basin in their entirety -- every cell of both,
#' not just the cells immediately straddling the boundary -- since once air
#' can spill across, the two basins behave as a single connected pooling
#' unit. Merges chain transitively: if basin A merges with B, and B
#' separately merges with C somewhere else, A/B/C all end up as one basin
#' even if A and C never directly touch.
#'
#' @details
#' For every pair of basins that share a boundary, the pour point (the
#' lowest point at which you could cross from one basin into the other) is
#' found by scanning every adjacent cell pair across the shared edge (all 8
#' neighbour directions), taking the higher of the two cells' elevations at
#' each crossing (the elevation you'd have to reach to cross there), and
#' then taking the minimum of that over the whole boundary. Two basins are
#' merged when this pour point sits within `boundary` metres of the lower
#' of the two basins' own floor elevations (each basin's minimum `dtm`
#' value) -- i.e. the basin only has to fill by a small amount before it
#' spills into its neighbour. A boundary that runs along a high ridge but
#' happens to be locally smooth somewhere along its length is correctly
#' *not* treated as shallow, since what matters is the height of the
#' lowest crossing above the basin floors, not how gently the terrain
#' happens to slope at any one point along the edge. Every basin is then
#' merged with every other basin it's connected to this way, directly or
#' transitively (e.g. via a third basin acting as a stepping stone), and
#' the merged groups are renumbered sequentially. `NA` cells (e.g. sea) are
#' left `NA` and never considered part of a basin or a barrier.
#'
#' @param dtm A `SpatRaster` of elevation values (metres).
#' @param bsn A `SpatRaster` of basin IDs (integer), on the same grid as
#'   `dtm` -- typically the output of [basindelin()].
#' @param boundary Numeric. How close (metres) the pour point between two
#'   basins must be to the lower of their two floor elevations for the
#'   basins to be treated as connected and merged. Default `0.25`.
#'
#' @return A `SpatRaster` of basin IDs (integer), on the same grid as
#'   `dtm`/`bsn`, with any basins connected by a shallow-enough barrier
#'   merged into a single ID. Basin ID *numbers* are assigned in whatever
#'   order merged groups happen to be discovered internally and carry no
#'   meaning beyond grouping.
#' @export
#'
#' @seealso [basindelin()], [fillsinks()]
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' bsn <- basindelin(dtm)
#' bsn_merged <- basinmerge(dtm, bsn, boundary = 0.25)
#' plot(bsn_merged)
#' }
basinmerge <- function(dtm, bsn, boundary = 0.25) {
  stopifnot(inherits(dtm, "SpatRaster"))
  stopifnot(inherits(bsn, "SpatRaster"))
  .check_compatible(dtm, bsn, "dtm", "bsn")
  bm <- .is(bsn)
  dm <- .is(dtm)
  nr <- dim(bm)[1]
  nc <- dim(bm)[2]

  ids   <- sort(unique(bm[is.na(bm) == FALSE]))
  n_ids <- length(ids)
  if (n_ids == 0) return(.rast(bm, dtm))   # nothing to merge

  # Put a 1-cell buffer around both elevation and basin IDs, so every real
  # cell has 8 well-defined neighbours without special-casing the raster
  # edge -- same convention basindelin() uses: NA elevation becomes a 9999
  # sentinel (basinmerge_cpp() treats 9999 as "no real elevation here"),
  # while basin id stays a genuine NA.
  dm2 <- .pad1(dm, 9999, na_as_fill = TRUE)
  bm2 <- .pad1(bm, NA_real_)

  # ---- pour-point-based merge (basin floor elevations, the boundary-
  # crossing scan, and the union-find merge) done in basinmerge_cpp() for
  # speed -- see src/wetness.cpp for the full algorithm explanation. The
  # result is basin ids relabelled by merged group: a sparse/gappy subset
  # of ids at this point, not yet sequential. ----
  bmerged <- basinmerge_cpp(dm2, bm2, boundary)
  dd <- dim(bmerged)
  bmerged <- bmerged[2:(dd[1] - 1), 2:(dd[2] - 1)]
  if (class(bmerged)[1] != 'matrix') bmerged <- matrix(bmerged, ncol = nc, nrow = nr)

  # ---- final step: renumber sequentially so basin IDs come out gap-free ----
  bmerged <- .renumber_seq(bmerged)

  r <- .rast(bmerged, dtm)
  return(r)
}

# ============================================================================ #
# ~~~~~~~~~ Sink/depression filling worker functions here ~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================ #
#' Fill sinks (depressions) in a digital terrain model
#'
#' Real-world DTMs are full of small spurious pits -- noise, a bridge deck
#' cutting across a stream, a data artefact -- that would otherwise trap
#' flow at that cell and stop it going any further downhill. This raises
#' every such pit just high enough that a continuous, never-uphill path
#' exists from it all the way down to an outlet: either the edge of the
#' raster, an `NA`/no-data cell (this package's convention for sea or
#' other missing ground, also used by [basindelin()]/[basinmerge()]), or
#' now optionally a real water body supplied via `waterbody` (see below).
#' Intended as a DTM-conditioning step to run *before* flow-routing
#' functions such as flow accumulation, which otherwise have nowhere for
#' water in a pit to go.
#'
#' @details
#' Works outward from every outlet at once: every cell already touching an
#' outlet is left untouched (it's already valid to drain there), and every
#' other cell is raised to whatever the *lowest* possible elevation is
#' that still keeps a continuous downhill-or-flat path open to some
#' outlet -- never any higher than that true pour point, and never at all
#' if a low enough path already exists without raising it. A pit tucked
#' inside a pit inside a pit fills in stages, each one only as much as it
#' actually needs to escape into the depression immediately outside it.
#'
#' `method = "fill"` (the only option currently implemented) floods each
#' depression flat up to its pour-point elevation, the same way ArcGIS's
#' `Fill` tool or `whitebox`'s `FillDepressions` do. A `"breach"` option
#' (carving a descending channel through the barrier instead of flooding
#' the whole depression flat, preserving more of the original terrain) is
#' planned but not yet implemented.
#'
#' Without `waterbody`, this treats every depression as an artefact to be
#' filled away -- fine for noise, but wrong for a real pond or lake, which
#' would otherwise get flattened up to whatever pour point the algorithm
#' finds rather than kept at its own true water level. When `waterbody` is
#' supplied, every connected water body it identifies (8-connected: cells
#' touching only at a corner still count as one water body) is instead
#' flattened to the *minimum* `dtm` elevation found anywhere within it --
#' surface-return DTMs already record the water surface reasonably well,
#' and taking the min rather than the mean avoids letting noisy near-shore
#' or through-water values pull the plane up -- and every cell of it is
#' fed into the fill as an already-resolved seed at that elevation, the
#' same way an `NA`/sea cell or the raster edge already seeds the fill
#' elsewhere. Ordinary pits everywhere else in `dtm` still fill exactly as
#' before, toward whichever outlet (edge, `NA`, or now also a water body)
#' turns out to be lowest/closest -- water bodies aren't a separate step,
#' they're just one more kind of outlet the same algorithm can drain
#' towards. This also means each water body's true pour point is found for
#' free, with no separate search needed: whichever neighbouring land cell
#' has the lowest-cost path in is naturally the first one the fill reaches.
#' A water body that turns out to have no cells at all (e.g. `water_value`
#' matches nothing present in `waterbody`) is silently ignored, identical
#' to not supplying `waterbody`.
#'
#' @param dtm A `SpatRaster` of elevation values (metres).
#' @param method Character. Currently only `"fill"` is implemented.
#' @param waterbody Optional `SpatRaster`, on the same grid as `dtm`,
#'   marking real water bodies (ponds, lakes, reservoirs) that should be
#'   flattened to their own true water level rather than treated as
#'   artefacts to fill away -- see Details. `NULL` (default): no water
#'   bodies, every depression is filled as an artefact.
#' @param water_value Which value(s) in `waterbody` mark a cell as water.
#'   Default `NA`: cells where `waterbody` is `NA` count as water,
#'   matching this package's `landsea`/`basindelin()`/`basinmerge()`
#'   convention elsewhere. Supply one or more explicit values instead
#'   (e.g. a land-cover code for open water) if `waterbody` uses a
#'   different convention; `NA` and explicit values can be combined in the
#'   same vector.
#'
#' @return A `SpatRaster` of elevation values (metres), on the same grid
#'   as `dtm`, with every interior pit raised just enough to drain, and
#'   (if `waterbody` was supplied) every real water body flattened to its
#'   own minimum elevation. `NA` cells are left `NA`.
#' @export
#'
#' @seealso [basindelin()], [basinmerge()], [flowacc()]
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' dtm_filled <- fillsinks(dtm)
#'
#' # A known pond/lake mask (NA = not water, the default convention) kept
#' # at its own water level instead of being filled away as an artefact
#' dtm_filled2 <- fillsinks(dtm, waterbody = pond_mask)
#'
#' # Same, but the mask uses an explicit land-cover code for water
#' dtm_filled3 <- fillsinks(dtm, waterbody = landcover, water_value = 80)
#' }
fillsinks <- function(dtm, method = "fill", waterbody = NULL, water_value = NA) {
  method <- match.arg(method, c("fill"))
  dm <- .is(dtm)
  nr <- dim(dm)[1]
  nc <- dim(dm)[2]

  # Same padded-matrix convention as basindelin()/basinmerge(): a 1-cell
  # border of the 9999 NA/no-data sentinel, doubling as "every cell on the
  # raster's own edge already touches an outlet" for free.
  dm2 <- .pad1(dm, 9999, na_as_fill = TRUE)

  # ---- optional water body handling ----
  # wseed2 mirrors dm2's padded shape: a real (non-NA) value at a cell
  # means "this is a resolved water-body cell, already at this flattened
  # elevation" and is passed straight through to fillsinks_cpp(), which
  # seeds it into the fill exactly like an outlet (see that function's own
  # comment for how the seeding/growth actually uses it). Defaults to
  # all-NA ("no water bodies") when 'waterbody' isn't supplied.
  wseed2 <- array(NA_real_, dim = c(nr + 2, nc + 2))

  if (!is.null(waterbody)) {
    stopifnot(inherits(waterbody, "SpatRaster"))
    wb <- .is(waterbody)
    if (!identical(dim(wb), dim(dm)))
      stop("'waterbody' must be on the same grid as 'dtm'.", call. = FALSE)

    # Which raster values count as water: NA (default) and/or explicit
    # value(s) supplied via water_value. Built from is.na()/==/"|" only
    # (rather than %in%) so this works identically regardless of terra
    # version. The explicit-value branch is deliberately coerced to FALSE
    # (not NA) wherever waterbody itself is NA, so a cell with no
    # water-body data of its own reads unambiguously as "not water"
    # rather than leaving an indeterminate NA in the mask.
    na_flag  <- is.na(water_value)
    explicit <- water_value[!na_flag]

    wmask <- NULL
    if (any(na_flag)) wmask <- is.na(waterbody)
    for (v in explicit) {
      m <- !is.na(waterbody) & (waterbody == v)
      wmask <- if (is.null(wmask)) m else (wmask | m)
    }

    n_water <- terra::global(wmask, "sum", na.rm = TRUE)[1, 1]
    if (!is.na(n_water) && n_water > 0) {
      wmask <- terra::ifel(wmask, 1, NA)

      # Label each connected water body (8-connectivity, matching every
      # other neighbourhood search in this package), then assign every
      # cell of it the minimum dtm elevation found anywhere in that water
      # body. Done directly on plain matrices with tapply() rather than
      # terra::zonal()/classify() -- an earlier version used that pair
      # instead and silently produced mostly-NA output, most likely from
      # a zonal()/classify() column-order or value-matching mismatch that
      # couldn't be verified without a live R session in this cloud
      # session. tapply() on plain matrices sidesteps that uncertainty
      # entirely: every step here is base R, fully inspectable.
      patches_r <- terra::patches(wmask, directions = 8)
      pm  <- .is(patches_r)
      dmr <- .is(dtm)
      water_idx <- !is.na(pm) & !is.na(dmr)

      if (any(water_idx)) {
        zmin <- tapply(dmr[water_idx], pm[water_idx], min, na.rm = TRUE)

        ws <- matrix(NA_real_, nrow = nr, ncol = nc)
        ws[water_idx] <- zmin[as.character(pm[water_idx])]
        wseed2[2:(nr+1),2:(nc+1)] <- ws
      }
    }
  }

  filled <- fillsinks_cpp(dm2, wseed2)
  dd <- dim(filled)
  filled <- filled[2:(dd[1] - 1), 2:(dd[2] - 1)]
  if (class(filled)[1] != 'matrix') filled <- matrix(filled, ncol = nc, nrow = nr)

  # Restore genuine NA (the 9999 sentinel is purely an internal marker --
  # it must never leak into the returned elevation values).
  filled[is.na(.is(dtm))] <- NA

  r <- .rast(filled, dtm)
  return(r)
}

#' Flow accumulation
#'
#' Computes, for every cell, the total upstream contributing weight (area,
#' by default one unit per cell) that drains through it, following one of
#' three standard flow-routing methods. Real-world DTMs are full of small
#' pits that would otherwise trap flow before it reaches here -- run
#' [fillsinks()] first if you want continuous flow paths rather than flow
#' stopping dead at every unfilled depression (an unfilled pit simply
#' becomes a terminal accumulation point: its own weight plus whatever
#' has drained into it, with nothing propagating further downhill).
#'
#' @param dtm A `SpatRaster` of elevation values (metres).
#' @param method Character. `"d8"` (default): each cell's entire flow
#'   goes to its single steepest downhill neighbour. `"mfd"`: flow splits
#'   across every downhill neighbour, weighted by slope raised to a fixed
#'   convergence exponent (1.1, following Freeman 1991) -- not exposed as
#'   an argument, since a knob nobody can calibrate by feel is worse than
#'   no knob. `"dinf"`: Tarboton's (1997) D-infinity method, which finds
#'   the steepest of the 8 triangular facets around each cell and splits
#'   flow between that facet's two bounding neighbours in proportion to
#'   the facet's flow angle.
#' @param weight Optional `SpatRaster` of per-cell weights (e.g. rainfall
#'   depth, or a fractional area), on the same grid as `dtm`. Defaults to
#'   1 for every cell, so the returned value is a plain cell count (a
#'   proxy for contributing area) unless you provide something else.
#' @param bsn Optional `SpatRaster` of basin IDs (integer), on the same
#'   grid as `dtm` -- typically the output of [basindelin()] or
#'   [basinmerge()]. When supplied, every basin ID is treated as its own
#'   isolated catchment: flow can never cross from one basin into
#'   another, even where the terrain itself would otherwise let it. A
#'   cell right on a basin boundary treats its neighbour in the other
#'   basin exactly like an `NA`/no-data cell -- i.e. as a place flow
#'   cannot go, not as a valid escape route -- so a basin with no
#'   reachable low point of its own (under this restriction) simply
#'   leaves its flow stuck there rather than draining into its neighbour.
#'
#' @return A `SpatRaster` of accumulated flow (same units as `weight`,
#'   or a plain cell count if `weight` is left at its default), on the
#'   same grid as `dtm`. Every cell's value includes its own weight, so a
#'   headwater cell with nothing draining into it still reads its own
#'   weight (1 by default) rather than 0. `NA` cells are left `NA`.
#'
#' @details
#' Flow is routed cell-by-cell from the highest elevation down to the
#' lowest, so that every cell's upstream contributors have already been
#' resolved by the time it forwards its own running total onward.
#'
#' Flat areas (most commonly a depression flooded flat by [fillsinks()])
#' need special handling, since plain slope-based routing has nothing to
#' work with when every neighbour sits at the same elevation. Each flat
#' cell is instead routed toward whichever neighbour of identical
#' elevation is measurably closer -- in flat-only steps -- to an actual
#' downhill exit or the edge of the modelled area, found via a quick
#' preliminary sweep across the flat before routing begins. This lets
#' flow drain all the way across a filled depression to its true pour
#' point without ever needing to alter any elevations.
#'
#' @export
#'
#' @seealso [basindelin()], [basinmerge()], [fillsinks()], [twi()]
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' dtm_filled <- fillsinks(dtm)
#' acc <- flowacc(dtm_filled, method = "mfd")
#' plot(acc)
#'
#' # restrict accumulation to stay within each delineated basin
#' bsn <- basindelin(dtm_filled)
#' acc_by_basin <- flowacc(dtm_filled, bsn = bsn)
#' }
flowacc <- function(dtm, method = c("d8", "mfd", "dinf"), weight = NULL, bsn = NULL) {
  stopifnot(inherits(dtm, "SpatRaster"))
  method <- match.arg(method)
  method_code <- match(method, c("d8", "mfd", "dinf")) - 1L

  dm <- .is(dtm)
  nr <- dim(dm)[1]
  nc <- dim(dm)[2]
  navalues <- is.na(dm)

  if (is.null(weight)) {
    wt <- array(1, dim = c(nr, nc))
  } else {
    stopifnot(inherits(weight, "SpatRaster"))
    .check_compatible(dtm, weight, "dtm", "weight")
    wt <- .is(weight)
  }
  wt[navalues]  <- 0
  wt[is.na(wt)] <- 0

  if (is.null(bsn)) {
    # no basin restriction requested: every real cell is treated as
    # belonging to the same single dummy basin, so the same-basin test
    # inside flowacc_cpp() is always true and flow routes unrestricted.
    bm <- array(NA_integer_, dim = c(nr, nc))
    bm[navalues == FALSE] <- 1L
  } else {
    stopifnot(inherits(bsn, "SpatRaster"))
    .check_compatible(dtm, bsn, "dtm", "bsn")
    bm <- .is(bsn)
    storage.mode(bm) <- "integer"
  }

  # Same padded-matrix convention as basindelin()/basinmerge()/
  # fillsinks(): a 1-cell border doubling as "every cell on the raster's
  # own edge already touches an outlet" for free.
  dm2 <- .pad1(dm, 9999, na_as_fill = TRUE)
  wt2 <- .pad1(wt, 0)
  bm2 <- .pad1(bm, NA_integer_)

  acc <- flowacc_cpp(dm2, wt2, bm2, method_code)
  dd  <- dim(acc)
  acc <- acc[2:(dd[1] - 1), 2:(dd[2] - 1)]
  if (class(acc)[1] != 'matrix') acc <- matrix(acc, ncol = nc, nrow = nr)

  acc[navalues] <- NA

  r <- .rast(acc, dtm)
  return(r)
}

#' Topographic wetness index
#'
#' Estimates how likely each point in the landscape is to be wet, based
#' purely on its position in the terrain: flatter spots that collect
#' drainage from a large area upslope score higher (wetter), while steep
#' ground or high points with little upslope area draining into them score
#' lower (drier). This is a simple, terrain-only proxy for soil moisture -
#' it doesn't know anything about rainfall, soil type, or vegetation - but
#' it's widely used as a quick first look at where water is likely to
#' pool or drain away. Three closely related formulations are supported to
#' match different conventions used in the literature.
#'
#' @details
#' The standard TWI (Beven & Kirkby 1979) is:
#' \deqn{TWI = \ln\!\left(\frac{a}{\tan\beta}\right)}
#' where \eqn{a} is the specific catchment area (the upslope area draining
#' through each point, per unit of contour length - in effect, "how much
#' water is likely funnelling through here") and \eqn{\beta} is the local
#' slope steepness. A larger contributing area or a flatter slope both push
#' the index up (wetter); a steeper slope pushes it down (drier).
#'
#' **Methods** (the TWI formula itself -- this is a completely separate
#' choice from how the underlying flow accumulation is *routed*, see
#' `route` below):
#' - `"standard"`: the classical formulation above; completely flat cells
#'   (slope = 0) produce `Inf`, since the formula divides by `tan(slope)`.
#' - `"modified"` (default): avoids that `Inf` problem by placing a floor
#'   under how flat a slope is allowed to be (`min_slope`) before doing the
#'   division (Gruber & Peckham 2009).
#' - `"saga"`: the SAGA GIS wetness index (Boehner & Selige 2006), a variant
#'   that estimates the contributing area a little differently and also
#'   applies the same slope floor. Its contributing-area term is also
#'   capped internally at `res * 1e6` (an arbitrarily large but finite
#'   ceiling, in the same units as the specific catchment area itself) so
#'   that an extreme, unbounded accumulation value can't dominate the
#'   result; this only matters for catchments of an unusually large scale.
#'
#' The upslope contributing area is worked out internally via
#' [flowacc()] unless you supply your own already-computed accumulation
#' via `flow_acc`. `route` picks which of [flowacc()]'s routing methods to
#' use for that internal call; left at its default (`NULL`), it
#' auto-selects `"mfd"` when `method = "saga"` (matching SAGA's own
#' literature-typical multiple-flow-direction routing, rather than a
#' plain single-direction D8 path) and `"d8"` for `"standard"`/`"modified"`.
#' Set `route` explicitly to override this for any `method`. `weight` and
#' `bsn` are passed straight through to that internal [flowacc()] call
#' (same meaning as there: a per-cell weight raster, and/or restricting
#' accumulation to stay within each basin of a supplied basin-ID raster).
#' All three (`route`/`weight`/`bsn`) are simply ignored if `flow_acc` is
#' supplied directly, since in that case every routing decision has
#' already been made by whoever computed it.
#'
#' @param dtm A `SpatRaster` of elevation values (metres) - your terrain
#'   model.
#' @param method Character. Which TWI formula to use: `"standard"`,
#'   `"modified"` (default), or `"saga"` - see Details.
#' @param flow_acc Optional `SpatRaster` of flow accumulation, used as a
#'   stand-in for the contributing area `a` above -- typically the output
#'   of [flowacc()]. Expected to be *self-inclusive* (every cell's own
#'   weight already counted, so a headwater cell reads its own weight
#'   rather than 0), matching [flowacc()]'s own convention -- if you
#'   import a flow accumulation layer computed elsewhere (ArcGIS,
#'   whitebox, etc.) that instead follows the older self-excluding
#'   convention (headwater cells read 0), this is detected automatically
#'   (a minimum value of exactly 0 anywhere in the raster) and corrected
#'   for internally, so either convention works without you needing to
#'   know which one your own raster follows. If `NULL`, worked out from
#'   `dtm` automatically via [flowacc()] -- see `route`/`weight`/`bsn`
#'   below.
#' @param route Optional character, one of `"d8"`, `"mfd"`, `"dinf"` --
#'   which [flowacc()] routing method to use when `flow_acc` isn't
#'   supplied. Left at `NULL` (default), auto-selects `"mfd"` for
#'   `method = "saga"` and `"d8"` otherwise -- see Details. Ignored if
#'   `flow_acc` is supplied directly.
#' @param weight Optional `SpatRaster` of per-cell weights, passed
#'   straight through to the internal [flowacc()] call (see its own
#'   `weight` argument). Ignored if `flow_acc` is supplied directly.
#' @param bsn Optional `SpatRaster` of basin IDs, passed straight through
#'   to the internal [flowacc()] call to restrict accumulation to stay
#'   within each basin (see its own `bsn` argument). Ignored if
#'   `flow_acc` is supplied directly.
#' @param slope Optional `SpatRaster` of slope (degrees). If `NULL`, worked
#'   out from `dtm` automatically.
#' @param min_slope Numeric. The minimum slope (degrees) allowed before
#'   dividing by `tan(slope)`, used by the `"modified"` and `"saga"`
#'   methods to avoid producing infinite values on flat ground. Default is
#'   `0.1`.
#' @param fillna Logical. If `TRUE`, `NA` cells in `dtm` (e.g. sea, or any
#'   other no-data gap), and the outer ring of `dtm` itself, are
#'   temporarily filled from their nearest real neighbour before slope is
#'   worked out, and the result is masked back to `dtm`'s own `NA` pattern
#'   afterwards - avoids slope (and the TWI itself) coming out `NA` at
#'   land cells that merely touch an `NA` cell, e.g. an entire coastline's
#'   worth of cells picking up `NA` purely for being next to the sea, and
#'   at the raster's own edge, where the same problem occurs because
#'   there's no neighbouring cell at all just past the boundary. Default
#'   `FALSE`. See [solarindex()]'s own `fillna` for the full rationale.
#'
#' @return A `SpatRaster` of TWI values.
#' @export
#'
#' @seealso [flowacc()], [fillsinks()]
#'
#' @references
#' Beven, K.J. & Kirkby, M.J. (1979). A physically based, variable contributing
#' area model of basin hydrology. *Hydrological Sciences Bulletin*, 24, 43-69.
#'
#' Gruber, S. & Peckham, S. (2009). Land-surface parameters and objects in
#' hydrology. *Developments in Soil Science*, 33, 171-194.
#'
#' Boehner, J. & Selige, T. (2006). Spatial prediction of soil attributes using
#' terrain analysis and climate regionalisation. *SAGA - Analysis and Modelling
#' Applications*, 115, 13-28.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' dtm_filled <- fillsinks(dtm)
#' w   <- twi(dtm_filled, method = "modified")
#' plot(w)
#'
#' # saga formula, auto-routed via mfd
#' w_saga <- twi(dtm_filled, method = "saga")
#'
#' # force D-infinity routing regardless of method
#' w_dinf <- twi(dtm_filled, method = "modified", route = "dinf")
#'
#' # DTM with NA cells at a coastline: avoid an NA halo on land
#' w_coast <- twi(dtm_filled, method = "modified", fillna = TRUE)
#' }
twi <- function(dtm, method = "modified", flow_acc = NULL, route = NULL,
                 weight = NULL, bsn = NULL, slope = NULL, min_slope = 0.1,
                 fillna = FALSE) {
  stopifnot(inherits(dtm, "SpatRaster"))
  method <- match.arg(method, c("standard", "modified", "saga"))
  if (!is.null(route)) route <- match.arg(route, c("d8", "mfd", "dinf"))

  # ---- Optional NA fill (see `fillna` below / solarindex()'s own version
  # of this for the full rationale) -- NA cells in `dtm`, and the outer
  # ring of `dtm` itself (whose focal window reaches past the raster's own
  # edge), are temporarily filled from their nearest real neighbour via
  # .fillna_dtm() so terra::terrain()'s slope calculation doesn't come out
  # NA at every valid land cell touching one (e.g. a coastline) or at the
  # raster's edge, then the result is masked back to `dtm`'s real NA
  # pattern below.
  dtm_na <- dtm
  fill_ext <- NULL
  if (isTRUE(fillna)) {
    fill_ext <- .fillna_dtm(dtm)
    dtm <- fill_ext$interior
  }

  res <- terra::res(dtm)[1]

  # When fillna = TRUE, computed on fill_ext$extended (1 cell wider than
  # `dtm` in every direction) so the outer ring gets a real value too,
  # then cropped back to `dtm`'s own extent.
  if (is.null(slope)) {
    slope <- if (isTRUE(fillna)) {
      terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm_na)
    } else terra::terrain(dtm, "slope", unit = "degrees")
  }

  if (is.null(flow_acc)) {
    # No precomputed accumulation supplied -- work it out via flowacc().
    # `route` picks the method; left NULL, auto-select a sensible default
    # per TWI `method` (see @details above for the reasoning).
    route_use <- route
    if (is.null(route_use)) route_use <- if (method == "saga") "mfd" else "d8"
    flow_acc <- flowacc(dtm, method = route_use, weight = weight, bsn = bsn)
  }

  fa_m  <- as.matrix(flow_acc, wide = TRUE)
  sl_m  <- as.matrix(slope,    wide = TRUE)

  # twi_cpp() expects flow_acc to be self-inclusive (see its own comment
  # in src/wetness.cpp) -- true of flowacc()'s own output, but a user
  # importing a flow accumulation layer from elsewhere (ArcGIS, whitebox,
  # etc.) may well be handing in the older self-EXCLUDING convention,
  # where a headwater cell reads 0 rather than 1. A minimum of exactly 0
  # anywhere in the raster is the tell for that -- flowacc()'s own output
  # can never be exactly 0 under its default weight=1, so this is a safe,
  # low-false-positive way to detect and silently correct for it, rather
  # than requiring the caller to know or remember which convention their
  # own raster follows.
  if (isTRUE(min(fa_m, na.rm = TRUE) == 0)) fa_m <- fa_m + 1

  result <- twi_cpp(fa_m, sl_m, res, method, min_slope)

  out <- terra::rast(dtm)
  terra::values(out) <- result
  names(out) <- paste0("twi_", method)
  if (isTRUE(fillna)) out <- terra::mask(out, dtm_na)
  out
}

