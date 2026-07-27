# Data download functions for terravars.
# All network access lives in this file.

# Web Mercator "sphere" radius used by the XYZ/slippy-map tile scheme
# (the spherical approximation the whole tile ecosystem is built on, not
# the WGS84 ellipsoid) -- metres.
.WEBMERC_R <- 6378137

# Reprojects r's extent to lon/lat (EPSG:4326), via its 4 corner points
# built as a SpatVector and run through terra::project()'s SpatVector
# method -- the most well-established form of project(). Errors clearly
# if the result isn't usable, rather than letting a stray NA silently
# propagate into the tile-index maths downstream.
.bbox_lonlat <- function(r) {
  e <- as.vector(terra::ext(r))
  pts <- terra::vect(cbind(x = c(e["xmin"], e["xmax"], e["xmin"], e["xmax"]),
                            y = c(e["ymin"], e["ymin"], e["ymax"], e["ymax"])),
                      crs = terra::crs(r))
  bb <- as.vector(terra::ext(terra::project(pts, "EPSG:4326")))
  if (anyNA(bb))
    stop("Could not reproject 'r's extent to longitude/latitude -- check ",
         "that 'r' has a valid coordinate reference system.", call. = FALSE)
  bb
}

# The single AWS "elevation-tiles-prod" zoom level that guarantees ground
# resolution at least as fine as `r`'s own resolution everywhere within
# `r`'s extent. Ground resolution at a fixed zoom varies with latitude
# (a cos(lat) term -- tiles cover less real ground the further from the
# equator), so the worst case (finest zoom needed) occurs wherever `r`'s
# extent comes closest to the equator -- found from `r`'s reprojected
# extent (see .bbox_lonlat()), cheap compared to reprojecting every cell.
# This deliberately fetches one zoom for the whole of `r` rather than a
# spatially-varying zoom per pixel: simpler to implement and verify
# correctly, at the cost of possibly fetching more resolution than
# strictly necessary wherever `r` extends well away from its
# closest-to-equator point. Capped at 14, the finest zoom at which the
# source tiles add new data. Tile size is 512px for the
# elevation-tiles-prod GeoTIFF tiles .aws_dem_mosaic() below fetches
# (confirmed from an actual downloaded tile's dimensions -- this
# differs from the 256px tile size of the same bucket's Terrarium PNG
# tiles, so it must be kept in step with whichever format
# .aws_dem_mosaic() is actually using).
.aws_dem_zoom <- function(r) {
  reso <- min(terra::res(r))
  e_ll <- .bbox_lonlat(r)
  lat_range <- c(e_ll["ymin"], e_ll["ymax"])
  lat0 <- if (lat_range[1] <= 0 && lat_range[2] >= 0) 0 else lat_range[which.min(abs(lat_range))]
  z <- ceiling(log2((cos(lat0 * pi / 180) * 2 * pi * .WEBMERC_R) / (512 * reso)))
  min(max(z, 0), 14)
}

# Downloads and mosaics every AWS "elevation-tiles-prod" GeoTIFF tile
# needed to cover `r`'s extent at zoom `z` -- the same
# ".../elevation-tiles-prod/geotiff/{z}/{x}/{y}.tif" endpoint elevatr's
# own get_aws_terrain() fetches from (as opposed to the same bucket's
# Terrarium PNG tiles, which need RGB decoding and manually-assigned
# extent/CRS). GeoTIFF tiles are self-georeferencing -- terra::rast()
# reads each tile's real extent and CRS (EPSG:3857, the tile scheme's
# native CRS) straight from the file, so there's no hand-rolled tile
# extent maths here for a bug to hide in; reprojecting to match `r` is
# left to the caller. A tile that fails to download or read is skipped
# with a warning rather than aborting the whole call. Returns NULL (with
# a message already issued) if not a single tile could be retrieved.
.aws_dem_mosaic <- function(r, z) {
  e_ll <- .bbox_lonlat(r)
  n <- 2^z

  # Web Mercator's own valid latitude range (beyond this the projection
  # is undefined) -- clamp before the tan()/log() tile-index formula.
  .clamp_lat <- function(lat) max(min(lat, 85.0511), -85.0511)

  # unname() lon/lat before building the c(x=, y=) result below -- lon/lat
  # arrive as single named elements (e.g. e_ll["xmax"], named "xmax"), and
  # c(x = <value already named "xmax">) concatenates both names into
  # "x.xmax" rather than plain "x". tl["x"]/br["x"] would then silently
  # return NA instead of erroring, exactly the bug that produced
  # `max(0, tl["x"]):min(n-1, br["x"]) : NA/NaN argument` here.
  tile_xy <- function(lon, lat) {
    lon <- unname(lon)
    lat <- unname(.clamp_lat(lat))
    c(x = floor((lon + 180) / 360 * n),
      y = floor((1 - log(tan(lat * pi / 180) + 1 / cos(lat * pi / 180)) / pi) / 2 * n))
  }

  tl <- tile_xy(e_ll["xmin"], e_ll["ymax"])  # smallest x, smallest y
  br <- tile_xy(e_ll["xmax"], e_ll["ymin"])  # largest x, largest y
  xs <- max(0, tl["x"]):min(n - 1, br["x"])
  ys <- max(0, tl["y"]):min(n - 1, br["y"])

  ntiles <- length(xs) * length(ys)
  if (ntiles > 500)
    warning(ntiles, " elevation tiles are needed to cover this area at zoom ",
            z, " -- this may take a while and use significant bandwidth.",
            call. = FALSE)

  tiles <- list()
  for (x in xs) {
    for (y in ys) {
      url <- sprintf("https://s3.amazonaws.com/elevation-tiles-prod/geotiff/%d/%d/%d.tif", z, x, y)
      f <- tempfile(fileext = ".tif")
      ok <- tryCatch({
        utils::download.file(url, f, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)

      if (ok && file.exists(f)) {
        # Re-assigning a SpatRaster's own values to itself forces them
        # into memory rather than leaving the object pointing at the
        # tempfile on disk -- necessary since that tempfile is
        # unlink()ed right below, before the mosaic step that actually
        # needs the data. (terra also has a function for this, but its
        # name has changed across versions -- readAll() in some,
        # toMemory() in others -- so this values<- idiom is used
        # instead, since it has worked unchanged across all terra
        # versions.)
        elev <- tryCatch({
          e <- terra::rast(f)
          terra::values(e) <- terra::values(e)
          e
        }, error = function(e) NULL)
        if (!is.null(elev)) {
          tiles[[length(tiles) + 1]] <- elev
        } else {
          warning("Tile z=", z, " x=", x, " y=", y, " could not be read; skipped.", call. = FALSE)
        }
      } else {
        warning("Tile z=", z, " x=", x, " y=", y, " could not be downloaded; skipped.", call. = FALSE)
      }
      unlink(f)
    }
  }

  if (length(tiles) == 0) return(NULL)
  if (length(tiles) == 1) return(tiles[[1]])
  terra::mosaic(terra::sprc(tiles), fun = "first")
}

# Downloads and mosaics every ESA WorldCover 10 m land-cover GeoTIFF
# tile needed to cover `r`'s extent, from the public, unauthenticated
# "esa-worldcover" Amazon S3 bucket (part of the AWS Registry of Open
# Data). Tiles are 3x3 degree, self-georeferencing (EPSG:4326) GeoTIFFs
# named by the lat/lon of their lower-left corner -- e.g. "S48E036" for
# the tile spanning 45S-48S, 36E-39E -- computed directly here from
# `r`'s extent (floor()ed to the nearest 3 degrees) rather than looked
# up via the bucket's own grid index geojson, the same way the DEM
# tiles' z/x/y indices are computed directly in .aws_dem_mosaic() above
# rather than fetched: avoids a second network round-trip and an sf
# dependency just to look up tile names. Each tile is cropped to `r`'s
# (buffered) extent immediately -- terra::crop() on a file-backed
# SpatRaster reads only the needed window, so this stays cheap even
# though a full tile is enormous (3x3 degrees at 10 m is roughly
# 33000 x 33000 cells). A tile that fails to download or read is
# skipped with a warning rather than aborting the whole call -- this
# is also how a genuinely all-ocean tile (not present in the bucket at
# all) is handled. Returns NULL (with a message already issued by the
# caller) if not a single tile could be retrieved.
.worldcover_mosaic <- function(r, year) {
  version <- if (year == 2020) "v100" else "v200"
  base_url <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com"

  e_ll <- .bbox_lonlat(r)
  buf <- 0.05
  lats <- seq(floor((e_ll["ymin"] - buf) / 3) * 3, floor((e_ll["ymax"] + buf) / 3) * 3, by = 3)
  lons <- seq(floor((e_ll["xmin"] - buf) / 3) * 3, floor((e_ll["xmax"] + buf) / 3) * 3, by = 3)

  ntiles <- length(lats) * length(lons)
  if (ntiles > 20)
    warning(ntiles, " ESA WorldCover tiles (each 3x3 degrees, ~33000x33000 ",
            "cells) are needed to cover this area -- this may take a while ",
            "and use significant bandwidth.", call. = FALSE)

  tile_code <- function(lat, lon) {
    sprintf("%s%02d%s%03d",
            if (lat < 0) "S" else "N", abs(lat),
            if (lon < 0) "W" else "E", abs(lon))
  }

  tiles <- list()
  for (lat in lats) {
    for (lon in lons) {
      code <- tile_code(lat, lon)
      url <- sprintf("%s/%s/%d/map/ESA_WorldCover_10m_%d_%s_%s_Map.tif",
                      base_url, version, year, year, version, code)
      f <- tempfile(fileext = ".tif")
      ok <- tryCatch({
        utils::download.file(url, f, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)

      if (ok && file.exists(f)) {
        lc <- tryCatch({
          x <- terra::rast(f)
          x <- terra::crop(x, terra::ext(e_ll["xmin"] - buf, e_ll["xmax"] + buf,
                                          e_ll["ymin"] - buf, e_ll["ymax"] + buf))
          # Force into memory before the tempfile is unlink()ed below --
          # see the identical comment in .aws_dem_mosaic() above for why
          # this values<- idiom is used rather than a named terra function.
          terra::values(x) <- terra::values(x)
          x
        }, error = function(e) NULL)
        if (!is.null(lc)) {
          tiles[[length(tiles) + 1]] <- lc
        } else {
          warning("WorldCover tile ", code, " could not be read; skipped.", call. = FALSE)
        }
      } else {
        warning("WorldCover tile ", code, " could not be downloaded (expected ",
                "for a tile that's entirely ocean, which the bucket doesn't ",
                "provide); skipped.", call. = FALSE)
      }
      unlink(f)
    }
  }

  if (length(tiles) == 0) return(NULL)
  if (length(tiles) == 1) return(tiles[[1]])
  terra::mosaic(terra::sprc(tiles), fun = "first")
}

#' @title Download ESA WorldCover Land Cover Data
#'
#' @description Downloads ESA WorldCover 10 m resolution global land-cover
#' data directly from its public Amazon Web Services (AWS) home - no
#' account, API key, or additional package required.
#' @param r a `SpatRaster` object used as a template for downloading data,
#'   with a coordinate reference system defined.
#' @param year optional numeric, either `2020` or `2021` (default),
#'   selecting which WorldCover release to fetch (`2020` = v100, `2021`
#'   = v200 - the more recent release).
#' @param watersplit optional logical (default `TRUE`) splitting code
#'   `80` ("Permanent water bodies") into two: `80` (sea) and `81`
#'   (freshwater) - see Details.
#' @return a `SpatRaster` of land-cover class codes, on the same grid
#'   (resolution, extent, coordinate reference system) as `r`.
#' @details The data returned matches the resolution, extent and
#'   coordinate reference system of `r`. Because land-cover class codes
#'   are categorical, resampling onto `r`'s grid uses modal ("mode")
#'   resampling - whichever class is most common within each output
#'   cell - rather than the interpolation used for continuous data like
#'   elevation.
#'
#'   Data is fetched directly from the public, unauthenticated
#'   `esa-worldcover` Amazon S3 bucket (part of the AWS Registry of Open
#'   Data - see https://registry.opendata.aws/esa-worldcover-vito/) as
#'   3x3 degree GeoTIFF tiles and mosaicked internally.
#'
#'   Returned values follow the ESA WorldCover 10 m classification:
#'
#'   | Code | Class | Code | Class |
#'   |---|---|---|---|
#'   | 10 | Tree cover | 70 | Snow and ice |
#'   | 20 | Shrubland | 80 | Permanent water bodies (sea) |
#'   | 30 | Grassland | 81 | Permanent water bodies (freshwater) |
#'   | 40 | Cropland | 90 | Herbaceous wetland |
#'   | 50 | Built-up | 95 | Mangroves |
#'   | 60 | Bare / sparse vegetation | 100 | Moss and lichen |
#'
#'   ESA WorldCover's own code `80` doesn't distinguish sea from a large
#'   inland lake or river - both are "Permanent water bodies", per its
#'   Product User Manual ("lakes, reservoirs, and rivers... either
#'   fresh or salt-water bodies"). When `watersplit = TRUE` (default),
#'   code `80` is split into two after fetching: connected water
#'   patches (8-connectivity) that touch `r`'s own edge are presumed
#'   connected to open water beyond this extent and kept as `80` (sea);
#'   patches fully enclosed within the interior, never touching the
#'   edge, become `81` (freshwater) - a new code, not part of ESA's own
#'   scheme, chosen from the unused gap between `80` and `90` so it
#'   can't collide with a real WorldCover class. This is a heuristic
#'   based purely on whether a water body touches `r`'s edge: a lake
#'   that happens to sit right at the edge of `r`'s extent (but isn't
#'   actually the sea) would be called `80`, and a tidal inlet fully
#'   enclosed within a large `r` would be called `81`. Set
#'   `watersplit = FALSE` to keep ESA's original undivided code `80`
#'   for all water.
#'
#'   The result is a natural fit for [dem_download()]'s `mask` argument
#'   (e.g. `mask_value = 80` to mask out sea only, leaving inland water
#'   bodies alone) - a 10 m land-cover classification distinguishes
#'   land from water at far finer detail than a generalised global
#'   coastline vector.
#' @source ESA WorldCover 10 m 2020 v100, Zenodo
#'   https://doi.org/10.5281/zenodo.5571936; ESA WorldCover 10 m 2021
#'   v200, Zenodo https://doi.org/10.5281/zenodo.7254221.
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(dtm100m)
#' lc <- landcover_download(r)
#' plot(lc)
#'
#' # Keep ESA's original undivided water code
#' lc_raw <- landcover_download(r, watersplit = FALSE)
#' }
#' @export
landcover_download <- function(r, year = 2021, watersplit = TRUE) {
  stopifnot(inherits(r, "SpatRaster"))
  stopifnot(is.numeric(year), length(year) == 1L, year %in% c(2020, 2021))
  stopifnot(is.logical(watersplit), length(watersplit) == 1L, !is.na(watersplit))

  if (!curl::has_internet()) {
    message("Please connect to the internet and try again.")
    return(NULL)
  }

  mosaic <- .worldcover_mosaic(r, year)
  if (is.null(mosaic)) {
    message("No land-cover tiles could be downloaded for this area.")
    return(NULL)
  }

  lc <- terra::project(mosaic, r, method = "mode")

  if (isTRUE(watersplit)) {
    wmask <- terra::ifel(!is.na(lc) & (lc == 80), 1, NA)
    n_water <- terra::global(wmask, "sum", na.rm = TRUE)[1, 1]
    if (!is.na(n_water) && n_water > 0) {
      # Label each connected water patch (8-connectivity, matching every
      # other neighbourhood search in this package -- see fillsinks()'s
      # own waterbody handling for the same pattern), then split code 80
      # into sea (touches r's own edge) vs freshwater (doesn't) -- see
      # the roxygen block above for the reasoning and its limits.
      patches_r <- terra::patches(wmask, directions = 8)
      pm <- .is(patches_r)
      nr <- nrow(pm); nc <- ncol(pm)
      border_ids <- unique(c(pm[1, ], pm[nr, ], pm[, 1], pm[, nc]))
      border_ids <- border_ids[!is.na(border_ids)]
      on_border <- matrix(pm %in% border_ids, nrow = nr, ncol = nc)
      fresh <- !is.na(pm) & !on_border

      lcm <- .is(lc)
      lcm[fresh] <- 81
      lc <- .rast(lcm, lc)
    }
  }

  lc
}

#' @title Download Digital Elevation Data
#'
#' @description Downloads Digital Elevation Data directly from the Amazon
#' Web Services (AWS) public "Terrain Tiles" dataset - no account, API
#' key, or additional elevation-data package required.
#' @param r a `SpatRaster` object used as a template for downloading data,
#'   with a coordinate reference system defined.
#' @param mask Optional `SpatRaster`, on the same grid as `r`, used to
#'   mask cells out of the downloaded elevation data (see Details) --
#'   this can be anything, e.g. a land-cover raster from
#'   [landcover_download()] (to mask out water), or a mask of your own.
#'   Not supplied (`NULL`, the default) leaves the raw downloaded
#'   elevation data untouched.
#' @param mask_value Which value(s) in `mask` mark a cell to be masked
#'   out. Default `NA`: cells where `mask` is `NA` are masked out. Set
#'   to a specific value or vector of values if `mask` uses a
#'   real-valued classification instead - e.g. `80` for sea when `mask`
#'   is a [landcover_download()] raster with its default
#'   `watersplit = TRUE` (its code `81`, freshwater, is then left
#'   alone; pass `c(80, 81)` to mask out all water instead). Ignored if
#'   `mask` isn't supplied. Follows the same convention as
#'   [fillsinks()]'s own `waterbody`/`water_value` arguments.
#' @param mask_out Either `"NA"` (default) or `"0"`, controlling what
#'   masked-out cells become in the returned raster: `NA`, or `0` (sea
#'   level). Ignored if `mask` isn't supplied.
#' @param belowsea optional logical (default `FALSE`). When `FALSE`, any
#'   cell in the returned raster still below 0 (sea level) is clamped to
#'   `0` -- see Details. Set to `TRUE` to keep genuine below-sea-level
#'   land (e.g. a polder) unmodified.
#' @return a `SpatRaster` of digital elevation data, on the same grid
#'   (resolution, extent, coordinate reference system) as `r`.
#' @details The data returned matches the resolution, extent and
#'   coordinate reference system of `r`. The resolution of the underlying
#'   source tiles varies globally, so in some instances the data returned
#'   after resampling to match `r` may be derived by interpolation rather
#'   than reflecting the tiles' own native resolution. A map showing the
#'   resolution of available data is at
#'   https://github.com/tilezen/joerd/blob/master/docs/images/footprints-preview.png.
#'
#'   The elevation data available from AWS blends digital elevation and
#'   bathymetry data, the latter generally only available at coarser
#'   resolution, which frequently leaves coastal land cells with
#'   spurious below-sea-level (negative) values instead of the small
#'   positive elevations they should have. Two related fixes are applied
#'   for this. First, whenever `mask` is supplied, any cell that is not
#'   itself masked out but is immediately adjacent (8-connectivity) to a
#'   cell that is masked out -- i.e. a land cell right at the coast --
#'   is clamped to `0` if it comes out negative. This runs automatically
#'   whenever `mask` is given, regardless of `belowsea`, since a
#'   negative value immediately next to the sea is essentially always
#'   this artefact rather than real below-sea-level terrain. Second,
#'   `belowsea` controls a blanket version of the same fix applied to
#'   the *whole* returned raster rather than just cells next to `mask`:
#'   by default (`belowsea = FALSE`), every remaining cell below `0` is
#'   also clamped to `0`. Set `belowsea = TRUE` to skip this blanket
#'   step and keep genuinely below-sea-level land elsewhere in `r`
#'   (a polder, for example) at its true elevation -- the coastal fix
#'   above still runs in that case if `mask` is supplied.
#'
#'   Data is fetched directly from the public, unauthenticated
#'   `elevation-tiles-prod` Amazon S3 bucket (part of the AWS Registry of
#'   Open Data - see https://registry.opendata.aws/terrain-tiles/) as
#'   GeoTIFF tiles (the same source and format the `elevatr` package's
#'   own `get_aws_terrain()` uses) and mosaicked internally. The zoom
#'   level used is whichever single level guarantees resolution at
#'   least as fine as `r`'s own resolution everywhere within `r`'s extent
#'   (capped at zoom 14, the finest zoom at which the source tiles add
#'   new data) - this can mean more resolution than strictly necessary is
#'   fetched for an `r` spanning a wide range of latitudes, trading some
#'   extra download size for a simpler, single-zoom mosaic.
#'
#'   Supplying `mask` (e.g. a [landcover_download()] raster, at 10 m
#'   resolution - far finer than a generalised coastline vector, with
#'   `mask_value = 80` to mask out its "Permanent water bodies" class)
#'   also lets sea cells be masked out entirely instead of relying on
#'   the elevation data's own land/sea classification. This function does
#'   not fetch `mask` itself; call [landcover_download()] (or build
#'   your own mask) and pass its result in explicitly, so the same
#'   fetch can be reused across multiple `dem_download()` calls rather
#'   than repeated on every one.
#' @source Terrain Tiles was accessed from
#'   https://registry.opendata.aws/terrain-tiles.
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(dtm100m)
#' dtm <- dem_download(r)
#' plot(dtm)
#'
#' # Mask out water using a land-cover raster fetched separately, so it
#' # can be reused across more than one dem_download() call
#' lc <- landcover_download(r)
#' dtm_clean <- dem_download(r, mask = lc, mask_value = 80)
#'
#' # Same, but sea level (0) instead of NA
#' dtm_clean0 <- dem_download(r, mask = lc, mask_value = 80, mask_out = "0")
#'
#' # Keep genuine below-sea-level land away from the coast unmodified
#' dtm_polder <- dem_download(r, mask = lc, mask_value = 80, belowsea = TRUE)
#' }
#' @export
dem_download <- function(r, mask = NULL, mask_value = NA, mask_out = c("NA", "0"),
                          belowsea = FALSE) {
  stopifnot(inherits(r, "SpatRaster"))
  mask_out <- match.arg(mask_out)
  stopifnot(is.logical(belowsea), length(belowsea) == 1L, !is.na(belowsea))
  if (!is.null(mask)) {
    stopifnot(inherits(mask, "SpatRaster"))
    if (!identical(dim(.is(mask)), dim(.is(r))))
      stop("'mask' must be on the same grid as 'r'.", call. = FALSE)
  }

  if (!curl::has_internet()) {
    message("Please connect to the internet and try again.")
    return(NULL)
  }

  z <- .aws_dem_zoom(r)
  mosaic <- .aws_dem_mosaic(r, z)
  if (is.null(mosaic)) {
    message("No elevation tiles could be downloaded for this area.")
    return(NULL)
  }

  dtm <- terra::project(mosaic, r)

  if (!is.null(mask)) {
    # Same NA-flag + explicit-value(s) mask-building convention as
    # fillsinks()'s own 'waterbody'/'water_value' handling -- built from
    # is.na()/==/"|" only (not %in%) so this works identically regardless
    # of terra version.
    na_flag  <- is.na(mask_value)
    explicit <- mask_value[!na_flag]
    excl <- NULL
    if (any(na_flag)) excl <- is.na(mask)
    for (v in explicit) {
      m <- !is.na(mask) & (mask == v)
      excl <- if (is.null(excl)) m else (excl | m)
    }

    # Coastal bathymetry-blend fix: any cell not itself masked out but
    # directly touching (8-connectivity) a masked-out cell -- i.e. right
    # at the coast -- is clamped to 0 if it's come out negative. Dilating
    # 'excl' by one cell with a 3x3 max focal window is the same
    # adjacency test as the border-detection logic in
    # landcover_download()'s watersplit, just done with a neighbourhood
    # filter rather than terra::patches() since only direct adjacency,
    # not connected-component membership, is needed here. This runs
    # regardless of 'belowsea', which governs cells away from 'mask'.
    excl_num <- terra::ifel(excl, 1, 0)
    touches_masked <- terra::focal(excl_num, w = 3, fun = "max", na.rm = TRUE) == 1
    coastal_neg <- touches_masked & !excl & !is.na(dtm) & (dtm < 0)
    dtm <- terra::ifel(coastal_neg, 0, dtm)

    fillval <- if (mask_out == "0") 0 else NA
    dtm <- terra::ifel(excl, fillval, dtm)
  }

  if (!isTRUE(belowsea)) {
    dtm <- terra::ifel(!is.na(dtm) & (dtm < 0), 0, dtm)
  }

  dtm
}
