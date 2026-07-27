# ============================================================
# data.R
#
# Documentation for package example datasets (Lizard peninsula, Cornwall,
# UK). Every dataset here shares the same grid (resolution, extent, CRS)
# as `dtm100m`, the example DTM used throughout this package's examples.
#
# Stored as *wrapped* `SpatRaster` objects (`terra::wrap()`) so they can
# be serialised to `data/*.rda` -- a plain `SpatRaster` holds an external
# C++ pointer that cannot survive `save()`/`load()`. Unwrap with
# `terra::rast()` (or `terra::unwrap()`) before use; see each dataset's
# examples below.
# ============================================================

#' Digital terrain model (Lizard, Cornwall)
#'
#' Example digital terrain model (DTM, elevation in metres) raster for the
#' Lizard peninsula, Cornwall, UK - the primary example dataset used
#' throughout this package's documentation and examples (every function
#' that takes a `dtm` argument uses this raster in its `@examples`).
#' Defines the grid (resolution, extent, CRS) that every other bundled
#' dataset ([pai], [hgt], [vegx], [landcover]) shares. `NA` cells (e.g.
#' open sea beyond the coastline) follow this package's usual land/sea
#' convention, used directly by e.g. [fillsinks()] and [coastalexposure()].
#'
#' @format A wrapped `SpatRaster` (`terra::PackedSpatRaster`), single
#'   layer, 100 m resolution, projected CRS (transverse Mercator, airy
#'   ellipsoid - OSGB36 / British National Grid). Unwrap with
#'   `terra::rast(dtm100m)` before use.
#' @source Elevation data for the Lizard peninsula, Cornwall, UK,
#'   resampled to this package's own 100 m example grid.
#' @seealso [pai], [hgt], [vegx], [landcover], [dtmc]
#' @examples
#' \dontrun{
#' library(terra)
#' dtm <- rast(dtm100m)
#' plot(dtm)
#' }
"dtm100m"

#' Plant area index (Lizard, Cornwall)
#'
#' Example plant area index (PAI, m^2 m^-2) raster - a measure of how much
#' leaf, stem and branch material sits above each point in the landscape,
#' higher values meaning denser vegetation - for the Lizard peninsula,
#' Cornwall, UK. For use with the `pai` argument of [solarindex()],
#' [skyview()], [twostream()] and [canopy3d()]. Same grid (resolution,
#' extent, CRS) as `dtm100m`.
#'
#' Derived from Sentinel-based leaf area index (LAI), obtained via WEkEO
#' (the Copernicus Data and Information Access Service); used here as an
#' approximation of plant area index (PAI includes non-leaf plant material
#' such as stems and branches, which LAI alone does not capture).
#'
#' @format A wrapped `SpatRaster` (`terra::PackedSpatRaster`), single
#'   layer, 100 m resolution. Unwrap with `terra::rast(pai)` before use.
#' @source Leaf area index derived from Sentinel data via WEkEO
#'   (<https://www.wekeo.eu>).
#' @seealso [hgt], [vegx], [landcover], [solarindex()], [skyview()], [twostream()], [canopy3d()]
#' @examples
#' \dontrun{
#' library(terra)
#' pai_r <- rast(pai)
#' plot(pai_r)
#' }
"pai"

#' Canopy height (Lizard, Cornwall)
#'
#' Example canopy height (metres above local ground) raster for the Lizard
#' peninsula, Cornwall, UK - how tall the vegetation is at each point. For
#' use with the `hgt` argument of [solarindex()], [skyview()] and
#' [canopy3d()], which switches those functions on to a more realistic
#' calculation that accounts for the actual 3D shape of nearby vegetation
#' rather than assuming it forms a uniform layer. Same grid (resolution,
#' extent, CRS) as `dtm100m`.
#'
#' Derived from the global canopy height model of Lang et al. (2023).
#'
#' @format A wrapped `SpatRaster` (`terra::PackedSpatRaster`), single
#'   layer, 100 m resolution. Unwrap with `terra::rast(hgt)` before use.
#' @source Lang, N., Jetz, W., Schindler, K., & Wegner, J.D. (2023). A
#'   high-resolution canopy height model of the Earth. *Nature Ecology &
#'   Evolution*, 1--12.
#' @seealso [pai], [vegx], [landcover], [solarindex()], [skyview()], [canopy3d()]
#' @examples
#' \dontrun{
#' library(terra)
#' hgt_r <- rast(hgt)
#' plot(hgt_r)
#' }
"hgt"

#' Leaf angle distribution parameter (Lizard, Cornwall)
#'
#' Example leaf angle distribution raster for the Lizard peninsula,
#' Cornwall, UK, for use with the `x` argument of [solarindex()],
#' [skyview()], [twostream()] and [canopy3d()]. This describes how the
#' leaves in each point's canopy are typically angled, which affects how
#' much sunlight passing through gets blocked versus slipping through gaps:
#' values close to 0 mean canopies of mostly upright/vertical leaves (e.g.
#' many grasses and some conifers), 1 means a "typical" mixed canopy with
#' leaves at all angles, and large values (up to `Inf`) mean canopies of
#' mostly flat, horizontal leaves (e.g. many broadleaf trees). (Technically,
#' the ratio of leaf area projected onto a horizontal plane versus a
#' vertical one, following the Campbell & Norman (1998) ellipsoidal
#' convention.) Same grid (resolution, extent, CRS) as `dtm100m`.
#'
#' Assigned per land-cover class via a lookup table applied to [landcover],
#' this package's own 10 m resolution European Space Agency land cover
#' dataset (resampled to this same grid).
#'
#' @format A wrapped `SpatRaster` (`terra::PackedSpatRaster`), single
#'   layer, 100 m resolution. Unwrap with `terra::rast(vegx)` before use.
#' @source Derived from [landcover] via a per-class leaf angle
#'   distribution lookup table.
#' @seealso [pai], [hgt], [landcover], [solarindex()], [skyview()], [twostream()], [canopy3d()]
#' @examples
#' \dontrun{
#' library(terra)
#' vegx_r <- rast(vegx)
#' plot(vegx_r)
#' }
"vegx"

#' Habitat / land cover classification (Lizard, Cornwall)
#'
#' Example habitat/land cover raster for the Lizard peninsula, Cornwall,
#' UK, using the 21-class habitat scheme in the table below. This is the
#' land cover data [vegx]'s leaf angle distribution lookup table was built
#' from. Same grid (resolution, extent, CRS) as `dtm100m`. Classes 13
#' (Saltwater) and 14 (Freshwater) make this a ready-made example for
#' [fillsinks()]'s `waterbody` argument (e.g. `water_value = c(13, 14)`).
#'
#' Source data is the European Space Agency's 10 m resolution land cover
#' product, resampled to `dtm100m`'s 100 m grid using modal ("mode")
#' resampling -- the appropriate choice for categorical data, since it
#' keeps whichever class is most common within each 100 m cell rather than
#' blending neighbouring classes together the way a numeric resampling
#' method would -- reprojected/aligned to `dtm100m`'s own CRS and extent,
#' and masked to match it (`NA` wherever `dtm100m` is `NA`, e.g. open sea
#' beyond the coastline, this package's usual convention).
#'
#' | Code | Habitat                 |
#' |-----:|:-------------------------|
#' | 1    | Broadleaved woodland     |
#' | 2    | Coniferous woodland      |
#' | 3    | Arable and horticulture  |
#' | 4    | Improved grassland       |
#' | 5    | Neutral grassland        |
#' | 6    | Calcareous grassland     |
#' | 7    | Acid grassland           |
#' | 8    | Fen, marsh and swamp     |
#' | 9    | Heather                  |
#' | 10   | Heather grassland        |
#' | 11   | Bog                      |
#' | 12   | Inland rock              |
#' | 13   | Saltwater                |
#' | 14   | Freshwater               |
#' | 15   | Supralittoral rock       |
#' | 16   | Supralittoral sediment   |
#' | 17   | Littoral rock            |
#' | 18   | Littoral sediment        |
#' | 19   | Saltmarsh                |
#' | 20   | Urban                    |
#' | 21   | Suburban                 |
#'
#' @format A wrapped `SpatRaster` (`terra::PackedSpatRaster`), single
#'   layer, 100 m resolution, integer class codes 1--21 (see Details for
#'   the full lookup table). Unwrap with `terra::rast(landcover)` before use.
#' @source European Space Agency 10 m resolution land cover, resampled
#'   (modal) to `dtm100m`'s grid.
#' @seealso [pai], [hgt], [vegx], [fillsinks()]
#' @examples
#' \dontrun{
#' library(terra)
#' landcover_r <- rast(landcover)
#' plot(landcover_r)
#'
#' # Use as a water body mask for fillsinks() -- Saltwater and Freshwater
#' # are classes 13 and 14
#' dtm <- rast(dtm100m)
#' dtm_filled <- fillsinks(dtm, waterbody = landcover_r, water_value = c(13, 14))
#' }
"landcover"

#' Coarser regional elevation, for extending coastal searches (Lizard, Cornwall region)
#'
#' A list of two coarser, larger-extent elevation rasters surrounding `dtm100m`'s own area, for use
#' with the `coarse` argument of [coastalexposure()] -- lets a search based on `dtm100m`'s fine 100 m
#' grid reach much further out to sea (or further inland) than `dtm100m`'s own extent would otherwise
#' allow, without needing one enormous fine-resolution raster. `NA` cells represent sea and any
#' non-`NA` value (elevation, metres) represents land, matching [coastalexposure()]'s
#' `landsea`/`coarse` convention.
#'
#' \describe{
#'   \item{`dtm1`}{500 m resolution, covering `dtm100m`'s own extent buffered outward by 50 km on
#'     every side.}
#'   \item{`dtm2`}{2500 m resolution, covering `dtm1`'s extent buffered outward by a further 500 km
#'     on every side.}
#' }
#'
#' Built by downloading elevation data onto two successively coarser, larger template grids, each
#' produced by aggregating the previous level five-fold (100 m -> 500 m -> 2500 m) and then extending
#' its extent outward, so each level's resolution and extent were fixed by its own template grid
#' before any data was fetched for it. Land/sea masking was applied afterward using Natural Earth's
#' 1:10m physical land polygons, reprojected to `dtm100m`'s own CRS; both layers' CRS was then
#' explicitly set to match `dtm100m` exactly, rather than relying on whatever CRS the reprojection
#' step itself returned, avoiding spurious CRS-mismatch issues downstream even though the projection
#' is unchanged either way.
#'
#' @format A list of two wrapped `SpatRaster` objects (`terra::PackedSpatRaster`), named `dtm1` and
#'   `dtm2` (see Details for each one's resolution and extent). Unwrap each with `terra::rast()`
#'   before use.
#' @source Elevation data via the `microclimdata` package's `dem_download()` function, which uses
#'   the `elevatr` package to download global elevation data from AWS-hosted terrain tiles, matched
#'   to the extent and resolution of whatever template raster it is given; land/sea mask via Natural
#'   Earth (1:10m physical vectors, "land" category,
#'   <https://www.naturalearthdata.com>).
#' @seealso [landcover], [coastalexposure()]
#' @examples
#' \dontrun{
#' library(terra)
#' plot(rast(dtmc$dtm1))
#'
#' # Pass the whole list straight in as 'coarse' to extend a
#' # coastalexposure() search well beyond dtm100m's own extent -- its
#' # elements are wrapped (terra::PackedSpatRaster) and get unwrapped
#' # automatically, no need to call rast() on them first
#' landsea <- rast(dtm100m)
#' ce <- coastalexposure(landsea, wdir = 270, coarse = dtmc)
#' }
"dtmc"
