library(terra)

# Minimal flat DTM for deterministic tests
flat_dtm <- function(nrow = 10, ncol = 10, val = 100) {
  r <- rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = 1000,
            ymin = 0, ymax = 1000, crs = "EPSG:27700")
  values(r) <- val
  r
}

test_that("solarindex returns SpatRaster with values in [0, 1]", {
  dtm <- flat_dtm()
  tme <- as.POSIXct("2023-06-21 12:00", tz = "UTC")
  si  <- solarindex(dtm, tme)
  expect_s4_class(si, "SpatRaster")
  expect_true(all(values(si) >= 0, na.rm = TRUE))
  expect_true(all(values(si) <= 1, na.rm = TRUE))
})

test_that("solarindex annual mode returns a single-layer SpatRaster in [0, 1]", {
  dtm <- flat_dtm()
  si_ann <- solarindex(dtm, "annual")
  expect_s4_class(si_ann, "SpatRaster")
  expect_equal(terra::nlyr(si_ann), 1)
  expect_true(all(values(si_ann) >= 0 & values(si_ann) <= 1, na.rm = TRUE))
})

# The internal solar_index_cpp() takes an explicit zenith/azimuth for a
# single solar position directly, which is the right level to pin down the
# exact geometric edge cases (overhead sun, sun below the horizon) without
# depending on real solar geometry for a specific date/location.
test_that("solar_index_cpp is 1 on flat terrain at normal incidence", {
  slope  <- matrix(0, 5, 5)
  aspect <- matrix(0, 5, 5)
  si <- terravars:::solar_index_cpp(slope, aspect, 0, 0, matrix(nrow = 0, ncol = 0))
  expect_equal(mean(si), 1, tolerance = 1e-6)
})

test_that("solar_index_cpp is 0 when sun is below the horizon", {
  slope  <- matrix(0, 5, 5)
  aspect <- matrix(0, 5, 5)
  si <- terravars:::solar_index_cpp(slope, aspect, 95, 180, matrix(nrow = 0, ncol = 0))
  expect_true(all(si == 0))
})

test_that("solar_index_cpp matches cos(zenith) on flat terrain", {
  slope  <- matrix(0, 5, 5)
  aspect <- matrix(0, 5, 5)
  si <- terravars:::solar_index_cpp(slope, aspect, 30, 180, matrix(nrow = 0, ncol = 0))
  expect_equal(mean(si), cos(30 * pi / 180), tolerance = 1e-6)
})
