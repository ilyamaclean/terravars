library(terra)

test_that("skyview on flat terrain returns 1 everywhere", {
  dtm <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000,
              ymin = 0, ymax = 2000, crs = "EPSG:27700")
  values(dtm) <- 100
  svf <- skyview(dtm, ndir = 8L)
  expect_s4_class(svf, "SpatRaster")
  v <- values(svf)
  expect_true(all(!is.na(v)))
  expect_equal(mean(v), 1, tolerance = 1e-6)
})

test_that("skyview values are in [0, 1]", {
  set.seed(1)
  dtm <- rast(nrows = 15, ncols = 15, xmin = 0, xmax = 1500,
              ymin = 0, ymax = 1500, crs = "EPSG:27700")
  values(dtm) <- runif(225, 0, 500)
  svf <- skyview(dtm, ndir = 8L)
  v <- values(svf)
  expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
})
