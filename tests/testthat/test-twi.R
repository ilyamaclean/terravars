library(terra)

test_that("twi returns SpatRaster of correct dimensions", {
  dtm <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000,
              ymin = 0, ymax = 1000, crs = "EPSG:27700")
  values(dtm) <- runif(100, 50, 150)
  result <- twi(dtm, method = "modified")
  expect_s4_class(result, "SpatRaster")
  expect_equal(terra::nrow(result), 10)
  expect_equal(terra::ncol(result), 10)
})

test_that("twi method argument is validated", {
  dtm <- rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500,
              ymin = 0, ymax = 500, crs = "EPSG:27700")
  values(dtm) <- 100
  expect_error(twi(dtm, method = "nonsense"))
})

test_that("all three TWI methods run without error", {
  dtm <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000,
              ymin = 0, ymax = 1000, crs = "EPSG:27700")
  values(dtm) <- runif(100, 50, 200)
  for (m in c("standard", "modified", "saga")) {
    expect_s4_class(twi(dtm, method = m), "SpatRaster")
  }
})
