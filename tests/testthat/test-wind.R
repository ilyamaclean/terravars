library(terra)

test_that("wind_shelter returns one layer per direction", {
  dtm <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000,
              ymin = 0, ymax = 2000, crs = "EPSG:27700")
  values(dtm) <- runif(400, 50, 200)
  dirs <- c(0, 90, 180, 270)
  ws   <- wind_shelter(dtm, wdir = dirs)
  expect_s4_class(ws, "SpatRaster")
  expect_equal(terra::nlyr(ws), length(dirs))
})

test_that("wind_shelter layer names match directions", {
  dtm <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000,
              ymin = 0, ymax = 1000, crs = "EPSG:27700")
  values(dtm) <- 100
  ws <- wind_shelter(dtm, wdir = c(45, 225))
  expect_equal(names(ws), c("wind_shelter_45", "wind_shelter_225"))
})

test_that("wind_shelter is close to 1 (unsheltered) on flat terrain", {
  dtm <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000,
              ymin = 0, ymax = 1000, crs = "EPSG:27700")
  values(dtm) <- 100
  ws <- wind_shelter(dtm, wdir = 270)
  expect_true(all(values(ws) >= 0 & values(ws) <= 1, na.rm = TRUE))
  expect_equal(mean(values(ws), na.rm = TRUE), 1, tolerance = 1e-6)
})
