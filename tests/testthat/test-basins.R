library(terra)

# A flat 10x10 plain (elevation 100) with two pits: (row 2, col 2) at
# elevation 10, and (row 8, col 8) at elevation 5 (the global minimum,
# discovered/seeded first). Because every flat cell shares an identical
# elevation with its neighbours, the whole plain ties into whichever basin
# reaches it first (the one seeded at the lower pit, elevation 5), leaving
# the shallower pit (10) surrounded entirely by already-claimed cells --
# it can never grow and stays a 1-cell basin of its own. This gives a
# small, fully deterministic two-basin scenario to test against.
two_pit_dtm <- function() {
  m <- matrix(100, 10, 10)
  m[2, 2] <- 10
  m[8, 8] <- 5
  rast(m, extent = ext(0, 10, 0, 10), crs = "EPSG:27700")
}

test_that("basindelin() rejects non-SpatRaster input", {
  expect_error(basindelin(matrix(1, 5, 5)))
})

test_that("basindelin() finds two basins for two separated pits", {
  dtm <- two_pit_dtm()
  bsn <- basindelin(dtm)
  expect_s4_class(bsn, "SpatRaster")
  ids <- na.omit(as.vector(values(bsn)))
  expect_equal(length(unique(ids)), 2)
  # the shallower, surrounded pit ends up isolated as a 1-cell basin
  expect_equal(sort(as.vector(table(ids))), c(1, 99))
})

test_that("basindelin() assigns a single basin to a simple monotonic bowl", {
  m <- outer(1:10, 1:10, function(i, j) (i - 5)^2 + (j - 5)^2)
  dtm <- rast(m, extent = ext(0, 10, 0, 10), crs = "EPSG:27700")
  bsn <- basindelin(dtm)
  ids <- na.omit(as.vector(values(bsn)))
  expect_equal(length(unique(ids)), 1)
})

test_that("basindelin() leaves NA cells as NA", {
  m <- matrix(100, 10, 10)
  m[1:3, 1:3] <- NA
  dtm <- rast(m, extent = ext(0, 10, 0, 10), crs = "EPSG:27700")
  bsn <- basindelin(dtm)
  expect_true(all(is.na(as.matrix(bsn, wide = TRUE)[1:3, 1:3])))
})

test_that("basinmerge() rejects mismatched inputs", {
  dtm <- two_pit_dtm()
  bsn <- basindelin(dtm)
  expect_error(basinmerge(matrix(1, 5, 5), bsn))
  expect_error(basinmerge(dtm, matrix(1, 5, 5)))
  small <- rast(matrix(1, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  expect_error(basinmerge(dtm, small))
})

test_that("basinmerge() only merges once the barrier is within 'boundary'", {
  dtm <- two_pit_dtm()
  bsn <- basindelin(dtm)

  # Pour point between the two basins here sits at elevation 100 (the flat
  # plain separating them), vs. basin floors of 5 and 10 -- a barrier of
  # 95 m. A small boundary shouldn't bridge that.
  bm_small <- basinmerge(dtm, bsn, boundary = 1)
  expect_equal(length(unique(na.omit(as.vector(values(bm_small))))), 2)

  # A boundary comfortably larger than the barrier should merge them.
  bm_big <- basinmerge(dtm, bsn, boundary = 1000)
  expect_equal(length(unique(na.omit(as.vector(values(bm_big))))), 1)
})

test_that("fillsinks() removes interior pits", {
  set.seed(1)
  m <- matrix(runif(400, 0, 100), 20, 20)
  m[10, 10] <- -50   # deep artificial pit, well below its neighbours
  dtm <- rast(m, extent = ext(0, 20, 0, 20), crs = "EPSG:27700")
  filled <- fillsinks(dtm)
  expect_s4_class(filled, "SpatRaster")

  fm <- as.matrix(filled, wide = TRUE)
  nr <- nrow(fm); nc <- ncol(fm)
  # No interior cell should be strictly lower than every one of its 8
  # neighbours (i.e. no pits left to trap flow).
  has_pit <- FALSE
  for (i in 2:(nr - 1)) {
    for (j in 2:(nc - 1)) {
      nb <- fm[(i - 1):(i + 1), (j - 1):(j + 1)]
      nb[2, 2] <- NA
      if (all(fm[i, j] < nb, na.rm = TRUE)) has_pit <- TRUE
    }
  }
  expect_false(has_pit)
  # The filled pit should have been raised (never lowered further)
  expect_gt(fm[10, 10], m[10, 10])
})

test_that("fillsinks() flattens a water body to its own minimum elevation", {
  m  <- matrix(50, 10, 10)
  m[4:6, 4:6] <- c(10, 12, 11, 13, 9, 14, 10, 11, 12)
  wb <- matrix(NA, 10, 10)
  wb[4:6, 4:6] <- 1
  dtm <- rast(m,  extent = ext(0, 10, 0, 10), crs = "EPSG:27700")
  wbr <- rast(wb, extent = ext(0, 10, 0, 10), crs = "EPSG:27700")

  filled <- fillsinks(dtm, waterbody = wbr, water_value = 1)
  fm <- as.matrix(filled, wide = TRUE)
  expect_true(all(fm[4:6, 4:6] == min(m[4:6, 4:6])))
})

test_that("fillsinks() rejects a waterbody on a different grid", {
  dtm <- two_pit_dtm()
  small <- rast(matrix(1, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  expect_error(fillsinks(dtm, waterbody = small))
})
