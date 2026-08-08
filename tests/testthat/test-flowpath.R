library(terra)

test_that("flowpath() rejects non-SpatRaster input", {
  expect_error(flowpath(matrix(1, 5, 5), xy = c(0, 0)))
})

test_that("flowpath() rejects a malformed 'xy'", {
  dtm <- rast(matrix(100, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  expect_error(flowpath(dtm, xy = c(1, 2, 3)))
  expect_error(flowpath(dtm, xy = c(1, NA)))
})

test_that("flowpath() errors when 'xy' falls outside dtm's extent", {
  dtm <- rast(matrix(100, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  expect_error(flowpath(dtm, xy = c(-999, -999)))
})

test_that("flowpath() errors when 'xy' falls on an NA cell", {
  m <- matrix(100, 5, 5)
  m[1, 1] <- NA
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  expect_error(flowpath(dtm, xy = xyFromCell(dtm, cellFromRowCol(dtm, 1, 1))))
})

test_that("flowpath() traces the only downhill (diagonal-only) route under d8 but is stuck at a d4 pit", {
  # Centre cell (row 2, col 2) sits at 10; every orthogonal neighbour is a
  # higher wall (20), so under d4 it's a pit. Its NE diagonal neighbour
  # (row 1, col 3) sits at 5 -- the only way downhill, reachable only
  # under d8.
  m <- matrix(20, 3, 3)
  m[2, 2] <- 10
  m[1, 3] <- 5
  dtm <- rast(m, extent = ext(0, 3, 0, 3), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 2, 2))

  fp8 <- flowpath(dtm, xy = as.numeric(xy), method = "d8")
  fp4 <- flowpath(dtm, xy = as.numeric(xy), method = "d4")

  m8 <- as.matrix(fp8, wide = TRUE)
  m4 <- as.matrix(fp4, wide = TRUE)

  expect_equal(m8[2, 2], 2)
  expect_equal(m8[1, 3], 1)                 # reached diagonally under d8
  expect_equal(sum(m8 == 1), 1)

  expect_equal(m4[2, 2], 2)
  expect_equal(sum(m4 == 1), 0)             # d4: no orthogonal way down -- pit
})

test_that("flowpath() marks same-elevation neighbours of the start point as source cells and traces each", {
  # Row 3 is a flat shelf (elevation 47); every other row is strictly lower
  # the further south it is (elevation 50 - row). Starting at the centre of
  # row 3, its same-row orthogonal neighbours (row 3, cols 2 and 4) share
  # its own elevation exactly and should become sources too, each with its
  # own downslope trace.
  m <- matrix(0, 5, 5)
  for (i in 1:5) m[i, ] <- 50 - i
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 3, 3))

  fp <- flowpath(dtm, xy = as.numeric(xy), method = "d8")
  fm <- as.matrix(fp, wide = TRUE)

  # source row: start plus both same-elevation orthogonal neighbours
  expect_equal(fm[3, 2:4], c(2, 2, 2))
  # each traces straight downhill (south always beats a diagonal here,
  # since the drop is identical but the distance is shorter)
  expect_equal(fm[4, 2:4], c(1, 1, 1))
  expect_equal(fm[5, 2:4], c(1, 1, 1))
  # nothing north of the shelf, or outside the three traced columns, is touched
  expect_true(all(fm[1:2, ] == 0))
  expect_true(all(fm[4:5, c(1, 5)] == 0))
})

test_that("flowpath() with flat_source = FALSE uses only the start cell as source", {
  m <- matrix(0, 5, 5)
  for (i in 1:5) m[i, ] <- 50 - i   # same flat-shelf grid as above
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 3, 3))

  fp <- flowpath(dtm, xy = as.numeric(xy), flat_source = FALSE)
  fm <- as.matrix(fp, wide = TRUE)

  expect_equal(fm[3, 3], 2)
  expect_equal(fm[3, 2], 0)   # same-elevation neighbour, but not picked up this time
  expect_equal(fm[3, 4], 0)
  expect_equal(sum(fm == 2), 1)
})

test_that("flowpath() masks NA cells from dtm into the output", {
  m <- matrix(0, 5, 5)
  for (i in 1:5) m[i, ] <- 50 - i
  m[5, 5] <- NA
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 1, 1))

  fp <- flowpath(dtm, xy = as.numeric(xy))
  fm <- as.matrix(fp, wide = TRUE)
  expect_true(is.na(fm[5, 5]))
})

test_that("flowpath() names its output layer 'flowpath_<out>_<route>_<method>'", {
  dtm <- rast(matrix(100, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 3, 3))
  expect_equal(names(flowpath(dtm, xy = as.numeric(xy), method = "d8")),
               "flowpath_downstream_steepest_d8")
  expect_equal(names(flowpath(dtm, xy = as.numeric(xy), method = "d4")),
               "flowpath_downstream_steepest_d4")
  expect_equal(names(flowpath(dtm, xy = as.numeric(xy), route = "mfd", out = "upstream")),
               "flowpath_upstream_mfd_d8")
})

test_that("flowpath() delineates a real water body from 'waterbody'/'water_value' and traces from every cell of it", {
  # A 3x3 flat pond (elevation 47, same as row 3's shelf) sits in the
  # middle of the same flat-shelf grid used above; only 3 of its 9 cells
  # are immediately adjacent to the start point, so this must pick up more
  # source cells than flat_source's same-elevation-neighbour heuristic would.
  m <- matrix(0, 5, 5)
  for (i in 1:5) m[i, ] <- 50 - i
  m[2:4, 2:4] <- m[3, 3]
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")

  wb <- matrix(NA, 5, 5)
  wb[2:4, 2:4] <- 1
  wbr <- rast(wb, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")

  xy <- xyFromCell(dtm, cellFromRowCol(dtm, 3, 3))
  fp <- flowpath(dtm, xy = as.numeric(xy), waterbody = wbr, water_value = 1)
  fm <- as.matrix(fp, wide = TRUE)

  expect_true(all(fm[2:4, 2:4] == 2))   # the whole pond, not just xy's own neighbours
  expect_equal(sum(fm == 2), 9)
})

test_that("flowpath() falls back to a single-cell source with a warning when 'xy' isn't on a waterbody cell", {
  m <- matrix(0, 5, 5)
  for (i in 1:5) m[i, ] <- 50 - i
  m[2:4, 2:4] <- m[3, 3]
  dtm <- rast(m, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")

  wb <- matrix(NA, 5, 5)
  wb[2:4, 2:4] <- 1
  wbr <- rast(wb, extent = ext(0, 5, 0, 5), crs = "EPSG:27700")

  xy <- xyFromCell(dtm, cellFromRowCol(dtm, 1, 1))   # not on the pond
  expect_warning(fp <- flowpath(dtm, xy = as.numeric(xy), waterbody = wbr, water_value = 1))
  fm <- as.matrix(fp, wide = TRUE)
  expect_equal(sum(fm == 2), 1)
  expect_equal(fm[1, 1], 2)
})

test_that("flowpath() rejects a 'waterbody' on a different grid", {
  dtm <- rast(matrix(100, 5, 5), extent = ext(0, 5, 0, 5), crs = "EPSG:27700")
  small <- rast(matrix(1, 3, 3), extent = ext(0, 3, 0, 3), crs = "EPSG:27700")
  xy <- xyFromCell(dtm, cellFromRowCol(dtm, 3, 3))
  expect_error(flowpath(dtm, xy = as.numeric(xy), waterbody = small))
})

test_that("flowpath(route = 'mfd') splits water across tied downhill neighbours in proportion to their number, not always toward one", {
  # Centre cell (10) has two orthogonal downhill neighbours tied at
  # elevation 0 (E and S) -- both the same drop and the same distance, so
  # an even split is the only slope-consistent outcome regardless of the
  # weighting exponent.
  m <- matrix(20, 3, 3)
  m[2, 2] <- 10
  m[2, 3] <- 0   # E
  m[3, 2] <- 0   # S
  dtm <- rast(m, extent = ext(0, 3, 0, 3), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 2, 2))

  fp <- flowpath(dtm, xy = as.numeric(xy), route = "mfd")
  fm <- as.matrix(fp, wide = TRUE)

  expect_equal(fm[2, 2], 1)
  expect_equal(fm[2, 3], 0.5)
  expect_equal(fm[3, 2], 0.5)
  expect_true(all(fm >= 0 & fm <= 1))
})

test_that("flowpath(route = 'mfd') diffusion chains multiply out correctly (0.5 splitting again into two 0.25s)", {
  m <- matrix(1000, 9, 9)
  m[5, 5] <- 100   # source
  m[5, 6] <- 50    # E of source  -- ties with S
  m[6, 5] <- 50    # S of source  -- ties with E -> each gets 0.5
  m[7, 4] <- 0     # SW of (6,5)  -- ties with SE
  m[7, 6] <- 0     # SE of (6,5)  -- ties with SW -> each gets 0.25
  dtm <- rast(m, extent = ext(0, 9, 0, 9), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 5, 5))

  fp <- flowpath(dtm, xy = as.numeric(xy), route = "mfd")
  fm <- as.matrix(fp, wide = TRUE)

  expect_equal(fm[5, 5], 1)
  expect_equal(fm[5, 6], 0.5)
  expect_equal(fm[6, 5], 0.5)
  expect_equal(fm[7, 4], 0.25)
  expect_equal(fm[7, 6], 0.25)
})

test_that("flowpath(out = 'upstream', route = 'steepest') walks uphill, diagonal-only under d8, stuck at a d4 ridge", {
  # Mirror of the downstream d4-vs-d8 test: centre (row 2, col 2) sits at
  # 10; every orthogonal neighbour is a lower wall (5), so under d4 there's
  # no way further uphill. Its NE diagonal neighbour (row 1, col 3) sits at
  # 20 -- the only way up, reachable only under d8.
  m <- matrix(5, 3, 3)
  m[2, 2] <- 10
  m[1, 3] <- 20
  dtm <- rast(m, extent = ext(0, 3, 0, 3), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 2, 2))

  fp8 <- flowpath(dtm, xy = as.numeric(xy), method = "d8", out = "upstream")
  fp4 <- flowpath(dtm, xy = as.numeric(xy), method = "d4", out = "upstream")

  m8 <- as.matrix(fp8, wide = TRUE)
  m4 <- as.matrix(fp4, wide = TRUE)

  expect_equal(m8[2, 2], 2)
  expect_equal(m8[1, 3], 1)
  expect_equal(sum(m8 == 1), 1)

  expect_equal(m4[2, 2], 2)
  expect_equal(sum(m4 == 1), 0)
})

test_that("flowpath(out = 'upstream', route = 'mfd') grades each upslope cell by how much of its own water reaches the point", {
  # Same diffusion grid as the downstream chain test above, but now traced
  # upstream from one of the two leaves (row 7, col 4): only half of
  # (6, 5)'s water reaches that specific leaf (it splits the other half
  # toward (7, 6) instead), and only a quarter of the original source
  # (5, 5)'s water does (it first has to reach (6, 5) at all, itself only a
  # 50/50 split against the dead-end branch at (5, 6), which contributes
  # nothing to this particular leaf).
  m <- matrix(1000, 9, 9)
  m[5, 5] <- 100
  m[5, 6] <- 50
  m[6, 5] <- 50
  m[7, 4] <- 0
  m[7, 6] <- 0
  dtm <- rast(m, extent = ext(0, 9, 0, 9), crs = "EPSG:27700")
  xy  <- xyFromCell(dtm, cellFromRowCol(dtm, 7, 4))

  fp <- flowpath(dtm, xy = as.numeric(xy), route = "mfd", out = "upstream")
  fm <- as.matrix(fp, wide = TRUE)

  expect_equal(fm[7, 4], 1)     # the point itself
  expect_equal(fm[6, 5], 0.5)
  expect_equal(fm[5, 5], 0.25)
  expect_equal(fm[5, 6], 0)     # dead-end branch -- none of its water reaches this leaf
  expect_equal(fm[7, 6], 0)     # the other leaf -- its own water goes nowhere
})
