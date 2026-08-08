#include "utils.h"
#include <string>
#include <queue>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ============================================================
// wetness.cpp
//
// Soil moisture / hydrological functions:
//   - twi_cpp()         topographic wetness index (standard / modified / saga)
//   - basinCpp()        hydrological basin delineation
//   - renumberbasin()   renumber basins sequentially
// ============================================================

// ------------------------------------------------------------
// twi_cpp
//
// Topographic Wetness Index:  TWI = ln(a / tan(beta))
//   a    = specific catchment area (flow accumulation * res)
//   beta = local slope (radians)
//
// method:
//   "standard"  - classical Beven & Kirkby (1979)
//   "modified"  - uses a minimum slope floor to avoid ln(Inf)
//                 (Gruber & Peckham 2009)
//   "saga"      - SAGA wetness index: uses modified CA and slope
//                 (Boehner & Selige 2006); requires t_slope input
//
// Parameters:
//   flow_acc  - flow accumulation (cells) -- expected SELF-INCLUSIVE,
//               i.e. every cell's own weight already counted (matching
//               flowacc()'s convention: a headwater cell reads 1, never
//               0), so `a` here is simply flow_acc * res directly, no
//               +1 adjustment. (Older ArcGIS-style flow accumulation
//               excludes the cell itself, which is why an earlier
//               version of this function added 1 by hand -- now that
//               flowacc() is the package's own self-inclusive source,
//               that adjustment has moved into the convention itself
//               rather than being patched on here.)
//   slope_deg - slope in degrees
//   res       - cell resolution (metres)
//   method    - one of "standard", "modified", "saga"
//   min_slope - minimum slope (degrees) used in "modified"
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix twi_cpp(Rcpp::NumericMatrix flow_acc,
                             Rcpp::NumericMatrix slope_deg,
                             double res,
                             std::string method,
                             double min_slope) {

    int nrows = flow_acc.nrow();
    int ncols = flow_acc.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);

    double min_slope_r = d2r(min_slope);

    // Hoist the method string comparison out of the parallel loop
    bool is_modified = (method == "modified");
    bool is_saga     = (method == "saga");
    bool floor_slope = is_modified || is_saga;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            double sca = flow_acc(r, c) * res;
            double sl  = d2r(slope_deg(r, c));

            if (floor_slope && sl < min_slope_r) sl = min_slope_r;

            // "modified"/"saga" floor slope to min_slope_r above, so tan_sl
            // is normally already > 0 there -- this only guards the edge
            // case of a user-supplied min_slope of exactly 0. "standard"
            // has no floor at all: a genuinely flat cell (tan_sl == 0) is
            // meant to divide out to Inf (see twi()'s documented "standard"
            // behaviour in R/wetness.R), so it must NOT be clamped here.
            double tan_sl = std::tan(sl);
            if (floor_slope && tan_sl <= 0.0) tan_sl = 1e-6;

            if (is_saga && sca > res * 1e6) sca = res * 1e6;

            out(r, c) = std::log(sca / tan_sl);
        }
    }
    return out;
}

// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for delineating basins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// Basin delineation, one basin at a time, grown outward and upward from a
// local pit. Cells at exactly the same elevation as an already-claimed
// cell always merge into the same basin -- there's no well-defined
// "downhill direction" across genuinely flat ground (a plateau, a lake
// surface), so ties are never restricted by `method`. How a STRICTLY
// higher neighbour is allowed to join, though, depends on `method`:
//   - "any": it joins as soon as anything adjoining it is claimed and it's
//     equal-or-higher, regardless of whether that's the steepest way down
//     from the higher cell's own point of view. A basin can grow through
//     any non-decreasing connection.
//   - "steepest": it only joins through its own single steepest downhill
//     neighbour -- the one direction its water would actually flow, found
//     via the `sdrow`/`sdcol` precomputation below. A higher cell can no
//     longer be claimed by a basin that merely reaches it first; only the
//     one basin that owns its actual steepest-descent target ever can.
// A genuine local peak/ridge still blocks further growth either way, since
// nothing downhill of it can be claimed by definition.
//
// Cells are visited in strictly ascending elevation order, one basin at a
// time: a new basin is seeded at whichever unclaimed cell is currently
// lowest anywhere in the grid, and that basin is grown to complete
// exhaustion -- repeatedly claiming whichever unclaimed cell is lowest
// within it -- before the next basin is even considered. Under "any", this
// means whichever pit is discovered first (i.e. is globally lowest) gets
// first claim on any contested ground reachable from it, ahead of any
// other, not-yet-considered pit -- ties aside, "steepest" has no such
// contested ground to resolve, since a strictly-higher cell's fate is
// fixed by its own local gradient rather than by discovery order.
//
// Implemented via two min-heaps rather than the O(rows*cols) linear scans
// this originally used (giving O(N log N) instead of O(N^2) overall),
// without changing the result:
//   - `seedheap`: every real interior cell (see the exclusions below),
//     pushed once, up front -- finds "the next lowest not-yet-claimed cell
//     anywhere". Popped lazily: an entry for a cell that's since been
//     claimed by an earlier basin's growth is simply discarded and the
//     next entry tried. Each cell is pushed once and popped at most twice
//     (once here, once more below if it becomes part of a basin), so the
//     total work stays O(N log N).
//   - a fresh, LOCAL heap per basin, seeded with just that basin's seed
//     cell -- finds "the lowest still-pending cell within the current
//     basin". Pop the lowest pending cell, claim any still-unclaimed >=
//     neighbour, and push each newly-claimed neighbour onto the SAME local
//     heap so the basin keeps growing outward in ascending order until the
//     local heap runs dry, i.e. until the basin is completely finished.
// Only once a basin's local heap is empty does the outer loop return to
// `seedheap` for the next basin's seed, preserving the "finish one basin
// entirely before starting the next" priority described above.
//
// Two kinds of cell are excluded from `seedheap` entirely:
//   - the 1-cell padding border (row/col 0 and row rows-1/col cols-1) --
//     always 9999, stripped off by the R wrapper before the result comes
//     back, so there's nothing gained by visiting it;
//   - interior NA/no-data cells, also encoded as 9999 by the R wrapper
//     (basindelin() does `dm[is.na(dm)] <- 9999` before calling down here)
//     -- e.g. sea cells on a coastal DTM. Excluding these leaves them
//     NA/ignored in the result (nothing can claim a 9999 cell, per the
//     `!= 9999` guard below) rather than each becoming its own singleton
//     basin.
struct BasinHeapCell {
    double elev;
    int row;
    int col;
};
// std::priority_queue is a max-heap by default; reversing the comparison
// keeps the lowest elevation on top, i.e. a min-heap. Elevation ties are
// broken by row then column, ascending -- i.e. whichever tied cell would be
// encountered first scanning row by row, left to right. This matters
// because which pit is treated as "the seed" determines which basin gets
// first claim on any contested ground between two competing pits (see the
// discussion above); a tied pair of seed candidates needs a deterministic,
// well-defined tie-break rather than whatever order the heap's internal
// structure happens to produce.
struct BasinHeapCompare {
    bool operator()(const BasinHeapCell& a, const BasinHeapCell& b) const {
        if (a.elev != b.elev) return a.elev > b.elev;
        if (a.row  != b.row)  return a.row  > b.row;
        return a.col > b.col;
    }
};
// Function that does the basin delineation.
// `method`: "any" (default behaviour) or "steepest" -- see the discussion
// above.
// [[Rcpp::export]]
IntegerMatrix basinCpp(NumericMatrix& dm2, IntegerMatrix& bsn, std::string method) {
    int rows = dm2.nrow();
    int cols = dm2.ncol();
    int bn = 1;
    bool steepest = (method == "steepest");

    // Precompute each interior cell's single steepest downhill neighbour,
    // only needed for "steepest". sdrow/sdcol store its row/col, or -1 if
    // the cell has no strictly-lower real neighbour at all (a genuine pit,
    // or bordered only by 9999 cells) -- such a cell can still be claimed
    // via a tie (see the claim rule below), just never via this table.
    // Diagonal neighbours are ~1.41x farther away than orthogonal ones, so
    // their raw elevation drop is divided by sqrt(2) before comparing --
    // otherwise a diagonal neighbour with a merely-larger drop could beat
    // a closer orthogonal one that's actually steeper. Ties for steepest
    // (rare: two neighbours with the exact same distance-normalised drop)
    // resolve to whichever is hit first scanning the 8 neighbours in
    // reading order (top-left to bottom-right), since only a strictly
    // greater slope replaces the running best below.
    IntegerMatrix sdrow(rows, cols), sdcol(rows, cols);
    if (steepest) {
        std::fill(sdrow.begin(), sdrow.end(), -1);
        std::fill(sdcol.begin(), sdcol.end(), -1);
        for (int i = 1; i < rows - 1; ++i) {
            for (int j = 1; j < cols - 1; ++j) {
                if (dm2(i, j) == 9999) continue;
                double best_slope = -1.0; // any real downhill slope is > 0
                int best_ni = -1, best_nj = -1;
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        if (di == 0 && dj == 0) continue;
                        int ni = i + di, nj = j + dj;
                        if (dm2(ni, nj) == 9999) continue;
                        if (dm2(ni, nj) >= dm2(i, j)) continue; // not strictly lower
                        double dist = (di != 0 && dj != 0) ? std::sqrt(2.0) : 1.0;
                        double slope = (dm2(i, j) - dm2(ni, nj)) / dist;
                        if (slope > best_slope) {
                            best_slope = slope;
                            best_ni = ni;
                            best_nj = nj;
                        }
                    }
                }
                sdrow(i, j) = best_ni;
                sdcol(i, j) = best_nj;
            }
        }
    }

    std::priority_queue<BasinHeapCell, std::vector<BasinHeapCell>, BasinHeapCompare> seedheap;

    // Push every real interior cell once -- stands in for whichmin()'s
    // global "next lowest unclaimed cell" search. See exclusions above.
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (dm2(i, j) != 9999) {
                seedheap.push(BasinHeapCell{dm2(i, j), i, j});
            }
        }
    }

    while (!seedheap.empty()) {
        BasinHeapCell seed = seedheap.top();
        seedheap.pop();

        // Lazy deletion: this cell may already have been claimed by an
        // earlier basin's growth (below) since it was pushed. If so, it's
        // not a new seed -- skip straight to the next entry.
        if (!IntegerMatrix::is_na(bsn(seed.row, seed.col))) continue;

        // New basin, seeded here. Grow it to complete exhaustion via its
        // own local heap before this outer loop ever looks at seedheap
        // again -- see the "one basin at a time" note above.
        bsn(seed.row, seed.col) = bn;
        std::priority_queue<BasinHeapCell, std::vector<BasinHeapCell>, BasinHeapCompare> localheap;
        localheap.push(seed);

        while (!localheap.empty()) {
            BasinHeapCell c = localheap.top();
            localheap.pop();
            int i = c.row;
            int j = c.col;
            double d = dm2(i, j);
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    if (di == 0 && dj == 0) continue;
                    int ni = i + di;
                    int nj = j + dj;
                    if (dm2(ni, nj) == 9999 || !IntegerMatrix::is_na(bsn(ni, nj))) continue;

                    // Ties always merge regardless of method -- see the
                    // discussion above. A strictly higher neighbour merges
                    // under "any" unconditionally, but under "steepest"
                    // only if THIS cell is specifically its own
                    // steepest-descent target (sdrow/sdcol), i.e. only the
                    // one basin its water would actually flow into can
                    // claim it.
                    bool claim = (dm2(ni, nj) == d) ||
                        (dm2(ni, nj) > d &&
                         (!steepest || (sdrow(ni, nj) == i && sdcol(ni, nj) == j)));

                    if (claim) {
                        bsn(ni, nj) = bn;
                        localheap.push(BasinHeapCell{dm2(ni, nj), ni, nj});
                    }
                }
            }
        }
        bn++;
    }
    return bsn;
}

// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for merging basins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// Pour-point-based basin merging (see basinmerge() in R/wetness.R for the
// full explanation). Basin A and basin B are merged when the lowest point
// at which you could cross from one into the other (their pour point, or
// col) sits within `boundary` metres of the lower of the two basins' own
// floor elevations (each basin's minimum dtm value) -- i.e. the shallower
// basin only has to fill by a little before it spills into its neighbour.
// This deliberately does NOT look at how smoothly the terrain happens to
// slope at any one crossing point -- only the height of the lowest
// crossing relative to basin floor matters, since a high ridge that
// happens to be locally flat somewhere along its length is not actually a
// shallow barrier.
//
// Same padded-matrix convention as basinCpp(): dm2/bm2 already carry a
// 1-cell border (elevation 9999, basin id NA) added by the R wrapper, and
// the returned matrix keeps that border for the wrapper to strip.
//
// Basin ids need not be contiguous (e.g. bsn from a cropped raster), so
// they're mapped to a compact 0-based index internally. The output keeps
// each merged group's ids as this internal, still-gappy index (+1) --
// NOT yet renumbered sequentially -- since the R wrapper does that as an
// explicit final step via .renumber_seq(), the same as it would for a
// pure-R merge.
// [[Rcpp::export]]
IntegerMatrix basinmerge_cpp(NumericMatrix& dm2, IntegerMatrix& bm2, double boundary) {
    int rows = dm2.nrow();
    int cols = dm2.ncol();

    // ---- distinct basin ids present, mapped to a compact 0-based index ----
    std::unordered_map<int,int> id_index;
    std::vector<int> ids;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (IntegerMatrix::is_na(bm2(i, j))) continue;
            int v = bm2(i, j);
            if (id_index.find(v) == id_index.end()) {
                id_index[v] = (int)ids.size();
                ids.push_back(v);
            }
        }
    }
    int n_ids = (int)ids.size();
    IntegerMatrix out(rows, cols);
    if (n_ids == 0) {
        std::fill(out.begin(), out.end(), NA_INTEGER);
        return out;
    }

    // ---- basin floor elevation: each basin's own minimum dtm value ----
    std::vector<double> floor_by_id(n_ids, 9999.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (IntegerMatrix::is_na(bm2(i, j))) continue;
            int k = id_index[bm2(i, j)];
            double d = dm2(i, j);
            if (d < floor_by_id[k]) floor_by_id[k] = d;
        }
    }

    // ---- scan every boundary crossing between two different basins,
    // tracking the minimum "toll" (higher of the two adjoining cells'
    // elevations, i.e. the height you'd have to reach to cross there) per
    // unordered basin-id pair -- that minimum is the pair's pour point.
    // Only the "forward" half of the 8 neighbours (E, SE, S, SW) is
    // scanned from each cell, so every boundary crossing is visited
    // exactly once overall (its mirror image is covered from the other
    // side when that other cell has its own turn), rather than twice. ----
    std::unordered_map<long long, double> pour;
    static const int dOff[4][2] = {{0,1},{1,-1},{1,0},{1,1}};
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (IntegerMatrix::is_na(bm2(i, j))) continue;
            int v = bm2(i, j);
            int kself = id_index[v];
            double dself = dm2(i, j);
            for (auto& off : dOff) {
                int ni = i + off[0], nj = j + off[1];
                if (IntegerMatrix::is_na(bm2(ni, nj))) continue;
                int vn = bm2(ni, nj);
                if (vn == v) continue;
                int knbr = id_index[vn];
                double dnbr = dm2(ni, nj);
                double toll = std::max(dself, dnbr);
                int lo = std::min(kself, knbr), hi = std::max(kself, knbr);
                long long key = (long long)lo * n_ids + hi;
                auto it = pour.find(key);
                if (it == pour.end() || toll < it->second) pour[key] = toll;
            }
        }
    }

    // ---- union-find over basin indices, merging wherever the pour point
    // is shallow relative to the lower of the two basins' floors ----
    std::vector<int> parent(n_ids);
    for (int k = 0; k < n_ids; ++k) parent[k] = k;
    std::function<int(int)> find_root = [&](int k) {
        int root = k;
        while (parent[root] != root) root = parent[root];
        while (parent[k] != root) { int nxt = parent[k]; parent[k] = root; k = nxt; }
        return root;
    };
    for (auto& kv : pour) {
        int lo = (int)(kv.first / n_ids);
        int hi = (int)(kv.first % n_ids);
        double barrier = kv.second - std::min(floor_by_id[lo], floor_by_id[hi]);
        if (barrier < boundary) {
            int ra = find_root(lo), rb = find_root(hi);
            if (ra != rb) {
                int a = std::min(ra, rb), b = std::max(ra, rb);
                parent[b] = a;
            }
        }
    }

    // ---- relabel every cell by its merged group's root. Root values are
    // just whichever id-index happened to end up on top of each
    // union-find tree, so they're a sparse/gappy subset of 1:n_ids at this
    // point -- the R wrapper renumbers sequentially afterwards. ----
    std::vector<int> roots(n_ids);
    for (int k = 0; k < n_ids; ++k) roots[k] = find_root(k);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (IntegerMatrix::is_na(bm2(i, j))) { out(i, j) = NA_INTEGER; continue; }
            int k = id_index[bm2(i, j)];
            out(i, j) = roots[k] + 1;
        }
    }
    return out;
}

// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for filling sinks/depressions ~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// Priority-flood depression filling (Barnes et al. 2014 / Wang & Liu 2006):
// every interior pit is raised just high enough that a continuous
// non-increasing path exists from it down to an outlet -- either the
// raster edge, an NA/no-data cell (this package's "sea" convention, see
// basindelin()/basinmerge() above), or now optionally a real water body
// (see wseed2 below). A single min-heap grows inward from every outlet
// simultaneously: seed the heap with every cell that already touches an
// outlet (at its own, unmodified elevation -- it's already valid to
// drain there), then repeatedly pop the lowest resolved cell and resolve
// each of its unresolved neighbours to
// max(neighbour's own elevation, the popped cell's resolved elevation) --
// i.e. a neighbour only gets raised if it would otherwise sit below the
// cell it has to drain through, and only ever up to the true pour point
// of whatever depression it's part of, never further. This is simpler
// than basinCpp()'s two-tier heap scheme above (global seed heap + a
// fresh local heap grown to full completion per basin) -- that
// complexity exists only to reproduce basindelin()'s specific, historical
// basin-discovery-order semantics; there's no equivalent legacy
// behaviour to match here, so one heap growing outward from every
// outlet at once is both simpler and the standard textbook algorithm.
//
// Reuses BasinHeapCell/BasinHeapCompare (defined above, ahead of
// basinCpp()) -- exactly the min-heap-by-elevation-with-deterministic-
// tie-break this needs too, no reason to duplicate it.
//
// dm2 carries a 1-cell border (elevation 9999, this package's NA/no-data
// sentinel) added by the R wrapper -- this doubles as "every real cell on
// the raster's own edge already touches an outlet" for free, since the
// padding ring surrounds the whole raster. Interior 9999 cells (genuine
// no-data, e.g. a sea/lake pocket) work exactly the same way: any real
// cell touching one is already free to drain there, at whatever its own
// elevation is. The returned matrix keeps the border for the wrapper to
// strip; 9999 cells are returned unchanged (the wrapper converts them
// back to real NA).
//
// wseed2 is the same shape as dm2 (1-cell border included) and marks real
// water bodies: a real (non-NA) value at (i, j) means "this cell belongs
// to a real pond/lake, already resolved at this flattened elevation" --
// the R wrapper computes this once per connected water body (the minimum
// dtm elevation found anywhere in that water body) and paints it across
// every cell of it, so every water-body cell reads the same constant.
// Genuine NA (R's NA_real_, checked with R_IsNA -- ordinary NaN from
// arithmetic would not count, but the wrapper only ever writes NA_real_
// or a real elevation here) means "not a water body, resolve normally".
// Every wseed2 cell is seeded into the heap up front, before the usual
// outlet scan below, at its given elevation -- this covers interior water
// bodies that never touch the raster edge or an NA/sea cell at all, not
// just ones that do. Seeding the *entire* water body rather than just its
// rim is deliberate: it lets the ordinary priority-flood growth below
// discover each water body's true pour point for free, exactly the way it
// already discovers pour points for ordinary pits -- whichever
// neighbouring land cell has the lowest-cost path in is naturally the
// first one the growth phase reaches, with no separate pour-point search
// needed. A cell can't be both a 9999 (no-data) cell and a water-body
// seed in any input this wrapper produces, but the 9999 check still comes
// first defensively -- a cell with no real elevation has nothing to seed.
// [[Rcpp::export]]
NumericMatrix fillsinks_cpp(NumericMatrix& dm2, NumericMatrix& wseed2) {
    int rows = dm2.nrow();
    int cols = dm2.ncol();
    NumericMatrix filled(rows, cols);
    std::fill(filled.begin(), filled.end(), 9999.0);
    std::vector<char> resolved(rows * cols, 0);
    auto idx = [cols](int i, int j) { return i * cols + j; };

    std::priority_queue<BasinHeapCell, std::vector<BasinHeapCell>, BasinHeapCompare> heap;

    // Seed every real interior cell belonging to a water body directly at
    // its already-known, flattened elevation -- resolved outright, same
    // as an outlet, whether or not it happens to also touch the raster
    // edge or an NA/sea cell.
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (dm2(i, j) == 9999.0) continue;
            if (R_IsNA(wseed2(i, j))) continue;
            double we = wseed2(i, j);
            filled(i, j) = we;
            resolved[idx(i, j)] = 1;
            heap.push(BasinHeapCell{we, i, j});
        }
    }

    // Seed every remaining real interior cell that already touches an
    // outlet (the padding ring, or a genuine interior NA/no-data cell) at
    // its own, unmodified elevation.
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (dm2(i, j) == 9999.0) continue;
            if (resolved[idx(i, j)]) continue;
            bool touches_outlet = false;
            for (int di = -1; di <= 1 && !touches_outlet; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    if (di == 0 && dj == 0) continue;
                    if (dm2(i + di, j + dj) == 9999.0) { touches_outlet = true; break; }
                }
            }
            if (touches_outlet) {
                filled(i, j) = dm2(i, j);
                resolved[idx(i, j)] = 1;
                heap.push(BasinHeapCell{dm2(i, j), i, j});
            }
        }
    }

    // Standard priority-flood growth -- unchanged; treats every seed
    // (ordinary outlet or water-body cell) identically once in the heap.
    while (!heap.empty()) {
        BasinHeapCell c = heap.top();
        heap.pop();
        int i = c.row, j = c.col;
        double d = filled(i, j);
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                if (di == 0 && dj == 0) continue;
                int ni = i + di, nj = j + dj;
                if (dm2(ni, nj) == 9999.0) continue;
                if (resolved[idx(ni, nj)]) continue;
                double rv = std::max(dm2(ni, nj), d);
                filled(ni, nj) = rv;
                resolved[idx(ni, nj)] = 1;
                heap.push(BasinHeapCell{rv, ni, nj});
            }
        }
    }

    // Defensive: any real cell somehow never reached (shouldn't happen --
    // every real cell's connected component of non-outlet cells must
    // include at least one cell touching the padding ring or an interior
    // NA cell) keeps its original elevation rather than the 9999
    // unresolved marker.
    for (int i = 1; i < rows - 1; ++i)
        for (int j = 1; j < cols - 1; ++j)
            if (dm2(i, j) != 9999.0 && !resolved[idx(i, j)]) filled(i, j) = dm2(i, j);

    return filled;
}

// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for flow accumulation ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// flowacc_cpp(): D8 / MFD (Freeman 1991, convergence exponent fixed at
// 1.1 -- deliberately not exposed as a user argument, since nobody
// calibrating it by feel is safer than a knob nobody understands) /
// D-infinity (Tarboton 1997) flow accumulation, with an optional per-cell
// weight raster and optional basin isolation. This is a direct,
// numerically-verified port of an R/Python reference prototype built and
// tested cell-by-cell before writing any C++ (linear channels, symmetric
// and asymmetric splits, flat/plateau routing, basin isolation, and basin
// isolation combined with flats all independently checked).
//
// dm2/wt2/bm2 all carry the usual 1-cell padding border (elevation 9999 =
// this package's NA/no-data sentinel; weight 0; basin id NA_INTEGER = "no
// basin"). When the caller doesn't supply a basin raster, the R wrapper
// fills bm2 with a single dummy basin id everywhere a real cell exists,
// so no separate "basin restriction on/off" flag is needed here -- the
// same-basin test is simply always true in that case, matching
// unrestricted routing.
//
// method: 0 = d8, 1 = mfd, 2 = dinf.
//
// Two-phase algorithm:
//  1. Distance-to-exit BFS. A cell "has an exit" if it touches true NA
//     (always a real outlet, basin-independent -- leaving the modelled
//     domain entirely) or a same-basin neighbour at strictly lower
//     elevation. A neighbour in a *different* basin is excluded from
//     consideration entirely -- it must NOT itself count as an exit, or
//     basin isolation would leak at every basin boundary (a cell sitting
//     right on the boundary would wrongly get dist=0 without ever having
//     found a genuine lower or true-outlet path). Cells with an exit get
//     distance 0; distance then propagates outward across chains of
//     identical-elevation, same-basin neighbours (a flat/plateau, e.g. a
//     fillsinks() fill) toward whichever dist-0 cell is nearest, giving
//     every flat cell a well-defined route down to a real exit.
//  2. Process every real cell in descending-elevation order, and within a
//     tied elevation, descending distance-to-exit order (farthest-from-
//     exit first). This ordering is essential for flats: a cell must
//     forward its accumulated total only after every cell that could
//     feed it has already been resolved, otherwise a cell closer to the
//     exit gets finalised (and forwards its total onward) before a
//     farther cell's contribution has even arrived, silently losing it.
//     Each cell's running total (own weight plus whatever it has
//     received so far) is split among its downhill/flat receivers per
//     `method`; a cell with no valid receiver at all (an unfilled pit, or
//     a basin with no reachable exit under isolation) simply keeps its
//     total and propagates no further.
static const int FA_DI[8] = { 0, -1, -1, -1,  0,  1, 1, 1 };
static const int FA_DJ[8] = { 1,  1,  0, -1, -1, -1, 0, 1 };
static const double FA_MFD_EXPONENT = 1.1;

// [[Rcpp::export]]
NumericMatrix flowacc_cpp(NumericMatrix& dm2, NumericMatrix& wt2, IntegerMatrix& bm2, int method) {
    int rows = dm2.nrow();
    int cols = dm2.ncol();
    auto idx = [cols](int i, int j) { return i * cols + j; };
    auto valid = [&](int i, int j) { return dm2(i, j) != 9999.0; };
    auto same_basin = [&](int i, int j, int ni, int nj) {
        if (IntegerMatrix::is_na(bm2(i, j)) || IntegerMatrix::is_na(bm2(ni, nj))) return false;
        return bm2(i, j) == bm2(ni, nj);
    };

    // ---- phase 1: distance-to-exit BFS ----
    std::vector<int> dist(rows * cols, -1); // -1 = unresolved (no exit reachable)
    std::queue<std::pair<int, int>> q;
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            if (!valid(i, j)) continue;
            bool has_exit = false;
            for (int k = 0; k < 8; ++k) {
                int ni = i + FA_DI[k], nj = j + FA_DJ[k];
                if (!valid(ni, nj)) { has_exit = true; continue; }
                if (!same_basin(i, j, ni, nj)) continue; // excluded, not an exit
                if (dm2(ni, nj) < dm2(i, j)) has_exit = true;
            }
            if (has_exit) {
                dist[idx(i, j)] = 0;
                q.push({ i, j });
            }
        }
    }
    while (!q.empty()) {
        auto cur = q.front(); q.pop();
        int i = cur.first, j = cur.second;
        for (int k = 0; k < 8; ++k) {
            int ni = i + FA_DI[k], nj = j + FA_DJ[k];
            if (ni < 1 || ni >= rows - 1 || nj < 1 || nj >= cols - 1) continue;
            if (!valid(ni, nj)) continue;
            if (!same_basin(ni, nj, i, j)) continue;
            if (dm2(ni, nj) != dm2(i, j)) continue;
            if (dist[idx(ni, nj)] == -1 || dist[idx(ni, nj)] > dist[idx(i, j)] + 1) {
                dist[idx(ni, nj)] = dist[idx(i, j)] + 1;
                q.push({ ni, nj });
            }
        }
    }

    // ---- phase 2: processing order (elevation desc, then dist desc) ----
    std::vector<std::pair<int, int>> cells;
    cells.reserve((size_t)rows * cols);
    for (int i = 1; i < rows - 1; ++i)
        for (int j = 1; j < cols - 1; ++j)
            if (valid(i, j)) cells.push_back({ i, j });

    std::sort(cells.begin(), cells.end(), [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        double ea = dm2(a.first, a.second), eb = dm2(b.first, b.second);
        if (ea != eb) return ea > eb;
        int da = dist[idx(a.first, a.second)], db = dist[idx(b.first, b.second)];
        if (da == -1 || db == -1) return db == -1 && da != -1; // unresolved cells last, order among themselves irrelevant
        return da > db;
    });

    NumericMatrix acc(rows, cols);
    std::fill(acc.begin(), acc.end(), 0.0);
    for (auto& c : cells) acc(c.first, c.second) += wt2(c.first, c.second);

    for (auto& c : cells) {
        int i = c.first, j = c.second;
        double total = acc(i, j);
        if (total <= 0) continue;
        double d0 = dm2(i, j);

        // gather strictly-downhill same-basin receivers
        std::vector<int> dh_k;
        std::vector<double> dh_slope;
        for (int k = 0; k < 8; ++k) {
            int ni = i + FA_DI[k], nj = j + FA_DJ[k];
            if (!valid(ni, nj) || !same_basin(i, j, ni, nj)) continue;
            if (dm2(ni, nj) < d0) {
                dh_k.push_back(k);
                dh_slope.push_back(d0 - dm2(ni, nj));
            }
        }

        if (!dh_k.empty()) {
            if (method == 0) { // d8: steepest wins, first-encountered breaks ties
                int best = 0;
                for (size_t t = 1; t < dh_k.size(); ++t)
                    if (dh_slope[t] > dh_slope[best]) best = (int)t;
                int k = dh_k[best];
                acc(i + FA_DI[k], j + FA_DJ[k]) += total;
            } else if (method == 1) { // mfd: slope^exponent-weighted split
                double s = 0.0;
                std::vector<double> w(dh_k.size());
                for (size_t t = 0; t < dh_k.size(); ++t) {
                    w[t] = std::pow(std::max(dh_slope[t], 1e-9), FA_MFD_EXPONENT);
                    s += w[t];
                }
                for (size_t t = 0; t < dh_k.size(); ++t) {
                    int k = dh_k[t];
                    acc(i + FA_DI[k], j + FA_DJ[k]) += total * w[t] / s;
                }
            } else { // dinf: steepest of the 8 triangular facets, split between its 2 bounding directions
                double best_slope = -1e-12;
                int best_kc = -1, best_kd = -1;
                double best_fc = 0.0, best_fd = 0.0;
                for (int f = 0; f < 8; ++f) {
                    int kA = f, kB = (f + 1) % 8;
                    int kc, kd; // cardinal (distance 1), diagonal (distance 1 from cardinal, sqrt2 from centre)
                    if (kA % 2 == 0) { kc = kA; kd = kB; } else { kc = kB; kd = kA; }
                    int nic = i + FA_DI[kc], njc = j + FA_DJ[kc];
                    int nid = i + FA_DI[kd], njd = j + FA_DJ[kd];
                    if (!valid(nic, njc) || !valid(nid, njd)) continue;
                    if (!same_basin(i, j, nic, njc) || !same_basin(i, j, nid, njd)) continue;
                    double ec = dm2(nic, njc), ed = dm2(nid, njd);
                    double s1 = d0 - ec, s2 = ec - ed;
                    if (s1 < 0 && s2 < 0) continue; // facet entirely uphill
                    double slope, r;
                    if (s1 <= 0) {
                        if (ed >= d0) continue;
                        slope = (d0 - ed) / std::sqrt(2.0);
                        r = M_PI / 4.0;
                    } else {
                        r = std::atan2(s2, s1);
                        if (r < 0) { r = 0.0; slope = s1; }
                        else if (r > M_PI / 4.0) {
                            r = M_PI / 4.0;
                            slope = (ed < d0) ? (d0 - ed) / std::sqrt(2.0) : s1;
                        } else {
                            slope = std::sqrt(s1 * s1 + s2 * s2);
                        }
                    }
                    if (slope > best_slope) {
                        best_slope = slope;
                        double frac_diag = r / (M_PI / 4.0);
                        best_kc = kc; best_kd = kd;
                        best_fc = 1.0 - frac_diag; best_fd = frac_diag;
                    }
                }
                if (best_kc >= 0) {
                    if (best_fc > 0) acc(i + FA_DI[best_kc], j + FA_DJ[best_kc]) += total * best_fc;
                    if (best_fd > 0) acc(i + FA_DI[best_kd], j + FA_DJ[best_kd]) += total * best_fd;
                } else {
                    // no facet had both bounding neighbours available (e.g.
                    // one side missing/cross-basin) -- fall back to the
                    // single steepest downhill neighbour, same as d8.
                    int best = 0;
                    for (size_t t = 1; t < dh_k.size(); ++t)
                        if (dh_slope[t] > dh_slope[best]) best = (int)t;
                    int k = dh_k[best];
                    acc(i + FA_DI[k], j + FA_DJ[k]) += total;
                }
            }
            continue;
        }

        // flat handling: same-elevation, same-basin neighbours strictly
        // closer to a reachable exit than this cell.
        int d_here = dist[idx(i, j)];
        if (d_here < 0) continue; // no exit reachable at all -- terminal (unfilled pit / isolated basin)
        std::vector<int> flat_k;
        for (int k = 0; k < 8; ++k) {
            int ni = i + FA_DI[k], nj = j + FA_DJ[k];
            if (!valid(ni, nj) || !same_basin(i, j, ni, nj)) continue;
            if (dm2(ni, nj) == d0 && dist[idx(ni, nj)] >= 0 && dist[idx(ni, nj)] < d_here)
                flat_k.push_back(k);
        }
        if (!flat_k.empty()) {
            if (method == 0) {
                int best = flat_k[0];
                for (int k : flat_k)
                    if (dist[idx(i + FA_DI[k], j + FA_DJ[k])] < dist[idx(i + FA_DI[best], j + FA_DJ[best])]) best = k;
                acc(i + FA_DI[best], j + FA_DJ[best]) += total;
            } else {
                double share = total / (double)flat_k.size();
                for (int k : flat_k) acc(i + FA_DI[k], j + FA_DJ[k]) += share;
            }
        }
        // else: terminal sink (unfilled pit, or basin with no reachable exit under isolation).
    }

    return acc;
}

// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for flow path tracing ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// flowpath_cpp(): traces a steepest-descent downslope path from each of a
// set of starting ("source") cells -- the same steepest-descent rule
// flowacc_cpp()'s own method 0 ("d8") routes with (diagonal drop divided
// by sqrt(2) before comparing to an orthogonal drop, so a diagonal
// neighbour only wins if it's genuinely steeper, not just farther away),
// reusing this file's own FA_DI/FA_DJ neighbour-offset tables (even
// indices = the 4 cardinal directions, used alone for d4; all 8 used for
// d8). Deciding *which* cells are sources (the point itself; same-
// elevation neighbours; a whole land-cover-delineated water body) is done
// entirely on the R side (flowpath()), which can lean on terra for the
// water-body case -- this function only ever traces from whatever source
// list it's handed.
//
// Every source cell is marked 2, then traced separately; every other cell
// any of those traces passes through is marked 1. A trace stops as soon
// as it reaches a cell already marked (already covered by another
// source's route down -- the rest of the way is already known), an edge,
// an NA cell, or a genuine pit (no strictly lower neighbour).
//
// dm                 - elevation matrix, NA (R's NA_real_, i.e. an
//                       ordinary NaN once passed to C++ -- checked with
//                       std::isnan(), matching this package's usual
//                       convention for a plain elevation matrix)
//                       representing no-data.
// src_rows/src_cols  - 0-based row/col of every source cell (already
//                       validated as in-bounds and non-NA on the R side).
// d8                 - TRUE: 8-neighbour search; FALSE: 4-neighbour (only
//                       the cardinal directions).
// ------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix flowpath_cpp(NumericMatrix& dm, IntegerVector src_rows, IntegerVector src_cols, bool d8) {
    int rows = dm.nrow();
    int cols = dm.ncol();
    int nk = d8 ? 8 : 4;

    IntegerMatrix out(rows, cols);
    std::fill(out.begin(), out.end(), 0);

    // Traces a single steepest-descent path starting at (row, col) --
    // (row, col) itself is assumed already marked by the caller.
    auto trace = [&](int row, int col) {
        int i = row, j = col;
        while (true) {
            double d0 = dm(i, j);
            double best_slope = 0.0; // any real strictly-downhill slope is > 0
            int best_i = -1, best_j = -1;
            for (int kk = 0; kk < nk; ++kk) {
                int k = d8 ? kk : kk * 2; // d4: cardinal (even-index) offsets only
                int ni = i + FA_DI[k], nj = j + FA_DJ[k];
                if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
                if (std::isnan(dm(ni, nj)) || dm(ni, nj) >= d0) continue; // not strictly lower
                double dist  = (FA_DI[k] != 0 && FA_DJ[k] != 0) ? std::sqrt(2.0) : 1.0;
                double slope = (d0 - dm(ni, nj)) / dist;
                if (slope > best_slope) {
                    best_slope = slope;
                    best_i = ni;
                    best_j = nj;
                }
            }
            if (best_i < 0) break; // pit: no strictly lower neighbour -- path ends here
            i = best_i;
            j = best_j;
            if (out(i, j) != 0) break; // merges into an already-traced route
            out(i, j) = 1;
        }
    };

    // Marked in a first pass, all together, before any tracing starts --
    // so every trace below sees the *complete* source set as already
    // claimed (2), not just whichever sources happened to be processed
    // earlier.
    int n_src = src_rows.size();
    for (int s = 0; s < n_src; ++s) out(src_rows[s], src_cols[s]) = 2;
    for (int s = 0; s < n_src; ++s) trace(src_rows[s], src_cols[s]);

    return out;
}

// ------------------------------------------------------------
// flowpath_mfd_cpp(): like flowpath_cpp(), but instead of a single
// steepest-descent trace per source, water starting at each source cell
// is apportioned across *every* downhill neighbour at once, split in
// proportion to slope^FA_MFD_EXPONENT -- exactly flowacc_cpp()'s own MFD
// weighting (method 1), just run "downhill from a fixed source" instead
// of "uphill accumulation from every cell". Source cells start at 1.0
// (the full amount); every cell downhill of a source ends up with
// whatever fraction of that 1.0 actually reaches it, summed over every
// route in -- a "water diffusion" layer rather than a single line.
//
// Cells are processed in strictly descending elevation order (a valid
// topological order here, since every step in this package's flow model
// moves to a strictly lower cell -- no cycles), so by the time a cell is
// processed every contribution it could possibly receive from a higher
// cell has already arrived. A cell with no strictly-lower neighbour (a
// pit) or that sits at the domain edge simply keeps whatever it received
// and passes nothing further on -- same "run fillsinks() first if you
// want flow to continue through small artefact depressions" caveat as
// flowpath_cpp()'s single-path trace. No flat-crossing logic is applied
// here either, for the same reason.
//
// dm                - elevation matrix (NA/no-data as std::isnan(), as
//                      elsewhere in this package).
// src_rows/src_cols - 0-based row/col of every source cell.
// d8                - TRUE: 8-neighbour split; FALSE: 4-neighbour (only
//                      the cardinal directions).
// ------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix flowpath_mfd_cpp(NumericMatrix& dm, IntegerVector src_rows, IntegerVector src_cols, bool d8) {
    int rows = dm.nrow();
    int cols = dm.ncol();
    int nk = d8 ? 8 : 4;

    NumericMatrix val(rows, cols);
    std::fill(val.begin(), val.end(), 0.0);

    int n_src = src_rows.size();
    for (int s = 0; s < n_src; ++s) val(src_rows[s], src_cols[s]) = 1.0;

    // Every real (non-NA) cell, sorted by elevation descending.
    std::vector<int> order;
    order.reserve(rows * cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (!std::isnan(dm(i, j))) order.push_back(i * cols + j);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return dm(a / cols, a % cols) > dm(b / cols, b % cols);
    });

    double w[8];
    int wi[8], wj[8];
    for (int idx : order) {
        int i = idx / cols, j = idx % cols;
        double v = val(i, j);
        if (v <= 0.0) continue; // nothing arrived here (from a source or otherwise) -- nothing to pass on
        double d0 = dm(i, j);

        int nc = 0;
        double wsum = 0.0;
        for (int kk = 0; kk < nk; ++kk) {
            int k = d8 ? kk : kk * 2;
            int ni = i + FA_DI[k], nj = j + FA_DJ[k];
            if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
            if (std::isnan(dm(ni, nj)) || dm(ni, nj) >= d0) continue; // not strictly lower
            double dist  = (FA_DI[k] != 0 && FA_DJ[k] != 0) ? std::sqrt(2.0) : 1.0;
            double slope = (d0 - dm(ni, nj)) / dist;
            double wk = std::pow(slope, FA_MFD_EXPONENT);
            wi[nc] = ni; wj[nc] = nj; w[nc] = wk; wsum += wk;
            ++nc;
        }
        if (nc == 0) continue; // pit or edge -- water is lost here, nothing propagates further

        for (int t = 0; t < nc; ++t)
            val(wi[t], wj[t]) += v * (w[t] / wsum);
    }

    // Safety net: force every source cell back to exactly 1 (in the
    // ordinary case -- disjoint sources sharing the same elevation, as
    // both flowpath()'s same-elevation-neighbour and waterbody source
    // selection guarantee -- no source can ever be strictly downhill of
    // another, so this is a no-op; kept only to guarantee the "source
    // cell(s) = 1" contract exactly even if that assumption is ever
    // violated by a future source-selection option) and cap every other
    // cell at 1 (rounding at the edge of the domain could otherwise push
    // a value fractionally above 1).
    for (int s = 0; s < n_src; ++s) val(src_rows[s], src_cols[s]) = 1.0;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (val(i, j) > 1.0) val(i, j) = 1.0;

    return val;
}

// ------------------------------------------------------------
// flowpath_up_cpp(): the "upstream walk" -- literally flowpath_cpp() with
// the descent comparison flipped, tracing a single steepest-ASCENT path
// from each source cell instead of steepest-descent. Not a hydrological
// contributing-area calculation (that's flowpath_mfd_up_cpp() below) --
// just one representative line running uphill from the point, the same
// way flowpath_cpp() traces one representative line running downhill.
// Ends at a local peak/ridge (no strictly higher neighbour), the domain
// edge, an NA cell, or a cell already reached by another source's walk.
// ------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix flowpath_up_cpp(NumericMatrix& dm, IntegerVector src_rows, IntegerVector src_cols, bool d8) {
    int rows = dm.nrow();
    int cols = dm.ncol();
    int nk = d8 ? 8 : 4;

    IntegerMatrix out(rows, cols);
    std::fill(out.begin(), out.end(), 0);

    auto trace = [&](int row, int col) {
        int i = row, j = col;
        while (true) {
            double d0 = dm(i, j);
            double best_slope = 0.0; // any real strictly-uphill slope is > 0
            int best_i = -1, best_j = -1;
            for (int kk = 0; kk < nk; ++kk) {
                int k = d8 ? kk : kk * 2;
                int ni = i + FA_DI[k], nj = j + FA_DJ[k];
                if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
                if (std::isnan(dm(ni, nj)) || dm(ni, nj) <= d0) continue; // not strictly higher
                double dist  = (FA_DI[k] != 0 && FA_DJ[k] != 0) ? std::sqrt(2.0) : 1.0;
                double slope = (dm(ni, nj) - d0) / dist;
                if (slope > best_slope) {
                    best_slope = slope;
                    best_i = ni;
                    best_j = nj;
                }
            }
            if (best_i < 0) break; // local peak/ridge -- walk ends here
            i = best_i;
            j = best_j;
            if (out(i, j) != 0) break; // merges into an already-walked route
            out(i, j) = 1;
        }
    };

    int n_src = src_rows.size();
    for (int s = 0; s < n_src; ++s) out(src_rows[s], src_cols[s]) = 2;
    for (int s = 0; s < n_src; ++s) trace(src_rows[s], src_cols[s]);

    return out;
}

// ------------------------------------------------------------
// flowpath_mfd_up_cpp(): the full upstream contributing area, graded by
// how much of each cell's own water actually reaches the source/target --
// e.g. "pollution turned up at this point; which upslope cells could it
// plausibly have come from, and how much of a role did each play?" Exact
// reverse of flowpath_mfd_cpp(): same slope^FA_MFD_EXPONENT downhill
// split at every cell, but instead of *pushing* a source's value forward
// onto its downhill neighbours, each cell's own "influence" (fraction of
// its water that ends up at the target) is *pulled* from the influence
// values of its own downhill neighbours -- influence(x) = sum over x's
// downhill neighbours y of share(x -> y) * influence(y).
//
// Solved in a single pass by processing cells in strictly ASCENDING
// elevation order this time (the reverse of flowpath_mfd_cpp()'s
// descending order), since a cell's influence here depends on its
// strictly LOWER downhill neighbours' influence, which must therefore
// already be resolved. Source/target cells are fixed at 1 throughout (not
// recomputed from their own downhill split -- by definition, 100% of
// whatever starts there is already "at" the target).
// ------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix flowpath_mfd_up_cpp(NumericMatrix& dm, IntegerVector src_rows, IntegerVector src_cols, bool d8) {
    int rows = dm.nrow();
    int cols = dm.ncol();
    int nk = d8 ? 8 : 4;

    NumericMatrix val(rows, cols);
    std::fill(val.begin(), val.end(), 0.0);

    LogicalMatrix is_src(rows, cols);
    std::fill(is_src.begin(), is_src.end(), false);
    int n_src = src_rows.size();
    for (int s = 0; s < n_src; ++s) {
        val(src_rows[s], src_cols[s]) = 1.0;
        is_src(src_rows[s], src_cols[s]) = true;
    }

    // Every real (non-NA) cell, sorted by elevation ascending.
    std::vector<int> order;
    order.reserve(rows * cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (!std::isnan(dm(i, j))) order.push_back(i * cols + j);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return dm(a / cols, a % cols) < dm(b / cols, b % cols);
    });

    double w[8];
    int wi[8], wj[8];
    for (int idx : order) {
        int i = idx / cols, j = idx % cols;
        if (is_src(i, j)) continue; // fixed at 1 -- not recomputed

        double d0 = dm(i, j);
        int nc = 0;
        double wsum = 0.0;
        for (int kk = 0; kk < nk; ++kk) {
            int k = d8 ? kk : kk * 2;
            int ni = i + FA_DI[k], nj = j + FA_DJ[k];
            if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
            if (std::isnan(dm(ni, nj)) || dm(ni, nj) >= d0) continue; // not strictly lower
            double dist  = (FA_DI[k] != 0 && FA_DJ[k] != 0) ? std::sqrt(2.0) : 1.0;
            double slope = (d0 - dm(ni, nj)) / dist;
            double wk = std::pow(slope, FA_MFD_EXPONENT);
            wi[nc] = ni; wj[nc] = nj; w[nc] = wk; wsum += wk;
            ++nc;
        }
        if (nc == 0) continue; // pit/edge: this cell's water goes nowhere -- zero influence on the target

        double infl = 0.0;
        for (int t = 0; t < nc; ++t)
            infl += (w[t] / wsum) * val(wi[t], wj[t]);
        val(i, j) = infl;
    }

    for (int s = 0; s < n_src; ++s) val(src_rows[s], src_cols[s]) = 1.0;
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (val(i, j) > 1.0) val(i, j) = 1.0;

    return val;
}
