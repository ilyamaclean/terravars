#include "utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================
// solar.cpp
//
// All solar-radiation and canopy-radiation functions for the
// package (previously split across solar.cpp, clearsky.cpp and
// twostream.cpp -- consolidated into one file, 2026-07-17).
//
// Contents, in order:
//
//   1. Clear-sky radiation utilities
//        clearsky_beam()              internal helper
//        clearsky_cpp()                clear-sky beam irradiance
//        annual_solar_weights_cpp()    annual hourly solar position +
//                                       clear-sky weight
//
//   2. Core solar position / direct-beam index
//        solar_index_cpp()             direct beam index, single instant
//        solar_index_tme_cpp()         direct beam index averaged over
//                                       a time series
//        solar_pos_cpp()                solar position over a time vector
//
//   3. Terrain horizon / sky-view factor
//        horizon_angle_cpp()            terrain horizon elevation angle
//        sky_view_factor_cpp()          sky hemisphere visibility fraction
//
//   4. Two-stream canopy radiation model
//        canopy_kd()                    internal helper
//        twostream_tme_cpp()             full two-stream model, averaged
//                                       over a time vector
//        canopy_diffuse_frac_cpp()       two-stream diffuse transmission
//                                       fraction only
//
//   5. 3D ray-marched vegetation attenuation (used by canopy3d() and by
//      the pai+hgt path of solarindex()/skyview())
//        veg_path_atten_cpp()            direct-beam vegetation
//                                       transmission factor
//        veg_path_diffuse_cpp()          diffuse-sky-direction vegetation
//                                       transmission factor
// ============================================================


// ============================================================
// 1. Clear-sky radiation utilities
//
// Reference: Bird & Hulstrom (1981) simplified clear-sky model.
//            Air mass after Kasten (1965).
// ============================================================

// ------------------------------------------------------------
// clearsky_beam  (internal helper)
//
// Clear-sky direct-beam irradiance on a HORIZONTAL surface (W m⁻²)
// for a given solar zenith (radians).  Returns 0 when sun is below
// the horizon (zen_r >= pi/2).
//
// Inputs:
//   zen_r    solar zenith angle (radians)
//   temp_C   air temperature (°C)    — affects precipitable water
//   relhum   relative humidity (%)   — affects precipitable water
//   pres_hPa surface pressure (hPa)  — affects air mass
// ------------------------------------------------------------
double clearsky_beam(double zen_r,
                             double temp_C,
                             double relhum,
                             double pres_hPa) {

    if (zen_r >= 0.5 * MT_PI) return 0.0;   // night

    double cos_z = std::cos(zen_r);

    // ---- Air mass: Kasten (1965) pressure-corrected ----
    double m = (pres_hPa / 1013.25) *
               35.0 / std::sqrt(1224.0 * cos_z * cos_z + 1.0);

    // ---- Rayleigh transmittance (Bird & Hulstrom eq. 10) ----
    double mp  = std::pow(m, 0.84);
    double m1  = std::pow(m, 1.01);
    double Tr  = std::exp(-0.0903 * mp * (1.0 + m - m1));

    // ---- Ozone transmittance (3.5 cm-atm, typical midlatitude) ----
    double u_o  = 0.35;      // column ozone (cm-atm)
    double uo_m = u_o * m;
    double To   = 1.0 - (0.1611 * uo_m * std::pow(1.0 + 139.48 * uo_m, -0.3035)
                        - 0.002715 * uo_m /
                          (1.0 + 0.044 * uo_m + 3e-4 * uo_m * uo_m));

    // ---- Precipitable water (cm) from temp and relative humidity ----
    //      Saturation vapour pressure by August-Roche-Magnus
    double e_sat = 6.1078 * std::exp(17.27 * temp_C / (temp_C + 237.3));
    double e_act = (relhum / 100.0) * e_sat;      // actual vapour pressure (hPa)
    // Precipitable water column (cm) — Reitan (1963) approximation
    double W = 0.493 * e_act / (temp_C + 273.15);
    W = std::max(W, 0.01);   // floor to avoid log(0) or negative
    double wm  = W * m;
    double Tw  = 1.0 - 2.4959 * wm /
                  (std::pow(1.0 + 79.034 * wm, 0.6828) + 6.385 * wm);

    // ---- Aerosol transmittance (rural, visibility ~23 km) ----
    double Ta = std::pow(0.935, m);

    // ---- Solar constant corrected for Earth-Sun distance ----
    //      (annual mean — close enough for a mean annual index)
    double I0 = 1367.0 * 0.9751;

    // ---- Direct normal irradiance × cos(zen) = beam on horizontal ----
    double Idn = I0 * Tr * To * Tw * Ta;
    return Idn * cos_z;
}


// ------------------------------------------------------------
// clearsky_cpp
//
// Vectorised wrapper around clearsky_beam for use from R.
// Returns clear-sky direct-beam irradiance on a horizontal
// surface (W m⁻²) for each solar zenith angle supplied.
// Night-time inputs (zen >= 90°) return 0.
//
// Parameters:
//   zen_deg  - solar zenith angle(s) in degrees
//   temp_C   - air temperature (°C);       default 15
//   relhum   - relative humidity (%);      default 60
//   pres_hPa - surface pressure (hPa);     default 1013.25
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector clearsky_cpp(Rcpp::NumericVector zen_deg,
                                  double temp_C   = 15.0,
                                  double relhum   = 60.0,
                                  double pres_hPa = 1013.25) {

    int n = zen_deg.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i)
        out[i] = clearsky_beam(d2r(zen_deg[i]), temp_C, relhum, pres_hPa);
    return out;
}


// ------------------------------------------------------------
// annual_solar_weights_cpp
//
// Compute hourly solar zenith, azimuth, and clear-sky direct-beam
// weight (W m⁻²) over one calendar year at a single geographic
// location.  Night-time and below-horizon hours have weight = 0
// and are excluded from the output (daytime rows only).
//
// Parameters:
//   lat      geographic latitude  (degrees, +N)
//   lon      geographic longitude (degrees, +E)
//   year     calendar year to use (any non-leap year typical)
//   temp_C   air temperature for clear-sky model (°C)
//   relhum   relative humidity for clear-sky model (%)
//   pres_hPa surface pressure for clear-sky model (hPa)
//   dt_h     time step (hours); default 0.5 h for good bin coverage
//
// Returns: NumericMatrix with columns [zen_deg, azi_deg, weight]
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix annual_solar_weights_cpp(double lat,
                                              double lon,
                                              int    year      = 2023,
                                              double temp_C    = 15.0,
                                              double relhum    = 60.0,
                                              double pres_hPa  = 1013.25,
                                              double dt_h      = 0.5) {

    // Days in year (handle leap year)
    bool leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    int n_days = leap ? 366 : 365;

    // Days-per-month
    int days_in_month[12] = {31, leap ? 29 : 28, 31, 30, 31, 30,
                              31, 31, 30, 31, 30, 31};

    // Maximum possible rows: n_days * (24 / dt_h); shrink via day/night
    int max_rows = static_cast<int>(n_days * 24.0 / dt_h) + 10;

    std::vector<double> vzen, vazi, vwgt;
    vzen.reserve(max_rows / 2);
    vazi.reserve(max_rows / 2);
    vwgt.reserve(max_rows / 2);

    int month = 1, day = 1;

    for (int d = 0; d < n_days; ++d) {
        // Number of steps in this day
        int n_steps = static_cast<int>(24.0 / dt_h);

        for (int s = 0; s < n_steps; ++s) {
            double hr = s * dt_h;      // UTC hour (0.0 = midnight)

            solmodel pos = solpositionCpp(lat, lon, year, month, day, hr);

            if (pos.zend >= 90.0) continue;   // below horizon

            double w = clearsky_beam(pos.zenr, temp_C, relhum, pres_hPa);
            if (w <= 0.0) continue;

            vzen.push_back(pos.zend);
            vazi.push_back(pos.azid);
            vwgt.push_back(w);
        }

        // Advance day/month
        ++day;
        if (day > days_in_month[month - 1]) {
            day = 1;
            ++month;
        }
    }

    int n = static_cast<int>(vzen.size());
    Rcpp::NumericMatrix out(n, 3);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = vzen[i];
        out(i, 1) = vazi[i];
        out(i, 2) = vwgt[i];
    }
    Rcpp::colnames(out) = Rcpp::CharacterVector::create("zen", "azi", "weight");
    return out;
}


// ============================================================
// 2. Core solar position / direct-beam index
// ============================================================

// ------------------------------------------------------------
// solar_index_cpp
//
// Direct beam solar index for a single solar position.
// Returns cos(incidence angle) on a tilted surface, clamped
// to [0, 1]; 0 when sun is below the flat horizon.
//
// Optional terrain shading: when horizon has the same
// dimensions as slope/aspect, cells where the sun elevation
// is at or below the local terrain horizon angle are set to 0.
//
// Parameters (all angle matrices in degrees):
//   slope        - surface slope
//   aspect       - surface aspect (N=0, clockwise)
//   solar_zenith  - solar zenith angle (scalar)
//   solar_azimuth - solar azimuth angle (scalar, N=0, clockwise)
//   horizon      - horizon elevation angles in the solar azimuth
//                  direction; pass a 0x0 matrix to skip
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix solar_index_cpp(Rcpp::NumericMatrix slope,
                                     Rcpp::NumericMatrix aspect,
                                     double solar_zenith,
                                     double solar_azimuth,
                                     Rcpp::NumericMatrix horizon) {

    int nrows = slope.nrow();
    int ncols = slope.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

    if (solar_zenith >= 90.0) return out;

    bool use_horizon = (horizon.nrow() == nrows && horizon.ncol() == ncols);
    double sun_elev  = 90.0 - solar_zenith;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            if (use_horizon && sun_elev <= horizon(r, c)) continue;
            out(r, c) = solarindexCpp(slope(r, c), aspect(r, c),
                                      solar_zenith, solar_azimuth);
        }
    }
    return out;
}


// ------------------------------------------------------------
// solar_index_tme_cpp
//
// Direct beam solar index averaged over a time series.
// Solar position is computed per cell (lat/lon differ across
// the grid) for each timestep, then averaged.
//
// Night-time steps (sun below horizon) contribute 0 to the
// mean, so averaging over a full day gives the mean daytime
// index weighted by night-time zeros — pass only daytime
// steps if you want a daytime-only mean.
//
// Optional terrain shading: supply a precomputed horizon-angle
// array (nrows x ncols x n_az, column-major R array) built
// from evenly-spaced azimuths (step = 360 / n_az degrees).
// For each timestep the local horizon at the solar azimuth is
// obtained by linear interpolation between the two nearest
// array slices.  Pass an empty (length-0) vector to skip.
//
// Parameters:
//   slope, aspect  - surface geometry matrices (degrees)
//   lat, lon       - per-cell geographic coordinates (degrees)
//   years, months, days - date components for each timestep
//   hours          - UTC decimal hour for each timestep
//   horizon_cube   - precomputed horizon angles (degrees);
//                    R array with dim c(nrows, ncols, n_az)
//   az_step        - azimuth spacing of horizon cube (degrees)
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix solar_index_tme_cpp(Rcpp::NumericMatrix slope,
                                         Rcpp::NumericMatrix aspect,
                                         Rcpp::NumericMatrix lat,
                                         Rcpp::NumericMatrix lon,
                                         Rcpp::IntegerVector years,
                                         Rcpp::IntegerVector months,
                                         Rcpp::IntegerVector days,
                                         Rcpp::NumericVector hours,
                                         Rcpp::NumericVector horizon_cube,
                                         double az_step) {

    int nrows = slope.nrow();
    int ncols = slope.ncol();
    int n_t   = years.size();

    // Determine if terrain shading is active from the cube dimensions
    int n_az  = 0;
    bool shade = false;
    if (horizon_cube.size() > 0) {
        Rcpp::IntegerVector dims = horizon_cube.attr("dim");
        n_az  = dims[2];
        shade = (n_az > 0);
    }

    // Copy read-only vectors to plain arrays for thread-safe access
    // (Rcpp proxy objects are not guaranteed thread-safe for concurrent reads
    //  across different indices, but raw pointers are fine)
    const int*    y_ptr  = years.begin();
    const int*    mo_ptr = months.begin();
    const int*    d_ptr  = days.begin();
    const double* h_ptr  = hours.begin();
    const double* hc_ptr = (shade) ? horizon_cube.begin() : nullptr;

    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            double la  = lat(r, c);
            double lo  = lon(r, c);
            double slp = slope(r, c);
            double asp = aspect(r, c);
            double sum = 0.0;

            for (int t = 0; t < n_t; ++t) {
                solmodel pos = solpositionCpp(la, lo,
                                             y_ptr[t], mo_ptr[t], d_ptr[t],
                                             h_ptr[t]);

                double si = solarindexCpp(slp, asp, pos.zend, pos.azid);

                if (shade && si > 0.0) {
                    double az_frac = pos.azid / az_step;
                    int    k0      = static_cast<int>(std::floor(az_frac)) % n_az;
                    int    k1      = (k0 + 1) % n_az;
                    double alpha   = az_frac - std::floor(az_frac);
                    double hor = (1.0 - alpha) * hc_ptr[r + c * nrows + k0 * nrows * ncols]
                               +        alpha  * hc_ptr[r + c * nrows + k1 * nrows * ncols];
                    double sun_elev = 90.0 - pos.zend;
                    if (sun_elev <= hor) si = 0.0;
                }

                sum += si;
            }

            out(r, c) = (n_t > 0) ? sum / n_t : 0.0;
        }
    }
    return out;
}


// ------------------------------------------------------------
// solar_pos_cpp
//
// Solar position (zenith and azimuth) for a single geographic
// location over a time vector.  Convenience wrapper around
// solpositionCpp used by the R solarindex() function when
// cenloc = TRUE.
//
// Parameters:
//   lat, lon  - geographic coordinates (degrees, scalar)
//   years, months, days, hours  - UTC time components (vectors)
//
// Returns: NumericMatrix [n × 2], columns "zen" (degrees) and
//          "azi" (degrees, N = 0, clockwise).
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix solar_pos_cpp(double lat,
                                   double lon,
                                   Rcpp::IntegerVector years,
                                   Rcpp::IntegerVector months,
                                   Rcpp::IntegerVector days,
                                   Rcpp::NumericVector hours) {

    int n = years.size();
    Rcpp::NumericMatrix out(n, 2);

    for (int t = 0; t < n; ++t) {
        solmodel pos = solpositionCpp(lat, lon,
                                      years[t], months[t], days[t], hours[t]);
        out(t, 0) = pos.zend;
        out(t, 1) = pos.azid;
    }
    Rcpp::colnames(out) = Rcpp::CharacterVector::create("zen", "azi");
    return out;
}


// ============================================================
// 3. Terrain horizon / sky-view factor
// ============================================================

// ------------------------------------------------------------
// horizon_angle_cpp
//
// For each cell, find the maximum elevation angle (degrees)
// to the terrain horizon in a given azimuth direction.
//
// Two sampling modes:
//   n_steps == 0  (linear): checks every cell along the transect
//                 to the raster edge — most accurate.
//   n_steps  > 0  (quadratic): samples at step² cell distances
//                 (1, 4, 9, … n² cells) — faster, physically
//                 motivated (near-field blockers dominate), but
//                 misses terrain at intermediate distances.
//
// Parameters:
//   elev     - elevation matrix (metres)
//   res      - cell resolution (metres, assumed square)
//   azimuth  - azimuth direction to search (degrees, N=0, clockwise)
//   n_steps  - 0 = linear to raster edge; >0 = quadratic steps
//   obs_hgt  - observer height above ground (metres), added to each
//              cell's own elevation before searching -- e.g. for a
//              wind sensor or turbine hub mounted above the surface,
//              which "sees over" more of the upwind terrain than a
//              ground-level observer would. 0 (default) matches the
//              original ground-level behaviour exactly.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix horizon_angle_cpp(Rcpp::NumericMatrix elev,
                                       double res,
                                       double azimuth,
                                       int n_steps = 0,
                                       double obs_hgt = 0.0) {

    int nrows = elev.nrow();
    int ncols = elev.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

    double az_r   = d2r(azimuth);
    double sin_az = std::sin(az_r);
    double cos_az = std::cos(az_r);

    int max_cells = static_cast<int>(
        std::ceil(std::sqrt(static_cast<double>(nrows * nrows + ncols * ncols)))
    );

    // dynamic schedule: search length varies per cell (breaks at edge)
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            double z0      = elev(r, c);
            if (std::isnan(z0)) { out(r, c) = NA_REAL; continue; }
            z0 += obs_hgt;
            double max_tan = 0.0;

            if (n_steps == 0) {
                for (int step = 1; step <= max_cells; ++step) {
                    int rr = r - static_cast<int>(std::round(cos_az * step));
                    int cc = c + static_cast<int>(std::round(sin_az * step));
                    if (!in_grid(rr, cc, nrows, ncols)) break;
                    double dz = elev(rr, cc) - z0;
                    double tan_elev = dz / (step * res);
                    if (tan_elev > max_tan) max_tan = tan_elev;
                }
            } else {
                for (int step = 1; step <= n_steps; ++step) {
                    int dist_cells = step * step;
                    int rr = r - static_cast<int>(std::round(cos_az * dist_cells));
                    int cc = c + static_cast<int>(std::round(sin_az * dist_cells));
                    if (!in_grid(rr, cc, nrows, ncols)) break;
                    double dz = elev(rr, cc) - z0;
                    double tan_elev = dz / (dist_cells * res);
                    if (tan_elev > max_tan) max_tan = tan_elev;
                }
            }

            out(r, c) = r2d(std::atan(max_tan));
        }
    }
    return out;
}


// ------------------------------------------------------------
// sky_view_factor_cpp
//
// Sky-view factor as mean of cos²(horizon_angle) over N
// evenly-spaced azimuth directions.
// Expects a 3-D R array of horizon angles (rows x cols x n_az)
// in degrees, passed as a NumericVector with dim attribute.
//
//   SVF = (1/N) * sum_i [ cos²(h_i) ]
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix sky_view_factor_cpp(Rcpp::NumericVector horizon_cube) {

    Rcpp::IntegerVector dims = horizon_cube.attr("dim");
    int nrows = dims[0];
    int ncols = dims[1];
    int n_az  = dims[2];

    const double* hc_ptr = horizon_cube.begin();

    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            double sum = 0.0;
            for (int k = 0; k < n_az; ++k) {
                double h  = d2r(hc_ptr[r + c * nrows + k * nrows * ncols]);
                double ch = std::cos(h);
                sum += ch * ch;
            }
            out(r, c) = sum / n_az;
        }
    }
    return out;
}


// ============================================================
// 4. Two-stream canopy radiation model
// ============================================================


// ------------------------------------------------------------
// canopy_kd  (internal helper)
//
// Slope-corrected direct-beam canopy extinction coefficient.
// k = K_be(zenith, x) (see kBeCpp() in utils.h) re-expressed via
// G(zenith) = k * cos(zenith) and divided by the target surface's
// actual incidence cosine `si` instead of cos(zenith) -- G stays
// pinned to the true solar zenith (leaf orientation is referenced to
// true vertical), while the path-length elongation factor uses the
// slope-corrected incidence cosine.
//   zen_r  - solar zenith in radians
//   xv     - leaf angle distribution parameter (see kBeCpp())
//   si     - solar index on slope (cos of incidence angle)
// ------------------------------------------------------------
static double canopy_kd(double zen_r, double xv, double si) {
    double cos_z = std::cos(zen_r);
    double kd = (si > 1e-10) ? kBeCpp(zen_r, xv) * cos_z / si : 6000.0;
    return std::min(kd, 6000.0);
}


// ------------------------------------------------------------
// twostream_tme_cpp
//
// Two-stream model averaged (or single-step) over a time vector.
// Solar position is computed per cell from lat/lon matrices.
// Night-time steps contribute 0 to direct/scatter accumulators;
// the divisor is always n_t so averages are consistent with
// solarindex() convention.
//
// fdifdown is solar-angle-independent and returned as a single
// value per cell (no temporal averaging needed).
//
// A cell whose pai/lref/ltra/xmat/gref OR slope/aspect is NA gets NA in
// all three outputs -- notably this includes any land cell along a
// coastline (or beside any other NA/no-data gap) whose slope/aspect came
// out NA from terra::terrain()'s 3x3-neighbourhood requirement, when the
// caller hasn't pre-filled dtm's NA gaps via fillna. Such a cell's true
// geometry simply isn't known, so it's masked rather than silently
// computed from a fallback.
//
// Parameters:
//   slope, aspect  - terrain matrices (degrees)
//   pai            - plant area index matrix (m2 m-2)
//   lref, ltra     - leaf reflectance / transmittance matrices
//   xmat           - leaf angle distribution matrix (1 = spherical)
//   gref           - ground reflectance matrix
//   lat, lon       - per-cell geographic coordinates (degrees)
//   years, months, days, hours  - UTC time components
//   pai_horiz      - TRUE if `pai` is quoted per unit horizontal ground
//                    area (the usual LAI/PAI convention), FALSE if per unit
//                    sloped ground-surface area. Both streams are evaluated
//                    against a single shared canopy depth, pa_h, expressed
//                    per unit *horizontal* area (matching the diffuse
//                    two-stream solution's own horizontally-homogeneous-
//                    layer assumption): pa_h = pai when pai_horiz is TRUE,
//                    or pai / cos(slope) when FALSE (a fixed physical
//                    amount of leaf spread over a tilted patch is
//                    lower-density per sloped area than per horizontal area
//                    by cos(slope), so dividing a sloped-convention value
//                    by cos(slope) recovers its horizontal equivalent).
//                    The direct beam's own extinction coefficient, kd (from
//                    canopy_kd(), natively a per-*sloped*-area rate, via its
//                    cos(zen)/si path-length factor), is rescaled to
//                    kd_eff = kd * cos(slope) so that it decays the shared
//                    depth pa_h at the rate reproducing its native
//                    per-sloped behaviour -- kd_eff * pa_h is algebraically
//                    identical to the old per-stream-depth kd * pa_s either
//                    way, so fdirdown (S2) is numerically unchanged by this.
//                    What matters is that kd_eff, not the raw kd, is then
//                    used *consistently* throughout the rest of the
//                    scattered-flux particular solution (ss, sstr, sig, p8,
//                    v3, p9, p10) rather than only inside S2's own exponent.
//                    That consistency is required because that solution
//                    divides by sig = kd^2 - h^2, a *removable* singularity
//                    at kd == h only when both streams share one depth and
//                    one rate variable throughout -- feeding the two
//                    streams different depths (as an earlier version of
//                    this function did, converting only S2's own pai) turns
//                    that removable singularity into a real one, producing
//                    wildly inflated fdirscat values on any sloped cell
//                    whose solar geometry happens to put kd close to h.
//                    Leaf angles are unaffected -- G stays pinned to true
//                    solar zenith throughout.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List twostream_tme_cpp(Rcpp::NumericMatrix slope,
                               Rcpp::NumericMatrix aspect,
                               Rcpp::NumericMatrix pai,
                               Rcpp::NumericMatrix lref,
                               Rcpp::NumericMatrix ltra,
                               Rcpp::NumericMatrix xmat,
                               Rcpp::NumericMatrix gref,
                               Rcpp::NumericMatrix lat,
                               Rcpp::NumericMatrix lon,
                               Rcpp::IntegerVector years,
                               Rcpp::IntegerVector months,
                               Rcpp::IntegerVector days,
                               Rcpp::NumericVector hours,
                               bool pai_horiz) {

    int nrows = slope.nrow();
    int ncols = slope.ncol();
    int n_t   = years.size();

    Rcpp::NumericMatrix out_dir (nrows, ncols);   // zero-initialised
    Rcpp::NumericMatrix out_dif (nrows, ncols);
    Rcpp::NumericMatrix out_scat(nrows, ncols);

    // Raw pointers for thread-safe read access inside parallel region
    const int*    y_ptr  = years.begin();
    const int*    mo_ptr = months.begin();
    const int*    d_ptr  = days.begin();
    const double* h_ptr  = hours.begin();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {

            double pa = pai (r, c);
            double lr = lref(r, c);
            double lt = ltra(r, c);
            double xv = xmat(r, c);
            double gr = gref(r, c);
            double sl = slope (r, c);
            double as = aspect(r, c);

            // NA propagation -- includes slope/aspect. A land cell whose
            // 3x3 neighbourhood touches an NA/no-data cell (e.g. anywhere
            // along a coastline, whenever fillna = FALSE) gets NA slope
            // and aspect from terra::terrain() even though its own
            // elevation and pai are perfectly valid. Previously this fell
            // through to canopy_kd()'s grazing-incidence sentinel
            // (kd = 6000): that correctly drives fdirdown to ~0, but does
            // NOT drive fdirscat to 0 the same way, so cells with genuinely
            // unknown geometry were silently getting an inflated fdirscat
            // (worst at low pai) instead of a clean NA. Masking here
            // instead means all three outputs are NA at any cell whose
            // geometry isn't actually known -- callers who want these
            // cells filled in should pass fillna = TRUE, which supplies a
            // real (if approximate) slope/aspect instead of leaving NA.
            if (std::isnan(pa) || std::isnan(lr) || std::isnan(lt) ||
                std::isnan(xv) || std::isnan(gr) ||
                std::isnan(sl) || std::isnan(as)) {
                out_dir (r, c) = NA_REAL;
                out_dif (r, c) = NA_REAL;
                out_scat(r, c) = NA_REAL;
                continue;
            }

            // Shared canopy depth for both streams, per unit horizontal
            // area (see pai_horiz in the header). cos(slope) converts a
            // sloped-convention input up to its horizontal equivalent;
            // clamped so the /cosb path can't blow up on genuinely
            // near-vertical real terrain (slope/aspect NA is now handled
            // above, so this is no longer also standing in for that case).
            double cosb = std::cos(d2r(sl));
            if (cosb < 0.05) cosb = 0.05;
            double pa_h = pai_horiz ? pa : pa / cosb;

            // Time-invariant optical parameters
            double om   = lr + lt;
            double a    = 1.0 - om;
            double del_ = lr - lt;
            double mla  = 9.65 * std::pow(3.0 + xv, -1.65);
            if (mla > M_PI / 2.0) mla = M_PI / 2.0;
            double J    = std::cos(mla) * std::cos(mla);
            double gma  = 0.5 * (om + J * del_);
            double h    = std::sqrt(a * a + 2.0 * a * gma);
            double S1   = std::exp(-h * pa_h);

            // No canopy: direct fraction = fraction of sun-up steps,
            // diffuse fraction = 1, scatter = 0
            if (pa <= 0.0) {
                double acc_dir = 0.0;
                for (int t = 0; t < n_t; ++t) {
                    solmodel pos = solpositionCpp(lat(r, c), lon(r, c),
                                                  y_ptr[t], mo_ptr[t], d_ptr[t],
                                                  h_ptr[t]);
                    if (pos.zend < 90.0) acc_dir += 1.0;
                }
                out_dir (r, c) = acc_dir / n_t;
                out_dif (r, c) = 1.0;
                out_scat(r, c) = 0.0;
                continue;
            }

            // Time-invariant diffuse parameters
            double u2  = a + gma * (1.0 - gr);
            double D2  = (u2 + h) / S1 - (u2 - h) * S1;
            if (std::abs(D2) < 1e-10)
                D2 = (D2 >= 0.0) ? 1e-10 : -1e-10;

            // fdifdown = 2h/D2  (solar-angle independent)
            out_dif(r, c) = 2.0 * h / D2;

            double acc_dir  = 0.0;
            double acc_scat = 0.0;

            for (int t = 0; t < n_t; ++t) {
                solmodel pos = solpositionCpp(lat(r, c), lon(r, c),
                                              y_ptr[t], mo_ptr[t], d_ptr[t],
                                              h_ptr[t]);

                if (pos.zend >= 90.0) continue;

                double si     = solarindexCpp(sl, as, pos.zend, pos.azid);
                double zen_r  = d2r(pos.zend);
                double kd_raw = canopy_kd(zen_r, xv, si);
                // Rescaled to decay the shared depth pa_h at the rate that
                // reproduces kd_raw's native per-sloped behaviour (see
                // pai_horiz in the header). `kd` (not kd_raw) is what's
                // used everywhere below for the rest of this block.
                double kd    = kd_raw * cosb;
                double S2    = std::exp(-kd * pa_h);

                acc_dir += S2;

                double ss   = 0.5 * (om + J * del_ / kd) * kd;
                double sstr = om * kd - ss;

                double sig = kd * kd + gma * gma - (a + gma) * (a + gma);
                if (std::abs(sig) < 1e-10)
                    sig = (sig >= 0.0) ? 1e-10 : -1e-10;

                double p8  = sstr * (a + gma + kd) - gma * ss;
                double v3  = (sstr + gma * gr - (p8 / -sig) * (u2 - kd)) * S2;
                double p9  = (-1.0 / D2) * ((p8 / (-sig * S1)) * (u2 + h) + v3);
                double p10 = ( 1.0 / D2) * (((p8 * S1) / -sig) * (u2 - h) + v3);

                acc_scat += (p8 / -sig) * S2 + p9 * S1 + p10 / S1;
            }

            out_dir (r, c) = acc_dir  / n_t;
            out_scat(r, c) = acc_scat / n_t;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("fdirdown") = out_dir,
        Rcpp::Named("fdifdown") = out_dif,
        Rcpp::Named("fdirscat") = out_scat
    );
}


// ------------------------------------------------------------
// canopy_diffuse_frac_cpp
//
// Two-stream diffuse transmission fraction: the fraction of
// above-canopy isotropic sky diffuse radiation reaching the ground,
// per cell. Mirrors twostream_tme_cpp()'s fdifdown calculation exactly
// (kept in sync manually -- duplicated rather than factored out of
// twostream_tme_cpp so as not to disturb that already-verified
// per-timestep loop, which reuses the same intermediate quantities
// for its direct/scatter accumulators).
//
// Solar-angle-independent -- no lat/lon/time inputs needed. Also
// ignores terrain slope/aspect entirely (as does twostream()'s own
// fdifdown): it assumes a laterally-infinite, horizontally-homogeneous
// canopy layer, which is what makes the two-stream equations solvable
// in closed form. skyview() reintroduces terrain/horizon effects by
// multiplying this per-cell scalar against its existing (unrelated)
// direction-aware terrain-only sky-view factor -- see skyview() in
// R/solar.R for that combination and its rationale.
//
// Parameters:
//   pai            - plant area index matrix (m2 m-2)
//   lref, ltra     - leaf reflectance / transmittance matrices
//   xmat           - leaf angle distribution matrix (1 = spherical)
//   gref           - ground reflectance matrix
//
// Returns 1.0 (full transmission, no canopy) where pai <= 0; NA where
// any input is NA/NaN at that cell.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix canopy_diffuse_frac_cpp(Rcpp::NumericMatrix pai,
                                             Rcpp::NumericMatrix lref,
                                             Rcpp::NumericMatrix ltra,
                                             Rcpp::NumericMatrix xmat,
                                             Rcpp::NumericMatrix gref) {

    int nrows = pai.nrow();
    int ncols = pai.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {

            double pa = pai (r, c);
            double lr = lref(r, c);
            double lt = ltra(r, c);
            double xv = xmat(r, c);
            double gr = gref(r, c);

            // NA propagation
            if (std::isnan(pa) || std::isnan(lr) || std::isnan(lt) ||
                std::isnan(xv) || std::isnan(gr)) {
                out(r, c) = NA_REAL;
                continue;
            }

            // No canopy: full transmission
            if (pa <= 0.0) {
                out(r, c) = 1.0;
                continue;
            }

            double om   = lr + lt;
            double a    = 1.0 - om;
            double del_ = lr - lt;
            double mla  = 9.65 * std::pow(3.0 + xv, -1.65);
            if (mla > M_PI / 2.0) mla = M_PI / 2.0;
            double J    = std::cos(mla) * std::cos(mla);
            double gma  = 0.5 * (om + J * del_);
            double h    = std::sqrt(a * a + 2.0 * a * gma);
            double S1   = std::exp(-h * pa);

            double u2  = a + gma * (1.0 - gr);
            double D2  = (u2 + h) / S1 - (u2 - h) * S1;
            if (std::abs(D2) < 1e-10)
                D2 = (D2 >= 0.0) ? 1e-10 : -1e-10;

            // fdifdown = 2h/D2  (solar-angle independent)
            out(r, c) = 2.0 * h / D2;
        }
    }

    return out;
}


// ============================================================
// 5. 3D ray-marched vegetation attenuation
// ============================================================

// ------------------------------------------------------------
// veg_path_atten_cpp
//
// Direct-beam vegetation transmission factor accounting for a beam
// travelling diagonally across neighbouring pixels' canopies, not just
// straight down through the target pixel's own column (the
// approximation used by solarindex()'s pai-only, no-hgt path).
//
// Each pixel is modelled as a vertical "canopy box" from the ground
// (elev) up to elev + hgt, with leaf area distributed uniformly by
// height inside that box: u_L = pai / hgt per unit height. This is the
// vertically-uniform-foliage simplification of a true 3D voxel model --
// a pixel's canopy occupies only the column strictly above its own
// footprint, never extending sideways into a neighbour's airspace.
//
// The beam is ray-marched outward from each origin cell toward the sun
// using exactly the same step schedule as horizon_angle_cpp(): n_steps
// == 0 samples every cell out to the raster edge (exact, given this
// ray-sampling scheme); n_steps > 0 uses quadratic step^2 spacing for
// speed, at the cost of missing intermediate canopy further out. Given
// the vertically-uniform-column simplification is itself already
// approximate (no true 3D voxels), the coarser quadratic mode is a
// reasonable trade-off rather than a large extra loss of fidelity.
//
// Sample k = 0, 1, 2, ... corresponds to horizontal distance
// dist(k) = 0 for k = 0, else k (linear) or k^2 (quadratic) cell
// widths. The interval [dist(k), dist(k+1)] uses the *near* endpoint's
// cell (i.e. cell(k)) to represent that interval's canopy properties --
// so k = 0 is always the origin cell itself (dist(0) = 0), meaning the
// origin's own canopy is captured exactly regardless of n_steps, and
// only increasingly distant neighbours are coarsened under quadratic
// spacing.
//
// For each interval, the ray's vertical rise over that interval is
// intersected with the sampled cell's canopy box; the overlap (a true
// vertical distance) times that cell's own u_L and its own K_be(zenith,
// x) (see kBeCpp() in utils.h -- each traversed cell's own leaf angle
// distribution applies, not the origin's) gives that interval's optical
// depth contribution. Marching stops once the ray has risen above the
// highest canopy top anywhere in the raster (nothing further out can
// contribute), once it exits the grid, or once the raster-diagonal
// safety bound is reached.
//
// Missing data (NaN elevation/hgt/pai/x) at a *sampled neighbour* cell
// is treated as "no canopy known there" and skipped (contributes zero
// optical depth to that interval) rather than propagating NA through
// the whole accumulated sum. NaN elevation at the *origin* cell itself
// propagates NA to that cell's output, matching horizon_angle_cpp()'s
// convention.
//
// Parameters:
//   elev    - ground elevation (m)
//   hgt     - canopy height above local ground (m); <= 0 or NaN treated
//             as no canopy (u_L = 0) at that cell
//   pai     - plant area index (m2 m-2)
//   xmat    - leaf angle distribution parameter (Campbell & Norman
//             ellipsoidal convention; see kBeCpp())
//   res     - cell resolution (m, assumed square)
//   zenith  - solar zenith angle (degrees, scalar -- one solar position
//             per call, like solar_index_cpp())
//   azimuth - solar azimuth angle (degrees, N = 0, clockwise)
//   n_steps - 0 = linear stepping to the raster edge; > 0 = quadratic
//             step^2 spacing (same speed/accuracy tradeoff as
//             horizon())
//
// Returns: NumericMatrix of vegetation transmission factors in [0, 1]
// (exp(-total optical depth) per cell) -- multiply elementwise into the
// already-computed geometric solar index.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix veg_path_atten_cpp(Rcpp::NumericMatrix elev,
                                        Rcpp::NumericMatrix hgt,
                                        Rcpp::NumericMatrix pai,
                                        Rcpp::NumericMatrix xmat,
                                        double res,
                                        double zenith,
                                        double azimuth,
                                        int n_steps = 0) {

    int nrows = elev.nrow();
    int ncols = elev.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

    if (zenith >= 90.0) return out;   // sun below flat horizon: no beam at all

    double zen_r = d2r(zenith);
    double cos_z = std::cos(zen_r);
    double sin_z = std::sin(zen_r);
    // Vertical rise per unit horizontal distance travelled = cot(zenith).
    double cot_zen = (sin_z > 1e-10) ? cos_z / sin_z : 1e10;

    double az_r   = d2r(azimuth);
    double sin_az = std::sin(az_r);
    double cos_az = std::cos(az_r);

    // Raster-wide maximum canopy-top elevation, used to stop marching
    // early once a beam has risen above any possible canopy anywhere,
    // regardless of which cell it currently happens to be over.
    double max_top   = -std::numeric_limits<double>::infinity();
    bool   any_canopy = false;
    for (int rr = 0; rr < nrows; ++rr) {
        for (int cc = 0; cc < ncols; ++cc) {
            double e = elev(rr, cc);
            double h = hgt(rr, cc);
            if (!std::isnan(e) && !std::isnan(h) && h > 0.0) {
                double top = e + h;
                if (top > max_top) max_top = top;
                any_canopy = true;
            }
        }
    }

    int max_cells = static_cast<int>(
        std::ceil(std::sqrt(static_cast<double>(nrows * nrows + ncols * ncols)))
    );

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {

            double z0 = elev(r, c);
            if (std::isnan(z0)) { out(r, c) = NA_REAL; continue; }
            if (!any_canopy)    { out(r, c) = 1.0;     continue; }

            double optical_depth = 0.0;
            double prev_z         = z0;   // ray height accounted for so far

            for (int k = 0; ; ++k) {

                int this_step = (n_steps == 0) ? k : k * k;
                int rr = r - static_cast<int>(std::round(cos_az * this_step));
                int cc = c + static_cast<int>(std::round(sin_az * this_step));
                if (!in_grid(rr, cc, nrows, ncols)) break;

                if (prev_z > max_top) break;          // above all possible canopy now
                if (n_steps > 0 && k >= n_steps) break;

                int next_step   = (n_steps == 0) ? (k + 1) : (k + 1) * (k + 1);
                double dist_next = next_step * res;
                double z_next    = z0 + dist_next * cot_zen;

                double eg = elev(rr, cc);
                double hg = hgt (rr, cc);
                double pg = pai (rr, cc);
                double xg = xmat(rr, cc);

                if (!std::isnan(eg) && !std::isnan(hg) && hg > 0.0 &&
                    !std::isnan(pg) && pg > 0.0 && !std::isnan(xg)) {

                    double top = eg + hg;
                    double lo  = std::max(prev_z, eg);
                    double hi  = std::min(z_next, top);
                    if (hi > lo) {
                        double overlap = hi - lo;
                        double u_L     = pg / hg;
                        double kbe     = kBeCpp(zen_r, xg);
                        optical_depth += kbe * u_L * overlap;
                    }
                }
                // else: unknown/no canopy at this sampled cell -- skip,
                // contributes zero optical depth (see docstring above)

                prev_z = z_next;
                if (next_step >= max_cells) break;
            }

            out(r, c) = std::exp(-optical_depth);
        }
    }

    return out;
}


// ------------------------------------------------------------
// veg_path_diffuse_cpp
//
// Diffuse-sky-direction vegetation transmission factor, for a single
// (elevation, azimuth) sky direction, accounting for a ray travelling
// diagonally across neighbouring pixels' canopies -- the diffuse
// counterpart of veg_path_atten_cpp(), sharing an identical ray-marching
// scheme (same vertically-uniform-canopy-box simplification, same
// n_steps step schedule, same near-endpoint interval convention -- see
// veg_path_atten_cpp()'s docstring for the full rationale, not repeated
// here). Deliberately a separate function rather than a shared/
// parameterised core with veg_path_atten_cpp(), so as not to disturb
// that already-verified ray-marcher; the two are kept structurally in
// sync manually.
//
// The only physical difference from veg_path_atten_cpp() is the
// per-cell extinction coefficient: instead of the zenith/leaf-angle
// dependent K_be(zenith, x) (kBeCpp()), this uses Jan Goudriaan's
// simplified Beer-Lambert diffuse approximation of the two-stream
// model, sqrt(alpha) * K_diffuse with K_diffuse = 1 and
// alpha = 1 - leaf_reflectance - leaf_transmittance (leaf absorptivity)
// -- a fixed per-cell value with no zenith-angle dependence, unlike the
// direct-beam K_be(zenith, x). `zenith` here is 90 - elevation for the
// sky direction being sampled, not a solar zenith angle; it still
// determines the ray's actual 3D geometry (how fast it rises per unit
// horizontal distance), just not the per-step extinction multiplier.
//
// Parameters:
//   elev         - ground elevation (m)
//   hgt          - canopy height above local ground (m); <= 0 or NaN
//                  treated as no canopy (u_L = 0) at that cell
//   pai          - plant area index (m2 m-2)
//   lref, ltra   - leaf reflectance / transmittance (0-1); alpha at each
//                  traversed cell is 1 - lref - ltra, clamped to >= 0
//   res          - cell resolution (m, assumed square)
//   zenith       - 90 - elevation of the sky direction being sampled
//                  (degrees, scalar -- one direction per call, like
//                  veg_path_atten_cpp())
//   azimuth      - azimuth of the sky direction (degrees, N = 0,
//                  clockwise)
//   n_steps      - 0 = linear stepping to the raster edge; > 0 =
//                  quadratic step^2 spacing (same as horizon())
//
// Returns: NumericMatrix of diffuse vegetation transmission factors in
// [0, 1] for this single sky direction -- see skyview() (R/solar.R) for
// how multiple directions are combined into a hemisphere-integrated
// diffuse sky-view factor.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix veg_path_diffuse_cpp(Rcpp::NumericMatrix elev,
                                          Rcpp::NumericMatrix hgt,
                                          Rcpp::NumericMatrix pai,
                                          Rcpp::NumericMatrix lref,
                                          Rcpp::NumericMatrix ltra,
                                          double res,
                                          double zenith,
                                          double azimuth,
                                          int n_steps = 0) {

    int nrows = elev.nrow();
    int ncols = elev.ncol();
    Rcpp::NumericMatrix out(nrows, ncols);   // zero-initialised

    if (zenith >= 90.0) return out;   // direction at/below horizontal: no sky there

    double zen_r = d2r(zenith);
    double cos_z = std::cos(zen_r);
    double sin_z = std::sin(zen_r);
    double cot_zen = (sin_z > 1e-10) ? cos_z / sin_z : 1e10;

    double az_r   = d2r(azimuth);
    double sin_az = std::sin(az_r);
    double cos_az = std::cos(az_r);

    double max_top   = -std::numeric_limits<double>::infinity();
    bool   any_canopy = false;
    for (int rr = 0; rr < nrows; ++rr) {
        for (int cc = 0; cc < ncols; ++cc) {
            double e = elev(rr, cc);
            double h = hgt(rr, cc);
            if (!std::isnan(e) && !std::isnan(h) && h > 0.0) {
                double top = e + h;
                if (top > max_top) max_top = top;
                any_canopy = true;
            }
        }
    }

    int max_cells = static_cast<int>(
        std::ceil(std::sqrt(static_cast<double>(nrows * nrows + ncols * ncols)))
    );

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {

            double z0 = elev(r, c);
            if (std::isnan(z0)) { out(r, c) = NA_REAL; continue; }
            if (!any_canopy)    { out(r, c) = 1.0;     continue; }

            double optical_depth = 0.0;
            double prev_z         = z0;

            for (int k = 0; ; ++k) {

                int this_step = (n_steps == 0) ? k : k * k;
                int rr = r - static_cast<int>(std::round(cos_az * this_step));
                int cc = c + static_cast<int>(std::round(sin_az * this_step));
                if (!in_grid(rr, cc, nrows, ncols)) break;

                if (prev_z > max_top) break;
                if (n_steps > 0 && k >= n_steps) break;

                int next_step   = (n_steps == 0) ? (k + 1) : (k + 1) * (k + 1);
                double dist_next = next_step * res;
                double z_next    = z0 + dist_next * cot_zen;

                double eg = elev(rr, cc);
                double hg = hgt (rr, cc);
                double pg = pai (rr, cc);
                double lr = lref(rr, cc);
                double lt = ltra(rr, cc);

                if (!std::isnan(eg) && !std::isnan(hg) && hg > 0.0 &&
                    !std::isnan(pg) && pg > 0.0 &&
                    !std::isnan(lr) && !std::isnan(lt)) {

                    double top = eg + hg;
                    double lo  = std::max(prev_z, eg);
                    double hi  = std::min(z_next, top);
                    if (hi > lo) {
                        double overlap = hi - lo;
                        double u_L     = pg / hg;
                        double alpha   = std::max(1.0 - lr - lt, 0.0);
                        double k_eff   = std::sqrt(alpha);   // * K_diffuse (= 1)
                        optical_depth += k_eff * u_L * overlap;
                    }
                }
                // else: unknown/no canopy at this sampled cell -- skip,
                // contributes zero optical depth (see docstring above)

                prev_z = z_next;
                if (next_step >= max_cells) break;
            }

            out(r, c) = std::exp(-optical_depth);
        }
    }

    return out;
}
