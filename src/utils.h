#pragma once
#include <Rcpp.h>
#include <cmath>

// ============================================================
// Mathematical constants
// ============================================================
static const double MT_PI  = 3.14159265358979323846;
static const double MT_D2R = MT_PI / 180.0;   // degrees -> radians
static const double MT_R2D = 180.0 / MT_PI;   // radians -> degrees

// ============================================================
// Inline helpers
// ============================================================

// Clamp a value to [lo, hi]
template <typename T>
inline T mt_clamp(T val, T lo, T hi) {
    return (val < lo) ? lo : (val > hi) ? hi : val;
}

// Convert degrees to radians
inline double d2r(double deg) { return deg * MT_D2R; }

// Convert radians to degrees
inline double r2d(double rad) { return rad * MT_R2D; }

// Safe arctan2 returning [0, 2*pi)
inline double atan2_pos(double y, double x) {
    double a = std::atan2(y, x);
    return (a < 0.0) ? a + 2.0 * MT_PI : a;
}

// ============================================================
// Grid helpers
// ============================================================

// Linear index into a row-major matrix (0-based)
inline int idx(int row, int col, int ncols) {
    return row * ncols + col;
}

// Check whether (row, col) is inside the grid
inline bool in_grid(int row, int col, int nrows, int ncols) {
    return row >= 0 && row < nrows && col >= 0 && col < ncols;
}

// ============================================================
// Solar position helpers
// Ported from microclima/NicheMapR — replicate exactly.
// ============================================================

struct solmodel {
    double zend;   // solar zenith angle (degrees)
    double zenr;   // solar zenith angle (radians)
    double azid;   // solar azimuth (degrees, N = 0 clockwise)
    double azir;   // solar azimuth (radians)
};

// Astronomical Julian day
inline int juldayCpp(int year, int month, int day) {
    double dd   = day + 0.5;
    int    madj = month + (month < 3) * 12;
    int    yadj = year  + (month < 3) * -1;
    double j    = std::trunc(365.25  * (yadj + 4716))
                + std::trunc(30.6001 * (madj + 1))
                + dd - 1524.5;
    int b  = 2 - static_cast<int>(std::trunc(yadj / 100.0))
               + static_cast<int>(std::trunc(std::trunc(yadj / 100.0) / 4.0));
    int jd = static_cast<int>(j + (j > 2299160.0) * b);
    return jd;
}

// Apparent solar time (decimal hours)
// lt  = UTC decimal time;  lon = longitude (degrees, +E)
inline double soltimeCpp(int jd, double lt, double lon) {
    double m   = 6.24004077 + 0.01720197 * (jd - 2451545.0);
    double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2.0 * m + 3.5932);
    return lt + (4.0 * lon + eot) / 60.0;
}

// Solar zenith and azimuth at a given location and UTC clock time
inline solmodel solpositionCpp(double lat, double lon,
                                int year, int month, int day,
                                double lt) {
    int    jd   = juldayCpp(year, month, day);
    double st   = soltimeCpp(jd, lt, lon);
    double latr = lat * MT_D2R;
    double tt   = 0.261799 * (st - 12.0);
    double dec  = (MT_PI * 23.5 / 180.0) *
                  std::cos(2.0 * MT_PI * ((jd - 159.5) / 365.25));

    // cos(zenith)
    double coh = std::sin(dec) * std::sin(latr)
               + std::cos(dec) * std::cos(latr) * std::cos(tt);
    coh = mt_clamp(coh, -1.0, 1.0);
    double z = std::acos(coh) * MT_R2D;

    // Solar azimuth
    // sin(zenith) = cos(solar_elevation) — guard against zenith = 0
    double sin_zen = std::sqrt(std::max(1.0 - coh * coh, 0.0));
    if (sin_zen < 1e-6) sin_zen = 1e-6;

    double sazi = mt_clamp(std::cos(dec) * std::sin(tt) / sin_zen, -1.0, 1.0);
    double num  = std::cos(dec) * std::sin(tt);
    double dnom = std::sin(latr) * std::cos(dec) * std::cos(tt)
                - std::cos(latr) * std::sin(dec);
    double cazi_mag = std::sqrt(num * num + dnom * dnom);
    double cazi = (cazi_mag > 1e-9) ? dnom / cazi_mag : 1.0;

    double sqt = std::max(1.0 - sazi * sazi, 0.0);
    double azi = 180.0 + (180.0 / MT_PI) * std::atan(sazi / std::sqrt(sqt + 1e-18));
    if (cazi < 0.0) {
        azi = (sazi < 0.0) ? 180.0 - azi : 540.0 - azi;
    }

    solmodel s;
    s.zend = z;
    s.zenr = z * MT_D2R;
    s.azid = azi;
    s.azir = azi * MT_D2R;
    return s;
}

// Per-cell direct-beam solar index: cos(incidence), clamped to [0, 1].
// Returns 0 when sun is below the flat horizon (zend >= 90).
// Matches the solarindexCpp formula from microclima/NicheMapR.
inline double solarindexCpp(double slope_deg, double aspect_deg,
                             double zend, double azid) {
    if (zend >= 90.0) return 0.0;
    double zr = zend    * MT_D2R;
    double sr = slope_deg  * MT_D2R;
    double si = (slope_deg == 0.0)
        ? std::cos(zr)
        : std::cos(zr) * std::cos(sr)
          + std::sin(zr) * std::sin(sr)
            * std::cos((azid - aspect_deg) * MT_D2R);
    return (si < 0.0) ? 0.0 : si;
}

// ============================================================
// Canopy extinction coefficient (shared by twostream.cpp and
// solar.cpp's veg_path_atten_cpp())
// ============================================================

// Direct-beam canopy extinction coefficient for a horizontal reference
// surface, K_be(zenith, x) = G(zenith) / cos(zenith), where G is the
// Ross-Nilson leaf-area projection function -- the Campbell & Norman
// (1998) ellipsoidal leaf angle distribution formulation.
//
//   zen_r - solar zenith angle (radians)
//   xv    - leaf angle distribution parameter: ratio of mean leaf-area
//           projection onto a horizontal plane to that onto a vertical
//           plane (1 = spherical/isotropic, 0 = vertical leaves,
//           Inf = horizontal leaves; Campbell & Norman 1998 convention)
//
// Returns a large sentinel value (6000) at grazing incidence (cos_z ~ 0)
// so exp(-K_be * L) attenuates fully to zero rather than producing Inf.
//
// This is the "k" used two ways elsewhere in the package: twostream.cpp's
// canopy_kd() additionally slope-corrects it via a target surface's
// incidence cosine (k * cos(zenith) / si); solar.cpp's
// veg_path_atten_cpp() uses it directly against true vertical path
// segments during 3D ray marching, with no slope correction needed since
// the ray geometry is already fully 3D.
inline double kBeCpp(double zen_r, double xv) {
    double cos_z = std::cos(zen_r);
    double sin_z = std::sin(zen_r);
    double tan_z = (cos_z > 1e-10) ? sin_z / cos_z : 1e10;

    if (xv == 1.0)          return (cos_z > 1e-10) ? 0.5 / cos_z : 6000.0;
    // Vertical leaves, azimuthally random: G(zenith) = (2/pi) * sin(zenith)
    // (mean of |cos(phi)| over a uniform leaf azimuth is 2/pi), so
    // K_be = G / cos(zenith) = (2/pi) * tan(zenith) -- NOT tan(zenith)
    // alone, which is the projection for a single leaf whose azimuth is
    // locked to face the sun exactly, not the azimuthally-averaged canopy
    // value assumed everywhere else here (G as a function of zenith only).
    if (xv == 0.0)          return (2.0 / MT_PI) * tan_z;
    if (!std::isfinite(xv)) return 1.0;
    return std::sqrt(xv * xv + tan_z * tan_z) /
           (xv + 1.774 * std::pow(xv + 1.182, -0.733));
}
