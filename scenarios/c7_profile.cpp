/*!
 * @file scenarios/c7_profile.cpp
 * @brief The tabulated Model C7 atmosphere (Avrett & Loeser 2008, Table 26) and
 *        its two interpolants.
 * @ingroup scenarios
 *
 * The release `model_column` IC samples its temperature and ionization profile
 * from here; the legacy model_c7 scenario family reads MODEL_C7 directly.
 * Moved verbatim out of the former scenarios/model_c7.cpp.
 */
#include "c7_profile.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace chromosphere {

// Model C7 chromosphere + lower transition region profile.
// Columns: height [km], n_e [m^-3], n_n (= n_HI) [m^-3], T [K].
//
// Rows 0–14 are the original chromospheric table (writeup Table 1); rows 15–31
// extend it up through the lower TR to h ≈ 2153 km, T ≈ 23000 K — just above the
// κ_e/κ_n conductivity crossover (2149 km) — so the upper boundary can sit at a
// physically-grounded TR base (RTV T_0 ≈ 2×10⁴ K) rather than mid-chromosphere.
// All rows are Avrett & Loeser (2008) Table 26 (Model C7, ApJS 175, 229, p. 248):
// n_e ← their n_e column, n_n ← their n_HI column, T ← their T column, with the
// cm⁻³→m⁻³ (×10⁶) conversion applied. n_HI drops ~2 orders of magnitude across
// the thin TR (2138→2153 km) as hydrogen ionizes.
// Rows 0..N_C7_BASE-1 are the standard model_c7 column (chromosphere + lower TR,
// top h ≈ 2.153 Mm, T ≈ 23 kK). Rows N_C7_BASE..N_C7-1 extend it through the real
// Avrett & Loeser (2008) Table 26 rows, up the transition region and into the
// CORONA — all the way to h ≈ 68 Mm, T ≈ 1.59 MK (the top of Table 26). These are
// used only when model_c7_ic(extended=true) (the corona-as-boundary gentle-
// evaporation configuration), which truncates the table at grid.corona_top_km so
// the upper boundary sits high in a resolved coronal VOLUME (default ≈ 10 Mm,
// T ≈ 0.86 MK). h increases / T increases monotonically down the column, so the
// cubic spline stays well-behaved.
const Real MODEL_C7[N_C7][4] = {
    {1.003e+03f, 1.903e+17f, 2.693e+19f, 6.225e+03f},
    {1.032e+03f, 2.021e+17f, 2.179e+19f, 6.315e+03f},
    {1.065e+03f, 2.091e+17f, 1.722e+19f, 6.400e+03f},
    {1.101e+03f, 2.104e+17f, 1.340e+19f, 6.474e+03f},
    {1.143e+03f, 1.995e+17f, 1.009e+19f, 6.531e+03f},
    {1.214e+03f, 1.760e+17f, 6.322e+18f, 6.576e+03f},
    {1.299e+03f, 1.489e+17f, 3.641e+18f, 6.598e+03f},
    {1.398e+03f, 1.438e+17f, 1.880e+18f, 6.610e+03f},
    {1.520e+03f, 1.423e+17f, 7.956e+17f, 6.623e+03f},
    {1.617e+03f, 1.267e+17f, 4.012e+17f, 6.633e+03f},
    {1.722e+03f, 1.027e+17f, 1.916e+17f, 6.643e+03f},
    {1.820e+03f, 7.969e+16f, 9.843e+16f, 6.652e+03f},
    {1.894e+03f, 6.450e+16f, 6.005e+16f, 6.660e+03f},
    {1.946e+03f, 5.526e+16f, 4.371e+16f, 6.667e+03f},
    {1.989e+03f, 4.826e+16f, 3.447e+16f, 6.674e+03f},
    // --- lower transition region (Avrett & Loeser 2008 Table 26, rows 87→68) ---
    {2.024e+03f, 4.307e+16f, 2.916e+16f, 6.680e+03f},
    {2.055e+03f, 3.883e+16f, 2.583e+16f, 6.686e+03f},
    {2.083e+03f, 3.498e+16f, 2.351e+16f, 6.694e+03f},
    {2.098e+03f, 3.293e+16f, 2.258e+16f, 6.700e+03f},
    {2.110e+03f, 3.128e+16f, 2.204e+16f, 6.706e+03f},
    {2.120e+03f, 2.985e+16f, 2.176e+16f, 6.718e+03f},
    {2.130e+03f, 2.826e+16f, 2.157e+16f, 6.768e+03f},
    {2.134e+03f, 2.741e+16f, 2.118e+16f, 6.870e+03f},
    {2.138e+03f, 2.608e+16f, 1.954e+16f, 7.248e+03f},
    {2.142e+03f, 2.432e+16f, 1.693e+16f, 7.950e+03f},
    {2.145e+03f, 2.325e+16f, 1.165e+16f, 9.115e+03f},
    {2.147e+03f, 2.274e+16f, 4.413e+15f, 1.098e+04f},
    {2.149e+03f, 2.046e+16f, 1.749e+15f, 1.320e+04f},
    {2.150e+03f, 1.761e+16f, 8.787e+14f, 1.576e+04f},
    {2.151e+03f, 1.550e+16f, 5.412e+14f, 1.814e+04f},
    {2.152e+03f, 1.383e+16f, 3.809e+14f, 2.051e+04f},
    {2.153e+03f, 1.237e+16f, 2.769e+14f, 2.310e+04f},
    // --- upper transition region (Avrett & Loeser 2008 Table 26, rows 67→20;
    // h ≈ 2.153→2.628 Mm, T ≈ 25 kK→0.29 MK). Used only by model_c7_ic(extended).
    {2.15348e+03f, 1.144e+16f, 2.210e+14f, 2.512e+04f},
    {2.15550e+03f, 9.013e+15f, 9.798e+13f, 3.226e+04f},
    {2.15757e+03f, 7.454e+15f, 4.068e+13f, 3.940e+04f},
    {2.15942e+03f, 6.565e+15f, 1.832e+13f, 4.514e+04f},
    {2.16126e+03f, 5.955e+15f, 8.397e+12f, 5.014e+04f},
    {2.16456e+03f, 5.222e+15f, 2.253e+12f, 5.779e+04f},
    {2.16999e+03f, 4.497e+15f, 4.054e+11f, 6.800e+04f},
    {2.17528e+03f, 4.049e+15f, 1.611e+11f, 7.633e+04f},
    {2.18608e+03f, 3.463e+15f, 7.053e+10f, 9.064e+04f},
    {2.20192e+03f, 2.953e+15f, 3.595e+10f, 1.078e+05f},
    {2.22004e+03f, 2.592e+15f, 2.165e+10f, 1.240e+05f},
    {2.24332e+03f, 2.288e+15f, 1.377e+10f, 1.416e+05f},
    {2.26969e+03f, 2.054e+15f, 9.482e+09f, 1.586e+05f},
    {2.29437e+03f, 1.894e+15f, 7.241e+09f, 1.726e+05f},
    {2.32578e+03f, 1.739e+15f, 5.502e+09f, 1.886e+05f},
    {2.37819e+03f, 1.554e+15f, 3.873e+09f, 2.120e+05f},
    {2.45790e+03f, 1.366e+15f, 2.632e+09f, 2.419e+05f},
    {2.53547e+03f, 1.240e+15f, 1.993e+09f, 2.667e+05f},
    {2.62832e+03f, 1.136e+15f, 1.546e+09f, 2.925e+05f},
    // --- corona (Avrett & Loeser 2008 Table 26, rows 19→1; h ≈ 2.69→68 Mm,
    // T ≈ 0.31→1.59 MK). The resolved coronal volume the evaporated chromospheric
    // mass fills. model_c7_ic(extended) truncates at grid.corona_top_km. n_e, n_HI
    // are the table's n_e, n_HI columns × 10⁶ (cm⁻³→m⁻³). ---
    {2.68810e+03f, 1.084e+15f, 1.353e+09f, 3.073e+05f},
    {2.76149e+03f, 1.030e+15f, 1.174e+09f, 3.240e+05f},
    {2.84440e+03f, 9.794e+14f, 1.022e+09f, 3.413e+05f},
    {2.97076e+03f, 9.165e+14f, 8.549e+08f, 3.650e+05f},
    {3.11379e+03f, 8.529e+14f, 7.082e+08f, 3.920e+05f},
    {3.42499e+03f, 7.695e+14f, 5.457e+08f, 4.329e+05f},
    {3.75271e+03f, 7.022e+14f, 4.361e+08f, 4.716e+05f},
    {4.29586e+03f, 6.242e+14f, 3.303e+08f, 5.243e+05f},
    {4.96864e+03f, 5.578e+14f, 2.560e+08f, 5.774e+05f},
    {5.76283e+03f, 5.021e+14f, 2.034e+08f, 6.293e+05f},
    {7.36082e+03f, 4.274e+14f, 1.455e+08f, 7.117e+05f},
    {8.97424e+03f, 3.768e+14f, 1.133e+08f, 7.783e+05f},
    {1.15962e+04f, 3.206e+14f, 8.348e+07f, 8.650e+05f},
    {1.53920e+04f, 2.669e+14f, 6.013e+07f, 9.637e+05f},
    {2.11331e+04f, 2.147e+14f, 4.163e+07f, 1.080e+06f},
    {2.66766e+04f, 1.805e+14f, 3.152e+07f, 1.170e+06f},
    {3.60795e+04f, 1.411e+14f, 2.164e+07f, 1.294e+06f},
    {4.70093e+04f, 1.108e+14f, 1.522e+07f, 1.410e+06f},
    {6.80844e+04f, 7.491e+13f, 8.860e+06f, 1.586e+06f}
};

// Photosphere → low-chromosphere rows (Avrett & Loeser 2008 Table 26, h = −10 …
// 971 km) that MODEL_C7 (which begins at h = 1003 km) does not carry. Same columns
// as MODEL_C7: height [km], n_e [m^-3], n_n (= n_HI) [m^-3], T [K]. Used by
// c7_full_profile to sample the COMPLETE atmosphere (photosphere → corona) — e.g.
// the single-fluid model_column corona testbed needs a realistic ionization
// fraction n_e/(n_e+n_HI) all the way down so the gas is neutral in the chromosphere
// and ionized in the corona (κ_e Spitzer). n_e has a minimum (~4.8×10^16 m^-3) near
// the temperature-minimum region (~750 km).
static const int N_C7_PHOTO = 29;
static const Real MODEL_C7_PHOTO[N_C7_PHOTO][4] = {
    {-1.000e+01f, 1.206e+20f, 1.222e+23f, 6.780e+03f},
    { 0.000e+00f, 8.397e+19f, 1.187e+23f, 6.583e+03f},
    { 1.000e+01f, 5.893e+19f, 1.148e+23f, 6.397e+03f},
    { 2.000e+01f, 4.258e+19f, 1.106e+23f, 6.231e+03f},
    { 3.500e+01f, 2.729e+19f, 1.040e+23f, 6.006e+03f},
    { 5.000e+01f, 1.909e+19f, 9.687e+22f, 5.826e+03f},
    { 7.500e+01f, 1.233e+19f, 8.434e+22f, 5.607e+03f},
    { 1.000e+02f, 8.748e+18f, 7.244e+22f, 5.431e+03f},
    { 1.250e+02f, 6.596e+18f, 6.150e+22f, 5.288e+03f},
    { 1.500e+02f, 5.132e+18f, 5.173e+22f, 5.165e+03f},
    { 1.750e+02f, 4.102e+18f, 4.299e+22f, 5.080e+03f},
    { 2.000e+02f, 3.308e+18f, 3.549e+22f, 5.010e+03f},
    { 2.500e+02f, 2.182e+18f, 2.379e+22f, 4.907e+03f},
    { 3.000e+02f, 1.433e+18f, 1.577e+22f, 4.805e+03f},
    { 3.500e+02f, 9.337e+17f, 1.034e+22f, 4.700e+03f},
    { 4.000e+02f, 6.034e+17f, 6.726e+21f, 4.590e+03f},
    { 4.500e+02f, 3.871e+17f, 4.322e+21f, 4.485e+03f},
    { 4.900e+02f, 2.718e+17f, 2.993e+21f, 4.435e+03f},
    { 5.250e+02f, 1.996e+17f, 2.155e+21f, 4.410e+03f},
    { 5.600e+02f, 1.467e+17f, 1.545e+21f, 4.400e+03f},
    { 6.150e+02f, 9.163e+16f, 9.058e+20f, 4.435e+03f},
    { 6.600e+02f, 6.424e+16f, 5.825e+20f, 4.510e+03f},
    { 7.000e+02f, 5.159e+16f, 3.914e+20f, 4.640e+03f},
    { 7.500e+02f, 4.762e+16f, 2.404e+20f, 4.840e+03f},
    { 8.000e+02f, 5.458e+16f, 1.495e+20f, 5.090e+03f},
    { 8.540e+02f, 8.055e+16f, 9.085e+19f, 5.430e+03f},
    { 9.000e+02f, 1.162e+17f, 6.086e+19f, 5.720e+03f},
    { 9.460e+02f, 1.513e+17f, 4.178e+19f, 5.969e+03f},
    { 9.710e+02f, 1.727e+17f, 3.428e+19f, 6.100e+03f},
};

// Sample the COMPLETE Model C7 atmosphere (photosphere → corona, h ≈ −10 km … 68 Mm)
// at height h_km: T linearly, n_e and n_HI log-linearly (they span many decades).
// Concatenates MODEL_C7_PHOTO (h < 1003 km) with MODEL_C7 (h ≥ 1003 km); clamps to
// the table ends. (model_c7 itself only ever uses MODEL_C7 from h = 1003 km up, so
// this helper does not change any existing scenario.)
void c7_full_profile(Real h_km, Real& T, Real& n_e, Real& n_HI) {
    auto row = [](int k) -> const Real* {
        return (k < N_C7_PHOTO) ? MODEL_C7_PHOTO[k] : MODEL_C7[k - N_C7_PHOTO];
    };
    const int N = N_C7_PHOTO + N_C7;
    auto unpack = [&](int k) { const Real* r = row(k); T = r[3]; n_e = r[1]; n_HI = r[2]; };
    if (h_km <= row(0)[0])     { unpack(0);     return; }
    if (h_km >= row(N - 1)[0]) { unpack(N - 1); return; }
    int i = 0;
    while (i < N - 1 && row(i + 1)[0] < h_km) ++i;
    const Real* a = row(i);
    const Real* b = row(i + 1);
    const Real t = (h_km - a[0]) / (b[0] - a[0]);
    T    = a[3] + t * (b[3] - a[3]);
    n_e  = std::exp(std::log(a[1]) + t * (std::log(b[1]) - std::log(a[1])));
    n_HI = std::exp(std::log(a[2]) + t * (std::log(b[2]) - std::log(a[2])));
}

double c7_full_temperature_pchip(double h_km) {
    // Fritsch-Carlson/Fritsch-Butland monotone cubic Hermite slopes, including
    // the standard one-sided endpoint limiter. The complete C7 temperature is
    // not globally monotone (it has a photospheric minimum), so slopes are
    // limited interval-by-interval and vanish at sign-changing extrema.
    struct TemperaturePchip {
        std::vector<double> h;
        std::vector<double> temperature;
        std::vector<double> slope;

        TemperaturePchip() {
            auto row = [](int k) -> const Real* {
                return (k < N_C7_PHOTO) ? MODEL_C7_PHOTO[k]
                                        : MODEL_C7[k-N_C7_PHOTO];
            };
            const int n = N_C7_PHOTO+N_C7;
            h.resize(n); temperature.resize(n); slope.assign(n, 0.0);
            for (int i = 0; i < n; ++i) {
                h[i] = row(i)[0];
                temperature[i] = row(i)[3];
            }
            std::vector<double> width(n-1), secant(n-1);
            for (int i = 0; i+1 < n; ++i) {
                width[i] = h[i+1]-h[i];
                secant[i] = (temperature[i+1]-temperature[i])/width[i];
            }
            for (int i = 1; i+1 < n; ++i) {
                if (secant[i-1]*secant[i] <= 0.0) continue;
                const double w1 = 2.0*width[i]+width[i-1];
                const double w2 = width[i]+2.0*width[i-1];
                slope[i] = (w1+w2)/(w1/secant[i-1]+w2/secant[i]);
            }
            auto endpoint = [](double h0, double h1, double d0, double d1) {
                double value = ((2.0*h0+h1)*d0-h0*d1)/(h0+h1);
                if (value*d0 <= 0.0) return 0.0;
                if (d0*d1 < 0.0 && std::abs(value) > 3.0*std::abs(d0))
                    return 3.0*d0;
                return value;
            };
            slope.front() = endpoint(width[0], width[1], secant[0], secant[1]);
            slope.back() = endpoint(width[n-2], width[n-3],
                                    secant[n-2], secant[n-3]);
        }
    };
    static const TemperaturePchip data;
    if (h_km <= data.h.front()) return data.temperature.front();
    if (h_km >= data.h.back()) return data.temperature.back();
    const std::size_t i = static_cast<std::size_t>(
        std::upper_bound(data.h.begin(), data.h.end(), h_km)-data.h.begin()-1);
    const double width = data.h[i+1]-data.h[i];
    const double s = (h_km-data.h[i])/width;
    const double s2 = s*s, s3 = s2*s;
    return (2.0*s3-3.0*s2+1.0)*data.temperature[i]
         + (s3-2.0*s2+s)*width*data.slope[i]
         + (-2.0*s3+3.0*s2)*data.temperature[i+1]
         + (s3-s2)*width*data.slope[i+1];
}

} // namespace chromosphere
