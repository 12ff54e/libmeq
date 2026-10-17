#ifndef MEQ_GFILE_TEST_DATA_H
#define MEQ_GFILE_TEST_DATA_H

#include <array>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

// Independent COCOS table: sigma_Bp, sigma_RphiZ, sigma_rhothetaphi.
inline constexpr std::array<std::array<int, 3>, 8> cocos_signs{{{1, 1, 1},
                                                                {1, -1, 1},
                                                                {-1, 1, -1},
                                                                {-1, -1, -1},
                                                                {1, 1, -1},
                                                                {1, -1, -1},
                                                                {-1, 1, 1},
                                                                {-1, -1, 1}}};

// A tilted quadratic equilibrium expressed in any COCOS convention. Signs
// refer to COCOS 2; the writer does not use the production transformation.
inline std::string analytic_gfile(const std::string& header = "CHEASE fixture",
                                  int cocos = 2,
                                  int current_sign = 1,
                                  int field_sign = 1,
                                  unsigned n = 17,
                                  const std::string& extra = "") {
    const auto [bp, phi, theta] = cocos_signs[(cocos % 10) - 1];
    const double flux_scale =
        -phi * bp * (cocos > 10 ? 2 * std::numbers::pi : 1);
    constexpr double axis_flux = .37, edge_flux = .045, det = .92;
    const double pressure0 = .04 * 4 / (4e-7 * std::numbers::pi);
    std::ostringstream out;
    out << std::left << std::setw(48) << header.substr(0, 48) << std::right
        << std::setw(4) << 0 << std::setw(4) << n << std::setw(4) << n << extra
        << '\n'
        << std::scientific << std::setprecision(16);
    auto values = [&](const std::vector<double>& data) {
        for (std::size_t i = 0; i < data.size(); ++i) {
            out << ' ' << data[i];
            if (i % 5 == 4 || i + 1 == data.size()) { out << '\n'; }
        }
    };
    values({1.4, 1.4, 1.5, .8, 0});
    values({1.5, 0, flux_scale * current_sign * axis_flux,
            flux_scale * current_sign * (axis_flux + edge_flux),
            -phi * field_sign * 2.});
    values({-phi * current_sign * 1e6, 0, 0, 0, 0});
    values({0, 0, flux_scale * current_sign * (axis_flux + edge_flux), 0, 0});
    std::vector<double> f(n), pressure(n), ff_prime(n), p_prime(n), q(n);
    for (unsigned i = 0; i < n; ++i) {
        const auto u = edge_flux * i / (n - 1);
        f[i] = -phi * field_sign * (3 + .2 * u);
        pressure[i] = pressure0 * (1 - .2 * u / edge_flux);
        ff_prime[i] = (3 + .2 * u) * .2 * current_sign / flux_scale;
        p_prime[i] = -pressure0 * .2 / edge_flux * current_sign / flux_scale;
        q[i] = theta * field_sign * current_sign * (3 + .2 * u) /
               std::sqrt(det * (2.25 - 2 * u * .8 / det));
    }
    values(f);
    values(pressure);
    values(ff_prime);
    values(p_prime);
    std::vector<double> flux;
    for (unsigned j = 0; j < n; ++j) {
        for (unsigned i = 0; i < n; ++i) {
            const auto x = -.7 + 1.4 * i / (n - 1);
            const auto z = -.7 + 1.4 * j / (n - 1);
            flux.push_back(
                flux_scale * current_sign *
                (axis_flux + .5 * (1.2 * x * x + .4 * x * z + .8 * z * z)));
        }
    }
    values(flux);
    values(q);
    constexpr unsigned boundary_size = 128;
    out << boundary_size << " 0\n";
    std::vector<double> boundary;
    for (unsigned i = 0; i < boundary_size; ++i) {
        const auto theta_geo = 2 * std::numbers::pi * i / boundary_size;
        const auto c = std::cos(theta_geo), s = std::sin(theta_geo);
        const auto radius =
            std::sqrt(2 * edge_flux / (1.2 * c * c + .4 * c * s + .8 * s * s));
        boundary.push_back(1.5 + radius * c);
        boundary.push_back(radius * s);
    }
    values(boundary);
    return out.str();
}

#endif  // MEQ_GFILE_TEST_DATA_H
