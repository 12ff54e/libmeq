#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

#define ZQ_TIMER_IMPLEMENTATION
#include "MagneticEquilibrium.h"
#undef ZQ_TIMER_IMPLEMENTATION
#include "test/Assertion.hpp"
#include "test/GFileTestData.h"

auto make_equilibrium(int cocos, int ip, int b, bool si) {
    std::istringstream stream(analytic_gfile("analytic", cocos, ip, b));
    GFileRawData<double> data;
    data.read(stream, Cocos(cocos));
    return MagneticEquilibrium(
        data, 8, 32, 8, .8, si,
        [](auto one) {
            return std::array{one.psi, one.f, one.q, one.p, one.b2j};
        },
        [](auto one, auto two) {
            const auto b2 =
                (one.f * one.f + two.dpdr * two.dpdr + two.dpdz * two.dpdz) /
                (two.r * two.r);
            return std::array{two.dpdr, two.dpdz, two.r, two.z, one.b2j / b2};
        });
}

void test_equivalent_conventions(Assertion& assertion) {
    auto close = [](double a, double b) {
        return std::isfinite(a) && std::isfinite(b) &&
               std::abs(a - b) < 2e-9 * std::max(1., std::abs(b));
    };
    for (const bool si : {true, false}) {
        for (const int ip : {-1, 1}) {
            for (const int b : {-1, 1}) {
                const auto reference = make_equilibrium(2, ip, b, si);
                const double length_unit = si ? 1. : 1.5;
                const double field_unit = si ? 1. : b * 2.;
                const double current_unit = length_unit * field_unit;
                const double flux_unit = length_unit * current_unit;
                const double pressure_unit =
                    si ? 1. : 4 / (4e-7 * std::numbers::pi);
                assertion(close(reference.psi_delta(),
                                ip * .8 * .045 / (7 * flux_unit)),
                          "Flux spacing retains its COCOS 2 sign and units.");
                assertion(
                    close(reference.axis_value_1d(1), b * 3 / current_unit) &&
                        close(reference.axis_value_1d(2),
                              ip * b * 2 / std::sqrt(.92)) &&
                        close(reference.axis_value_1d(4),
                              ip * 6 / (std::sqrt(.92) * current_unit)),
                    "Axis F, q and B squared J use signed COCOS 2 values.");
                assertion(
                    close(reference.axis_value_1d(3),
                          .04 * 4 / (4e-7 * std::numbers::pi) / pressure_unit),
                    "Axis pressure uses pascals in SI mode and mu0*p/B0 "
                    "squared otherwise.");
                for (std::size_t i = 0; i < 8; ++i) {
                    const auto u = ip * reference.data_1d(0)[i] * flux_unit;
                    assertion(
                        close(reference.data_1d(1)[i],
                              b * (3 + .2 * u) / current_unit) &&
                            close(reference.data_1d(2)[i],
                                  ip * b * (3 + .2 * u) /
                                      std::sqrt(.92 *
                                                (2.25 - 2 * u * .8 / .92))) &&
                            close(reference.data_1d(3)[i],
                                  .04 * 4 / (4e-7 * std::numbers::pi) *
                                      (1 - .2 * u / .045) / pressure_unit),
                        "Profiles stay ordered from axis to edge.");
                    assertion(reference.data_1d(4)[i] * ip * field_unit > 0,
                              "The resolved Jacobian retains the flux sign.");
                    for (std::size_t j = 0; j <= 32; ++j) {
                        const auto x =
                            reference.data_2d(2)(i, j) * length_unit - 1.5;
                        const auto z = reference.data_2d(3)(i, j) * length_unit;
                        assertion(
                            close(reference.data_2d(0)(i, j),
                                  ip * (1.2 * x + .2 * z) / current_unit) &&
                                close(reference.data_2d(1)(i, j),
                                      ip * (.2 * x + .8 * z) / current_unit),
                            "Returned flux derivatives match the analytic "
                            "field.");
                    }
                }
                for (const int offset : {0, 10}) {
                    for (int base = 1; base <= 8; ++base) {
                        const auto candidate =
                            make_equilibrium(base + offset, ip, b, si);
                        bool equal =
                            close(candidate.psi_delta(), reference.psi_delta());
                        for (std::size_t k = 0; k < 5; ++k) {
                            equal = equal &&
                                    close(candidate.axis_value_1d(k),
                                          reference.axis_value_1d(k)) &&
                                    close(candidate.axis_value_2d(k),
                                          reference.axis_value_2d(k));
                            for (std::size_t i = 0; i < 8; ++i) {
                                equal = equal && close(candidate.data_1d(k)[i],
                                                       reference.data_1d(k)[i]);
                                for (std::size_t j = 0; j <= 32; ++j) {
                                    equal = equal &&
                                            close(candidate.data_2d(k)(i, j),
                                                  reference.data_2d(k)(i, j));
                                }
                            }
                        }
                        assertion(
                            equal,
                            "Equivalent COCOS " +
                                std::to_string(base + offset) +
                                " gives the same signed Boozer fields (Ip=" +
                                std::to_string(ip) +
                                ", B=" + std::to_string(b) +
                                ", SI=" + std::to_string(si) + ").");
                    }
                }
            }
        }
    }
}

int main() {
    Assertion assertion;
    test_equivalent_conventions(assertion);
    return assertion.status();
}
