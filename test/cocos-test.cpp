#include <cmath>
#include <sstream>

#include "Cocos.h"
#include "GFileRawData.h"
#include "test/Assertion.hpp"
#include "test/GFileTestData.h"

void test_factors(Assertion& assertion) {
    for (const int offset : {0, 10}) {
        for (int base = 1; base <= 8; ++base) {
            const Cocos source(base + offset);
            const auto [bp, phi, theta] = cocos_signs[base - 1];
            const auto scale = offset ? 2 * std::numbers::pi : 1.;
            const auto factors = source.factors_to(Cocos(2));
            assertion(
                factors.flux == -phi * bp / scale &&
                    factors.derivative == -phi * bp * scale &&
                    factors.toroidal == -phi && factors.safety_factor == theta,
                "Conversion to COCOS 2 matches the published sign table.");
            // Physical R/Z components use sigma_RphiZ*sigma_Bp/scale.
            const double dpdr = .8, dpdz = -.3, r = 1.7;
            assertion(std::abs(phi * bp * dpdz / (scale * r) +
                               factors.flux * dpdz / r) < 1e-15 &&
                          std::abs(-phi * bp * dpdr / (scale * r) -
                                   factors.flux * dpdr / r) < 1e-15,
                      "Converting flux preserves the physical poloidal field.");
        }
    }
    for (const int invalid : {-12, 0, 9, 10, 19}) {
        bool rejected = false;
        try {
            const Cocos cocos(invalid);
        } catch (const std::invalid_argument&) { rejected = true; }
        assertion(rejected, "Invalid COCOS indices are rejected.");
    }
}

void test_metadata(Assertion& assertion) {
    for (const auto& header : {"CHEASE", "chease", "Made by ChEaSe 2026"}) {
        std::istringstream stream(analytic_gfile(header));
        GFileRawData<double> data;
        stream >> data;
        assertion(data.is_complete() && data.cocos.index() == 2,
                  "CHEASE detection ignores case and position in the header.");
        assertion(std::abs(data.flux_magnetic_axis - .37) < 1e-15 &&
                      data.f_pol.front() == 3,
                  "The reader retains native flux and profile values.");
    }
    for (const auto& header : {"EFIT", "efit", "unknown", ""}) {
        for (const int ip : {-1, 1}) {
            for (const int b : {-1, 1}) {
                // Native EFIT keeps flux increasing and q positive.
                const int cocos = ip > 0 ? (b > 0 ? 1 : 5) : (b > 0 ? 3 : 7);
                std::istringstream stream(
                    analytic_gfile(header, cocos, -ip, -b));
                GFileRawData<double> data;
                stream >> data;
                assertion(
                    data.is_complete() && data.cocos.index() == cocos,
                    "EFIT and unknown producers infer all four sign cases.");
            }
        }
    }
    // Other EFIT writers retain signed psi/q in fixed COCOS 1.
    for (const int ip : {-1, 1}) {
        for (const int b : {-1, 1}) {
            std::istringstream stream(analytic_gfile("EFIT", 1, ip, b));
            GFileRawData<double> data;
            stream >> data;
            assertion(data.cocos.index() == 1,
                      "Signed EFIT flux and q remain COCOS 1.");
        }
    }
    for (const auto& tag : {"COCOS=11", "cocos_out:11", "cocos 11"}) {
        std::istringstream stream(analytic_gfile("CHEASE", 11, 1, 1, 17, tag));
        GFileRawData<double> data;
        stream >> data;
        assertion(data.cocos.index() == 11,
                  "An explicit COCOS tag overrides the producer hint.");
    }
    for (const auto& tag : {"COCOS=0", "COCOS=19", "COCOS=2x", "COCOS=2.5",
                            "COCOS=", "COCOS=no"}) {
        std::istringstream stream(analytic_gfile(tag));
        GFileRawData<double> data;
        bool rejected = false;
        try {
            stream >> data;
        } catch (const std::invalid_argument&) { rejected = true; }
        assertion(rejected && !data.is_complete(),
                  "A malformed explicit COCOS tag cannot become a valid read.");
        stream.clear();
        stream.seekg(0);
        data.read(stream, Cocos(2));
        assertion(data.is_complete() && data.cocos.index() == 2,
                  "The API override takes precedence over header metadata.");
    }
    GFileRawData<double> data;
    std::istringstream first(analytic_gfile("EFIT", 1, 1, 1, 17, " COCOS=1"));
    first >> data;
    std::istringstream second(analytic_gfile(std::string(42, 'x') + "CHEASE"));
    second >> data;
    assertion(
        data.cocos.index() == 2 && data.extra_metadata.empty() &&
            data.f_pol.size() == data.nw && data.boundary.size() == 128 &&
            data.geometric_poloidal_angles.size() == 128,
        "Repeated reads replace all data and inspect all 48 header bytes.");
    std::istringstream broken("short\n");
    broken >> data;
    assertion(!data.is_complete(), "A failed reread clears completion state.");
}

int main() {
    Assertion assertion;
    test_factors(assertion);
    test_metadata(assertion);
    return assertion.status();
}
