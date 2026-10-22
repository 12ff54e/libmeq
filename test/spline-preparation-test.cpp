#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "BSplineInterpolation/include/bsintp/Interpolation.hpp"
#include "SplinePreparation.h"
#include "test/Assertion.hpp"

void test_workload(Assertion& assertion) {
    std::vector<double> values(65), coordinates(65);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const auto x = double(i) / (values.size() - 1);
        coordinates[i] = x * x;
        values[i] = std::sin(coordinates[i]) + coordinates[i];
    }
    intp::InterpolationFunction1D<5> function(
        intp::util::get_range(coordinates), intp::util::get_range(values));
    assertion(!util::prepare_polynomial(function, 10) &&
                  !function.spline().has_polynomial(),
              "A short workload keeps control-point evaluation.");
    const auto direct = function;
    assertion(util::prepare_polynomial(function, 10000) &&
                  function.spline().has_polynomial(),
              "A repeated workload prepares polynomial coefficients.");
    assertion(util::prepare_polynomial(function, 0),
              "A later short workload retains prepared coefficients.");
    for (const auto x :
         {0., .01, .5, std::nextafter(.5, 0.), std::nextafter(.5, 1.), 1.}) {
        assertion(std::abs(function(x) - direct(x)) < 2e-12 &&
                      std::abs(function.derivative({x}, {1}) -
                               direct.derivative({x}, {1})) < 2e-10,
                  "Preparation preserves nonuniform values and derivatives.");
    }
    intp::InterpolationFunction1D<5> empty;
    assertion(!util::prepare_polynomial(
                  empty, std::numeric_limits<std::size_t>::max()),
              "An empty spline never requests polynomial conversion.");
}

void test_grid_cost(Assertion& assertion) {
    for (const std::size_t size : {17u, 129u}) {
        intp::Mesh<double, 2> data(size, size);
        for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j) {
                const auto x = double(i) / (size - 1);
                const auto y = double(j) / (size - 1);
                data(i, j) = std::sin(x) * std::cos(y) + .1 * x * y;
            }
        intp::InterpolationFunction<double, 2, 5> function(
            std::move(data), std::pair{0., 1.}, std::pair{0., 1.});
        const auto direct = function;
        const auto prepared = util::prepare_polynomial(function, 10000);
        assertion(prepared == (size == 17),
                  "The same workload converts a small grid but not a large "
                  "one.");
        if (!prepared) continue;
        const auto knot = function.spline().knots_begin(0)[8];
        for (const auto x :
             {0., std::nextafter(knot, 0.), knot, std::nextafter(knot, 1.), 1.})
            for (const auto y : {0., .37, 1.}) {
                assertion(std::abs(function(x, y) - direct(x, y)) < 2e-12,
                          "Polynomial values agree at knots and boundaries.");
                for (const auto derivative :
                     {std::array<std::size_t, 2>{1, 0}, {0, 1}})
                    assertion(
                        std::abs(function.derivative({x, y}, derivative) -
                                 direct.derivative({x, y}, derivative)) < 2e-10,
                        "Polynomial first derivatives agree at knot sides.");
            }
    }
}

int main() {
    Assertion assertion;
    test_workload(assertion);
    test_grid_cost(assertion);
    return assertion.status();
}
