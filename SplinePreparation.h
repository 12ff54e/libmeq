#ifndef MEQ_SPLINE_PREPARATION_H
#define MEQ_SPLINE_PREPARATION_H

#include <cstddef>

namespace util {

// Build coefficients only when enough evaluations are expected to offset
// setup time. Require twice the degree-plus-one calls per active cell;
// cell-layout quintics needed roughly 6–9 calls/cell in measurements of
// values and first derivatives.
// Call before publishing the spline to other threads. Existing coefficients
// remain useful even when the next workload is small.
template <typename F>
bool prepare_polynomial(F& function, std::size_t evaluations) {
    static_assert(F::dim > 0 && F::dim <= 2);
    const auto& spline = function.spline();
    if (spline.has_polynomial()) { return true; }
    auto budget = evaluations / (2 * (F::order + 1));
    for (std::size_t dim = 0; dim < F::dim; ++dim) {
        const auto knots = spline.knots_num(dim);
        if (knots <= 2 * F::order + 1) { return false; }
        const auto cells = knots - 2 * F::order - 1;
        // Divide the budget instead of multiplying grid dimensions, so the
        // decision also works for grids whose cell count would overflow.
        if (budget < cells) { return false; }
        budget /= cells;
    }
    function.generate_polynomial_form();
    return true;
}

}  // namespace util

#endif  // MEQ_SPLINE_PREPARATION_H
