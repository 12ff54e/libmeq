#ifndef MEQ_CONTOUR_H
#define MEQ_CONTOUR_H

#include <algorithm>
#include <vector>

#include "GFileRawData.h"
#include "Vec.h"
#include "util.h"

template <typename T>
class Contour {
   public:
    using val_type = T;
    using pt_type = Vec<2, val_type>;

   private:
    const val_type flux_;
    std::vector<pt_type> pts_;
    const GFileRawData<val_type>& g_file_;

   public:
    template <typename F>
    Contour(val_type psi, const F& flux, const auto& g_file)
        : flux_(psi), g_file_(g_file) {
        pts_.reserve(g_file.boundary.size());
        for (size_t i = 0; i < g_file.boundary.size(); ++i) {
            pts_.emplace_back(util::vec_field_find_root(
                flux, g_file.magnetic_axis, g_file.boundary[i], psi));
        }
#ifdef _DEBUG
        size_t count{};
        for (size_t i = 0; i < g_file.boundary.size(); ++i) {
            auto pt0 = pts_[i];
            auto pt1 = pts_[(i + 1) % pts_.size()];
            if (util::abs((flux(.5 * (pt0 + pt1)) - psi) / psi) > 1.e-3) {
                ++count;
            }
        }
        if (count > 0) {
            std::cout << "Contour \\psi = " << psi << " has " << count << "/"
                      << pts_.size() << " not-so-accurate segments\n";
        }
#endif
    }
    // properties

    size_t size() const noexcept { return pts_.size(); }

    val_type flux() const noexcept { return flux_; }

    // element access

    const Vec<2, val_type>& operator[](std::size_t i) const { return pts_[i]; }

    // numerics method

    template <typename Field, typename Measure>
    val_type definite_integrate_along

        (Field f, Measure s) {
        // sort contour pts according to measure
        std::sort(pts_.begin(), pts_.end(),
                  [&](const pt_type& p1, const pt_type& p2) {
                      return s(p1) < s(p2);
                  });

        std::vector<val_type> abscissa;
        std::vector<val_type> ordinate;
        abscissa.reserve(size());
        ordinate.reserve(size());
        for (size_t i = 0; i < size(); ++i) {
            auto& pt = operator[](i);
            abscissa.emplace_back(s(pt.x(), pt.y()));
            ordinate.emplace_back(f(pt.x(), pt.y()));
        }

        // last point being identical to the first one required for periodic
        // interpolation

        abscissa.emplace_back(abscissa.front() + 2 * M_PI);
        ordinate.emplace_back(ordinate.front());

        intp::InterpolationFunction<val_type, 1, 3> f_interp(
            true, std::make_pair(abscissa.begin(), abscissa.end()),
            std::make_pair(ordinate.begin(), ordinate.end()));

        return util::integrate(f_interp, 0., 2 * M_PI);
    }

    template <typename Field>
    intp::InterpolationFunction1D<3, val_type> indefinite_integrate_along(
        Field f,
        bool normalization) {
        std::vector<val_type> ordinate;
        ordinate.reserve(size() + 1);

        for (auto& pt : pts_) { ordinate.emplace_back(f(pt)); }

        // last point being identical to the first one required for periodic
        // interpolation

        ordinate.push_back(ordinate.front());
        std::vector<val_type> abscissa;
        if (std::fpclassify(g_file_.geometric_poloidal_angles.front()) ==
            FP_ZERO) {
            abscissa.reserve(size() + 1);
        } else {
            abscissa.reserve(size() + 2);
            abscissa.push_back(0);
        }
        abscissa.insert(abscissa.end(),
                        g_file_.geometric_poloidal_angles.begin(),
                        g_file_.geometric_poloidal_angles.end());
        abscissa.push_back(g_file_.geometric_poloidal_angles.front() +
                           2 * M_PI);

        intp::InterpolationFunction1D<3, val_type> f_interp(
            std::make_pair(std::next(abscissa.begin()), abscissa.end()),
            std::make_pair(ordinate.begin(), ordinate.end()), true);

        // do the integration segment by segment

        abscissa.back() = 2 * M_PI;
        std::vector<val_type> integral;
        integral.reserve(abscissa.size());
        integral.push_back(0);
        for (size_t i = 0; i < abscissa.size() - 1; ++i) {
            integral.push_back(
                util::integrate_coarse(f_interp, abscissa[i], abscissa[i + 1]) +
                integral.back());
        }

        // normalization

        if (normalization) {
            const val_type coef = 2 * M_PI / integral.back();
            for (auto& v : integral) { v *= coef; }
        }

        intp::InterpolationFunction1D<3, val_type> integral_interp(
            std::make_pair(abscissa.begin(), abscissa.end()),
            std::make_pair(integral.begin(), integral.end()), true);

        return integral_interp;
    }
};

#endif  // MEQ_CONTOUR_H
