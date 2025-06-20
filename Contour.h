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
    val_type flux_;
    std::vector<pt_type> pts_;

   public:
    Contour() = default;

    template <typename F>
    Contour(val_type psi, const F& flux, const auto& g_file) : flux_(psi) {
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
};

#endif  // MEQ_CONTOUR_H
