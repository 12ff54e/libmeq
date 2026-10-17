#ifndef MEQ_COCOS_H
#define MEQ_COCOS_H

#include <numbers>
#include <stdexcept>

// Sauter and Medvedev, CPC 184 (2013), Table I and Appendix C.
// Transform coordinates without reversing the physical field or current.
class Cocos {
   public:
    explicit constexpr Cocos(int index = 1) : index_(index) {
        if (!((index >= 1 && index <= 8) || (index >= 11 && index <= 18))) {
            throw std::invalid_argument("COCOS index must be 1-8 or 11-18.");
        }
    }

    constexpr int index() const { return index_; }
    constexpr int sigma_bp() const {
        const auto base = index_ % 10;
        return base == 1 || base == 2 || base == 5 || base == 6 ? 1 : -1;
    }
    constexpr int sigma_r_phi_z() const { return index_ % 2 ? 1 : -1; }
    constexpr int sigma_r_theta_phi() const {
        const auto base = index_ % 10;
        return base == 1 || base == 2 || base == 7 || base == 8 ? 1 : -1;
    }

    template <typename T = double>
    constexpr T flux_scale() const {
        return index_ > 10 ? 2 * std::numbers::pi_v<T> : T{1};
    }

    template <typename T>
    struct Factors {
        T flux;        // psi, including its R and Z derivatives
        T derivative;  // dp/dpsi and F*dF/dpsi
        int toroidal;  // F, B_phi and I_p
        int safety_factor;
    };

    template <typename T = double>
    constexpr Factors<T> factors_to(Cocos target) const {
        const int toroidal = sigma_r_phi_z() * target.sigma_r_phi_z();
        const T flux = toroidal * sigma_bp() * target.sigma_bp() *
                       target.flux_scale<T>() / flux_scale<T>();
        return {flux, 1 / flux, toroidal,
                sigma_r_theta_phi() * target.sigma_r_theta_phi()};
    }

   private:
    int index_;
};

#endif  // MEQ_COCOS_H
