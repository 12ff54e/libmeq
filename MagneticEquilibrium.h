#ifndef MEQ_MAGNETIC_EQUILIBRIUM_H
#define MEQ_MAGNETIC_EQUILIBRIUM_H

#include "BSplineInterpolation/src/include/Interpolation.hpp"
#include "Contour.h"
#include "GFileRawData.h"
#include "Timer.h"

#ifdef MEQ_ZERNIKE_SERIES_
#include "Zernike.h"
#endif

template <typename T>
class MagneticEquilibrium {
   public:
    // interpolation order of internal use
    constexpr static std::size_t ORDER = 5;
    constexpr static std::size_t ORDER_OUT = 2;

    using val_type = T;

   protected:
    static constexpr std::size_t FIELD_NUM_2D = 4;
    static constexpr std::size_t FIELD_NUM_1D = 6;

    // The contour mesh grid of data stored here, and the resolution should be
    // higher than the output spec
    struct MagneticEquilibriumRaw_ {
        // - magnetic_field
        // - r
        // - z
        // - jacobian
        std::array<intp::Mesh<val_type, 2>, FIELD_NUM_2D> data_2d;
        // - safety_factor
        // - poloidal_current
        // - toroidal_current
        // - pressure
        // - minor_radius
        // - toroidal_flux
        std::array<std::vector<val_type>, FIELD_NUM_1D> data_1d;

        std::array<val_type, FIELD_NUM_2D> axis_value_2d;
        std::array<val_type, FIELD_NUM_1D> axis_value_1d;

        val_type flux_unit;
    };

    struct MagneticEquilibriumIntp_ {
        std::vector<val_type> psi_sample_for_output;

        // - magnetic_field
        // - r
        // - z
        // - jacobian
#ifdef MEQ_ZERNIKE_SERIES_
        std::array<Zernike::Series<val_type>, FIELD_NUM_2D>
#else
        std::array<
            intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>,
            FIELD_NUM_2D>
#endif
            intp_2d;
        // - safety_factor
        // - poloidal_current
        // - toroidal_current
        // - pressure
        // - minor_radius
        // - toroidal_flux
        std::array<intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>,
                   FIELD_NUM_1D>
            intp_1d;

        template <typename F,
                  std::size_t... indices_2d,
                  std::size_t... indices_1d>
        MagneticEquilibriumIntp_(const MagneticEquilibriumRaw_& spdata_raw,
                                 const MagneticEquilibrium& spdata,
                                 const F& generate_psi_for_output,
                                 std::index_sequence<indices_2d...>,
                                 std::index_sequence<indices_1d...>)
            : psi_sample_for_output{generate_psi_for_output(spdata.psi_delta(),
                                                            spdata.lsp)},
              intp_2d{spdata.create_2d_spline_(spdata_raw.data_2d[indices_2d],
                                               psi_sample_for_output)...},
              intp_1d{spdata.create_1d_spline_(spdata_raw.data_1d[indices_1d],
                                               psi_sample_for_output)...} {}
    };

   private:
    MagneticEquilibriumRaw_ generate_boozer_coordinate_(
        const GFileRawData<val_type>& gfile_data,
        std::size_t radial_sample,
        val_type psi_ratio) {
        if (radial_sample % 2 != 0) {
            throw std::runtime_error(
                "[MagneticEquilibrium] Radial sample point must be even.");
        }

        auto& timer = Timer::get_timer();
        timer.start("Create Boozer grid");
        const val_type left_bd = gfile_data.r_left;
        const val_type right_bd = gfile_data.r_left + gfile_data.dim.x();
        const val_type top_bd = gfile_data.z_mid + .5 * gfile_data.dim.y();
        const val_type bottom_bd = gfile_data.z_mid - .5 * gfile_data.dim.y();
        intp::InterpolationFunction<val_type, 2u, ORDER, val_type>
            flux_function(gfile_data.flux, std::make_pair(left_bd, right_bd),
                          std::make_pair(bottom_bd, top_bd));

        // check poloidal flux value at magnetic axis
        const auto psi_ma_intp = flux_function(gfile_data.magnetic_axis);
        if (std::abs((psi_ma_intp - gfile_data.flux_magnetic_axis) /
                     (gfile_data.flux_LCFS - gfile_data.flux_magnetic_axis)) >
            1.e-4) {
            std::cout << "The poloidal flux of magnetic axis given in gfile "
                         "deviates from interpolated value.\n"
                      << "  \\psi_p in gfile: " << gfile_data.flux_magnetic_axis
                      << "\n  \\psi_p from interpolation: " << psi_ma_intp
                      << '\n';
        } else {
            std::cout << "The poloidal flux of magnetic axis is "
                      << gfile_data.flux_magnetic_axis << '\n';
        }

        val_type psi_boundary_min = 10. * std::pow(gfile_data.r_center, 2) *
                                    std::abs(gfile_data.b_center);
        for (const auto& pt : gfile_data.boundary) {
            psi_boundary_min = std::min(psi_boundary_min, flux_function(pt));
        }
        std::cout << "The poloidal flux of last closed flux surface is "
                  << gfile_data.flux_LCFS << '\n'
                  << "Minimum of interpolated value at boundary points is "
                  << psi_boundary_min << '\n';

        const val_type psi_bd = gfile_data.flux_LCFS - psi_ma_intp;
        val_type psi_wall = psi_ratio * psi_bd;
        if (psi_wall > psi_boundary_min - psi_ma_intp) {
            psi_wall = psi_boundary_min - psi_ma_intp;
            std::cout << "Interpolated flux value at boundary is too small, so "
                         "psi_wall is set to this value.\n";
        }
        psi_delta_ = psi_wall / static_cast<val_type>(lsp - 1);

        // contours are from \\Delta\\psi to LCFS
        std::vector<Contour<val_type>> contours;
        contours.reserve(radial_sample);
        for (std::size_t i = 0; i < radial_sample; ++i) {
            contours.emplace_back(
                util::lerp(psi_delta_, psi_wall,
                           static_cast<val_type>(i) /
                               static_cast<val_type>(radial_sample - 1)) +
                    psi_ma_intp,
                flux_function, gfile_data);
        }

        constexpr val_type magnetic_constant = 4.e-7 * M_PI;

        // safety factor, on shifted psi
        intp::InterpolationFunction1D<ORDER, val_type, val_type>
            safety_factor_intp{std::make_pair(0., psi_bd),
                               intp::util::get_range(gfile_data.safety_factor)};

        // poloidal current, on shifted psi
        intp::InterpolationFunction1D<ORDER, val_type, val_type>
            poloidal_current_intp{std::make_pair(0., psi_bd),
                                  intp::util::get_range(gfile_data.f_pol)};

        // pressure, on shifted psi
        intp::InterpolationFunction1D<ORDER, val_type, val_type> pressure_intp{
            std::make_pair(0., psi_bd),
            intp::util::get_range(gfile_data.pressure)};

        // this following function accepts shifted psi (0 at m.a.)

        auto b2j_field = [&](Vec<2, val_type> pt, val_type psi) {
            val_type dp_dr = flux_function.derivative(pt, {1, 0});
            val_type dp_dz = flux_function.derivative(pt, {0, 1});
            val_type r = pt.x();
            pt -= gfile_data.magnetic_axis;
            val_type r2 = pt.L2_norm_square_();
            val_type f = poloidal_current_intp(psi);

            return (f * f + dp_dr * dp_dr + dp_dz * dp_dz) * r2 /
                   (r * (dp_dr * pt.x() + dp_dz * pt.y()));
        };

        auto b_field = [&](Vec<2, val_type> pt, val_type psi) {
            val_type dp_dr = flux_function.derivative(pt, {1, 0});
            val_type dp_dz = flux_function.derivative(pt, {0, 1});
            val_type f = poloidal_current_intp(psi);

            return std::sqrt(f * f + dp_dr * dp_dr + dp_dz * dp_dz) / pt.x();
        };

        auto bp2j_field = [&](Vec<2, val_type> pt, val_type) {
            val_type dp_dr = flux_function.derivative(pt, {1, 0});
            val_type dp_dz = flux_function.derivative(pt, {0, 1});
            val_type r = pt.x();
            pt -= gfile_data.magnetic_axis;
            val_type r2 = pt.L2_norm_square_();

            return (dp_dr * dp_dr + dp_dz * dp_dz) * r2 /
                   ((dp_dr * pt.x() + dp_dz * pt.y()) * r);
        };

        constexpr val_type PI2 = 2 * M_PI;

        std::vector<val_type> poloidal_angles{
            gfile_data.geometric_poloidal_angles};
        poloidal_angles.push_back(poloidal_angles.front() + PI2);
        // \\theta range: \\theta_0, ..., \\theta_0 + 2\\pi
        intp::InterpolationFunctionTemplate1D<ORDER, val_type, val_type>
            poloidal_template{intp::util::get_range(poloidal_angles),
                              poloidal_angles.size(), true};

        if (std::fpclassify(poloidal_angles.front()) != FP_ZERO) {
            poloidal_angles.insert(poloidal_angles.begin(), 0);
        }

        poloidal_angles.back() = PI2;
        // \\theta range: 0, ..., 2\\pi
        intp::InterpolationFunctionTemplate1D<ORDER, val_type, val_type>
            poloidal_template_full{intp::util::get_range(poloidal_angles),
                                   poloidal_angles.size(), false};

        // output data

        intp::Mesh<val_type, 2> magnetic_boozer(radial_sample, lst + 1);
        intp::Mesh<val_type, 2> r_boozer(radial_sample, lst + 1);
        intp::Mesh<val_type, 2> z_boozer(radial_sample, lst + 1);
        intp::Mesh<val_type, 2> jacobian_boozer(radial_sample, lst + 1);

        std::vector<val_type> safety_factor, pol_current_n, tor_current_n,
            pressure_n, r_minor_n, tor_flux_n;

        safety_factor.reserve(radial_sample);
        pol_current_n.reserve(radial_sample);
        tor_current_n.reserve(radial_sample);
        pressure_n.reserve(radial_sample);
        r_minor_n.reserve(radial_sample);
        tor_flux_n.reserve(radial_sample);

        const val_type B0 = b_field(gfile_data.magnetic_axis, 0.);
        const val_type R0 = gfile_data.magnetic_axis.x();

        // This two basic unit determines the output spdata unit,
        // setting them to 1 means SI unit.
        const val_type length_unit = use_si_ ? 1. : R0;
        const val_type magnetic_field_unit = use_si_ ? 1. : B0;

        const val_type current_unit = length_unit * magnetic_field_unit;
        const val_type pressure_unit =
            magnetic_field_unit * magnetic_field_unit / magnetic_constant;
        const val_type flux_unit =
            length_unit * length_unit * magnetic_field_unit;

        // construct boozer coordinate, integrate B^2 * J along each contour

#define boozer_list() \
    X(r);             \
    X(z);             \
    X(b2j);           \
    X(bp2j)

        for (std::size_t ri = 0; ri < contours.size(); ++ri) {
            const val_type psi = contours[ri].flux() - psi_ma_intp;
            const std::size_t poloidal_size = contours[ri].size() + 1;
#define X(name)                       \
    std::vector<val_type> name##_geo; \
    name##_geo.reserve(poloidal_size)
            boozer_list();
#undef X

            // quantities on geometric grid
            for (size_t i = 0; i < poloidal_size; ++i) {
                const auto& pt = contours[ri][i % (poloidal_size - 1)];
                r_geo.push_back(pt.x());
                z_geo.push_back(pt.y());
                b2j_geo.push_back(b2j_field(pt, psi));
                bp2j_geo.push_back(bp2j_field(pt, psi));
            }
            // interpolation on geometric grid
#define X(name)            \
    auto name##_geo_intp = \
        poloidal_template.interpolate(intp::util::get_range(name##_geo))
            boozer_list();
#undef X

            // integrate

            std::vector<val_type> b2j_int;
            b2j_int.reserve(poloidal_angles.size());
            b2j_int.push_back(0);

            std::vector<val_type> bp2j_int;
            bp2j_int.reserve(poloidal_angles.size());
            bp2j_int.push_back(0);

            // Poloidal grid begins from \\theta = 0 and ends at \\theta = 2\\pi
            for (size_t i = 1; i < poloidal_angles.size(); ++i) {
                b2j_int.push_back(b2j_int.back() +
                                  util::integrate_coarse(b2j_geo_intp,
                                                         poloidal_angles[i - 1],
                                                         poloidal_angles[i]));
                bp2j_int.push_back(bp2j_int.back() + util::integrate_coarse(
                                                         bp2j_geo_intp,
                                                         poloidal_angles[i - 1],
                                                         poloidal_angles[i]));
            }
            const auto b2j_flux_avg = b2j_int.back() / PI2;
            tor_current_n.push_back(bp2j_int.back() / (PI2 * current_unit));
            // normalization
            for (auto& v : b2j_int) { v /= b2j_flux_avg; }
            auto boozer_geo_intp = poloidal_template_full.interpolate(
                intp::util::get_range(b2j_int));

            // calculate necessary values on a even-spaced boozer grid
            for (size_t i = 0; i <= lst; ++i) {
                const auto theta_boozer =
                    (static_cast<val_type>(i % lst) + .5) * theta_delta_;
                auto theta_geo = util::find_root(
                    [&](val_type t) {
                        return boozer_geo_intp(t) - theta_boozer;
                    },
                    val_type{}, PI2);
                auto r_grid = r_geo_intp(theta_geo);
                auto z_grid = z_geo_intp(theta_geo);

                // be careful of normalization
                const auto b = b_field({r_grid, z_grid}, psi);
                magnetic_boozer(ri, i) = b / magnetic_field_unit;
                r_boozer(ri, i) = r_grid / length_unit;
                // z value is shifted such that magnetic axis has z = 0
                z_boozer(ri, i) =
                    (z_grid - gfile_data.magnetic_axis.y()) / length_unit;
                jacobian_boozer(ri, i) =
                    b2j_flux_avg / (b * b) * magnetic_field_unit / length_unit;
            }

            safety_factor.push_back(safety_factor_intp(psi));
            pol_current_n.push_back(poloidal_current_intp(psi) / current_unit);
            pressure_n.push_back(pressure_intp(psi) / pressure_unit);
            tor_flux_n.push_back(
                (ri == 0 ? 0. : tor_flux_n.back()) +
                util::integrate_coarse(
                    safety_factor_intp,
                    static_cast<val_type>(
                        ri == 0 ? 0. : (contours[ri - 1].flux() - psi_ma_intp)),
                    psi) /
                    flux_unit);
            // r_minor defined as distance from magnetic axis at weak field side
            // this value is always normalized to R0
            r_minor_n.push_back(r_geo_intp(0.) / R0 - 1.);
        }

        timer.pause();

        // psi_delta_ is normalized after flux surface is fully constructed, and
        // should never be changed hereafter
        psi_delta_ /= flux_unit;

        const auto q0 = safety_factor_intp(0);
        const auto b0n = B0 / magnetic_field_unit;
        const auto g0n = poloidal_current_intp(0) / current_unit;
        const auto p0n = pressure_intp(0) / pressure_unit;
        return MagneticEquilibriumRaw_{
            {std::move(magnetic_boozer), std::move(r_boozer),
             std::move(z_boozer), std::move(jacobian_boozer)},
            {std::move(safety_factor), std::move(pol_current_n),
             std::move(tor_current_n), std::move(pressure_n),
             std::move(r_minor_n), std::move(tor_flux_n)},
            {b0n, R0 / length_unit, 0., q0 * g0n / (b0n * b0n)},
            {q0, g0n, 0., p0n, 0., 0.},
            flux_unit};
    }

   public:
    template <typename F>
    MagneticEquilibrium(const GFileRawData<val_type>& gfile_data,
                        std::size_t radial_grid_num,
                        std::size_t poloidal_grid_num,
                        bool use_si,
                        std::size_t radial_sample,
                        val_type psi_ratio,
                        const F& generate_psi_for_output)
        : lsp(radial_grid_num),
          lst(poloidal_grid_num),
          use_si_(use_si),
          psi_delta_(psi_ratio *
                     (gfile_data.flux_LCFS - gfile_data.flux_magnetic_axis) /
                     static_cast<val_type>(lsp - 1)),
          theta_delta_(2. * M_PI / static_cast<val_type>(lst)),
          spdata_raw_{generate_boozer_coordinate_(gfile_data,
                                                  radial_sample,
                                                  psi_ratio)},
          spdata_intp_{spdata_raw_, *this, generate_psi_for_output,
                       std::make_index_sequence<FIELD_NUM_2D>{},
                       std::make_index_sequence<FIELD_NUM_1D>{}} {}

    const std::size_t lsp, lst;

    val_type psi_delta() const { return psi_delta_; }
    val_type theta_delta() const { return theta_delta_; }

    const MagneticEquilibriumIntp_& intp_data() const { return spdata_intp_; }

    const std::vector<val_type>& psi_for_output() const {
        return intp_data().psi_sample_for_output;
    }

    const std::array<val_type, FIELD_NUM_2D>& axis_value_2d() const {
        return spdata_raw_.axis_value_2d;
    }
    const std::array<val_type, FIELD_NUM_1D>& axis_value_1d() const {
        return spdata_raw_.axis_value_1d;
    }

   private:
    const bool use_si_;
    val_type psi_delta_;
    const val_type theta_delta_;
    const MagneticEquilibriumRaw_ spdata_raw_;
    const MagneticEquilibriumIntp_ spdata_intp_;

#ifdef MEQ_ZERNIKE_SERIES_
    Zernike::Series<val_type>
#else
    intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>
#endif
    create_2d_spline_(const intp::Mesh<val_type, 2>& data,
                      const std::vector<val_type>& psi_sample) const {
#ifdef MEQ_ZERNIKE_SERIES_
        static_cast<void>(psi_sample);
        // The Zernike series is actually representing f(r, theta+delta/2)

        std::vector<val_type> r(spdata_raw_.data_1d[5]);
        const auto psi_w = r[r.size() - 1];
        for (auto& v : r) { v = std::sqrt(v / psi_w); }

        const auto zernike_order = static_cast<int>(
            lst / 5 > MEQ_MAX_ZERNIKE_POLAR_ORDER ? MEQ_MAX_ZERNIKE_POLAR_ORDER
                                                  : lst / 5);
        return {zernike_order, r.size(), lst, data, r};
#else
        // interpolate the even-spaced data
        intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type> data_intp(
            {false, true}, data,
            std::make_pair(psi_delta(),
                           psi_delta() * static_cast<val_type>(lsp - 1)),
            std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_));

        if (psi_sample.empty()) { return data_intp; }

        // resample on the interpolated function
        intp::Mesh<val_type, 2> data_resampled(lsp, lst + 1);
        for (std::size_t i = 0; i < lsp; ++i) {
            for (std::size_t j = 0; j <= lst; ++j) {
                data_resampled(i, j) =
                    data_intp(psi_sample[i],
                              (static_cast<val_type>(j) + .5) * theta_delta_);
            }
        }
        // construct the interpolation function for output
        return intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>(
            {false, true}, std::move(data_resampled),
            intp::util::get_range(psi_sample),
            std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_));
#endif
    }

    intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>
    create_1d_spline_(const std::vector<val_type>& data,
                      const std::vector<val_type>& psi_sample) const {
        // interpolate the even-spaced data
        intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type> data_intp(
            std::make_pair(psi_delta(),
                           psi_delta() * static_cast<val_type>(lsp - 1)),
            intp::util::get_range(data), false);

        if (psi_sample.empty()) { return data_intp; }

        // resample on the interpolated function
        std::vector<val_type> data_resampled;
        data_resampled.reserve(psi_sample.size());
        for (const auto& psi : psi_sample) {
            data_resampled.push_back(data_intp(psi));
        }

        return intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>{
            intp::util::get_range(psi_sample),
            intp::util::get_range(data_resampled)};
    }
};

#endif  // MEQ_MAGNETIC_EQUILIBRIUM_H
