#ifndef MEQ_MAGNETIC_EQUILIBRIUM_H
#define MEQ_MAGNETIC_EQUILIBRIUM_H

#include <limits>
#include <type_traits>
#include <utility>

#include "BSplineInterpolation/src/include/DedicatedThreadPool.hpp"
#include "BSplineInterpolation/src/include/Interpolation.hpp"
#include "Contour.h"
#include "GFileRawData.h"
#include "Timer.h"

#ifdef MEQ_ZERNIKE_SERIES_
#include "Zernike.h"
#endif

template <typename T, typename F1, typename F2>
class MagneticEquilibrium {
   public:
    // interpolation order of internal use
    constexpr static std::size_t ORDER = 5;

    using val_type = T;

    struct RawValues1D {
        val_type psi;  // poloidal flux
        val_type f;    // poloidal current
        val_type q;    // safety factor
        val_type p;    // pressure
        val_type b2j;  // b^2 * jacobi
    };
    struct RawValues2D {
        val_type dpdr;  // ... and its R derivative
        val_type dpdz;  // ... also Z derivative
        val_type r;     // major radius, i.e. R in cylindrical coord
        val_type z;     // Z in cylindrical coord
    };

    // F should be type of a function accepting `RawValues1D` and `RawValues2D`
    // (normalized, depends on `use_si_`), returning an array-like element
    // containing all the field needed to be constructed on boozer coordinate
    // grid (psi, theta)
    using field_1d_func_type = F1;
    using field_2d_func_type = F2;

    static constexpr std::size_t field_1d_num = std::tuple_size_v<
        std::invoke_result_t<field_1d_func_type, RawValues1D>>;
    static constexpr std::size_t field_2d_num = std::tuple_size_v<
        std::invoke_result_t<field_2d_func_type, RawValues1D, RawValues2D>>;

    using data_1d_type = std::array<std::vector<val_type>, field_1d_num>;
    using data_2d_type = std::array<intp::Mesh<val_type, 2>, field_2d_num>;

    const std::size_t radial_grid;    // For determining minimal psi_p.
    const std::size_t radial_sample;  // For sampling on radial direction. This
                                      // is the output data dimension.
    const std::size_t poloidal_grid;

    MagneticEquilibrium(const GFileRawData<val_type>& gfile_data,
                        std::size_t radial_grid_num,
                        std::size_t poloidal_grid_num,
                        std::size_t radial_sample_num,
                        val_type psi_ratio,
                        bool use_si,
                        field_1d_func_type field_1d_func,
                        field_2d_func_type field_2d_func)
        : radial_grid(radial_grid_num),
          radial_sample(radial_sample_num),
          poloidal_grid(poloidal_grid_num),
          use_si_(use_si),
          psi_ratio_(psi_ratio),
          psi_delta_{},
          theta_delta_(2. * M_PI / static_cast<val_type>(poloidal_grid)),
          axis_value_{},
          data{field_on_boozer_grid(gfile_data, field_1d_func, field_2d_func)} {
    }

    auto psi_delta() const { return psi_delta_; }
    auto theta_delta() const { return theta_delta_; }

    const auto& data_1d(std::size_t i) const { return data.first[i]; }
    const auto& data_2d(std::size_t i) const { return data.second[i]; }

    auto axis_value_1d(std::size_t i) const { return axis_value_.first[i]; }
    auto axis_value_2d(std::size_t i) const { return axis_value_.second[i]; }

   private:
    const bool use_si_;
    const val_type psi_ratio_;
    val_type psi_delta_;
    const val_type theta_delta_;
    std::pair<std::array<val_type, field_1d_num>,
              std::array<val_type, field_2d_num>>
        axis_value_;
    const std::pair<data_1d_type, data_2d_type> data;

    std::pair<data_1d_type, data_2d_type> field_on_boozer_grid(
        const GFileRawData<val_type>& gfile_data,
        field_1d_func_type field_1d_func,
        field_2d_func_type field_2d_func) {
        constexpr val_type magnetic_constant = 4.e-7 * M_PI;

        if (radial_sample % 2 != 0) {
            throw std::runtime_error(
                "[MagneticEquilibrium] Radial sample point must be even.");
        }

        std::cout << "Toroidal field direction (in cylindrical coordinate "
                     "(X,Z,phi)):\n";
        std::cout << "  Magnetic field: "
                  << (gfile_data.f_pol.front() > 0 ? '+' : '-') << '\n';
        std::cout << "  Current: "
                  << (gfile_data.flux_LCFS > gfile_data.flux_magnetic_axis
                          ? '+'
                          : '-')
                  << '\n';

        auto& timer = Timer::get_timer();
        timer.start(" - B-Spline of poloidal flux");
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

        auto psi_boundary_min = std::numeric_limits<val_type>::infinity();
        for (const auto& pt : gfile_data.boundary) {
            psi_boundary_min = std::min(psi_boundary_min, flux_function(pt));
        }
        std::cout << "The poloidal flux of last closed flux surface is "
                  << gfile_data.flux_LCFS << '\n'
                  << "Minimum of interpolated value at boundary points is "
                  << psi_boundary_min << '\n';

        const val_type psi_bd = gfile_data.flux_LCFS - psi_ma_intp;
        val_type psi_wall = psi_ratio_ * psi_bd;
        if (psi_wall > psi_boundary_min - psi_ma_intp) {
            psi_wall = psi_boundary_min - psi_ma_intp;
            std::cout << "Interpolated flux value at boundary is too small, so "
                         "psi_wall is set to this value.\n";
        }
        psi_delta_ = psi_wall / static_cast<val_type>(radial_grid - 1);

#ifdef MEQ_MULTITHREAD
        auto& thread_pool = intp::DedicatedThreadPool<void>::get_instance(4);
        std::vector<std::future<void>> tasks;
#endif
        timer.pause_last_and_start_next(" - Contour");
        // contours are from \\Delta\\psi to LCFS
        std::vector<Contour<val_type>> contours(radial_sample);
        constexpr std::size_t task_size = 8;
        for (std::size_t i = 0; i < (radial_sample + task_size - 1) / task_size;
             ++i) {
            const auto start = i * task_size;
            const auto finish = start + task_size > radial_sample
                                    ? radial_sample
                                    : start + task_size;
            auto construct_contour = [&, start, finish]() {
                for (std::size_t li = start; li < finish; ++li) {
                    contours[li] = Contour<val_type>(
                        util::lerp(
                            psi_delta_, psi_wall,
                            static_cast<val_type>(li) /
                                static_cast<val_type>(radial_sample - 1)) +
                            psi_ma_intp,
                        flux_function, gfile_data);
                }
            };
#ifdef MEQ_MULTITHREAD
            tasks.push_back(thread_pool.queue_task(construct_contour));
#else
            construct_contour();
#endif
        }
#ifdef MEQ_MULTITHREAD
        for (auto& res : tasks) { res.get(); }
        tasks.clear();
#endif
        timer.pause_last_and_start_next(" - Boozer grid");

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

        // b^2 * J used for boozer coordinate construction
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
        // This template is used for \theta_boozer(\theta_geo) thus not being
        // periodic

        const val_type R0 = gfile_data.magnetic_axis.x();
        const val_type B0 = gfile_data.f_pol.front() / R0;

        // This two basic unit determines the output spdata unit,
        // setting them to 1 means SI unit.
        const val_type length_unit = use_si_ ? 1. : R0;
        const val_type magnetic_field_unit = use_si_ ? 1. : B0;

        const val_type current_unit = length_unit * magnetic_field_unit;
        const val_type pressure_unit =
            magnetic_field_unit * magnetic_field_unit / magnetic_constant;
        const val_type flux_unit =
            length_unit * length_unit * magnetic_field_unit;

        // Output container
        data_1d_type data_1d;
        for (std::size_t fi = 0; fi < data_1d.size(); ++fi) {
            data_1d[fi].resize(radial_sample);
        }
        auto data_2d{([this]<std::size_t... idx>(
                          std::index_sequence<idx...>) -> data_2d_type {
            return {
                (static_cast<void>(idx),
                 intp::Mesh<val_type, 2>(radial_sample, poloidal_grid + 1))...};
        })(std::make_index_sequence<std::tuple_size_v<data_2d_type>>{})};

        for (std::size_t ri = 0;
             ri < (contours.size() + task_size - 1) / task_size; ++ri) {
            const auto start = ri * task_size;
            const auto finish = start + task_size > contours.size()
                                    ? contours.size()
                                    : start + task_size;
            const auto construct_boozer = [&, start, finish]() {
                for (std::size_t li = start; li < finish; ++li) {
                    const val_type psi = contours[li].flux() - psi_ma_intp;
                    const std::size_t poloidal_size = contours[li].size() + 1;

                    RawValues1D raw_values_1d{
                        psi / flux_unit,
                        poloidal_current_intp(psi) / current_unit,
                        safety_factor_intp(psi),
                        pressure_intp(psi) / pressure_unit, val_type{}};

                    std::vector<val_type> b2j_geo;
                    std::vector<val_type> rr_geo;
                    std::vector<val_type> zz_geo;
                    b2j_geo.reserve(poloidal_size);
                    rr_geo.reserve(poloidal_size);
                    zz_geo.reserve(poloidal_size);
                    for (size_t i = 0; i < poloidal_size; ++i) {
                        const auto& pt = contours[li][i % (poloidal_size - 1)];
                        b2j_geo.push_back(b2j_field(pt, psi));
                        rr_geo.push_back(pt.x());
                        zz_geo.push_back(pt.y());
                    }
                    auto b2j_geo_intp = poloidal_template.interpolate(
                        intp::util::get_range(b2j_geo));
                    auto rr_geo_intp = poloidal_template.interpolate(
                        intp::util::get_range(rr_geo));
                    auto zz_geo_intp = poloidal_template.interpolate(
                        intp::util::get_range(zz_geo));
                    std::vector<val_type> b2j_int;
                    b2j_int.reserve(poloidal_angles.size());
                    b2j_int.push_back(0);
                    for (size_t i = 1; i < poloidal_angles.size(); ++i) {
                        b2j_int.push_back(
                            b2j_int.back() +
                            util::integrate_coarse(b2j_geo_intp,
                                                   poloidal_angles[i - 1],
                                                   poloidal_angles[i]));
                    }
                    const auto b2j_flux_avg = b2j_int.back() / PI2;
                    for (auto& v : b2j_int) { v /= b2j_flux_avg; }
                    auto boozer_geo_intp = poloidal_template_full.interpolate(
                        intp::util::get_range(b2j_int));
                    raw_values_1d.b2j =
                        b2j_flux_avg / (magnetic_field_unit * length_unit);
                    // calculate necessary values on a even-spaced boozer
                    // grid
                    for (std::size_t i = 0; i <= poloidal_grid; ++i) {
                        const auto theta_boozer =
                            (static_cast<val_type>(i % poloidal_grid) + .5) *
                            theta_delta_;
                        const auto theta_geo = util::find_root(
                            [&](val_type t) {
                                return boozer_geo_intp(t) - theta_boozer;
                            },
                            val_type{}, PI2);
                        const auto rr = rr_geo_intp(theta_geo);
                        const auto zz = zz_geo_intp(theta_geo);

                        // invoke user-provided function
                        const auto field_2d_vals = field_2d_func(
                            raw_values_1d,
                            RawValues2D{
                                flux_function.derivative({rr, zz}, {1, 0}) /
                                    current_unit,
                                flux_function.derivative({rr, zz}, {0, 1}) /
                                    current_unit,
                                rr / length_unit,
                                (zz - gfile_data.magnetic_axis.y()) /
                                    length_unit});
                        for (std::size_t fi = 0; fi < data_2d.size(); ++fi) {
                            data_2d[fi](li, i) = field_2d_vals[fi];
                        }
                    }
                    const auto field_1d_vals = field_1d_func(raw_values_1d);
                    for (std::size_t fi = 0; fi < data_1d.size(); ++fi) {
                        data_1d[fi][li] = field_1d_vals[fi];
                    }
                }
            };
#ifdef MEQ_MULTITHREAD
            tasks.push_back(thread_pool.queue_task(construct_boozer));
#else
            construct_boozer();
#endif
        }
#ifdef MEQ_MULTITHREAD
        for (auto& res : tasks) { res.get(); }
#endif
        timer.pause();

        // calculate axis values
        ([&](RawValues1D r1d, RawValues2D r2d) {
            r1d.b2j = r1d.f * r1d.q;
            axis_value_.first = field_1d_func(r1d);
            axis_value_.second = field_2d_func(r1d, r2d);
        })({0., gfile_data.f_pol.front() / current_unit,
            gfile_data.safety_factor.front(),
            gfile_data.pressure.front() / pressure_unit, 0.},
           {0., 0., gfile_data.magnetic_axis.x() / length_unit, 0.});

        // psi_delta_ is normalized after flux surface is fully constructed, and
        // should never be changed hereafter
        psi_delta_ /= flux_unit;

        return {data_1d, data_2d};
    }

    // #ifdef MEQ_ZERNIKE_SERIES_
    //     Zernike::Series<val_type>
    // #else
    //     intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>
    // #endif
    //     create_2d_spline_(const intp::Mesh<val_type, 2>& data,
    //                       const std::vector<val_type>& psi_sample) const {
    // #ifdef MEQ_ZERNIKE_SERIES_
    //         static_cast<void>(psi_sample);
    //         // The Zernike series is actually representing f(r,
    //         theta+delta/2)
    //
    //         std::vector<val_type> r(spdata_raw_.data_1d[5]);
    //         const auto psi_w = r[r.size() - 1];
    //         for (auto& v : r) { v = std::sqrt(v / psi_w); }
    //
    //         const auto zernike_order = static_cast<int>(
    //             lst / 5 > MEQ_MAX_ZERNIKE_POLAR_ORDER ?
    //             MEQ_MAX_ZERNIKE_POLAR_ORDER
    //                                                   : lst / 5);
    //         return {zernike_order, r.size(), lst, data, r};
    // #else
    //         // interpolate the even-spaced data
    //         intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>
    //         data_intp(
    //             {false, true}, data,
    //             std::make_pair(psi_delta(),
    //                            psi_delta() * static_cast<val_type>(lsp - 1)),
    //             std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 *
    //             theta_delta_));
    //
    //         if (psi_sample.empty()) { return data_intp; }
    //
    //         // resample on the interpolated function
    //         intp::Mesh<val_type, 2> data_resampled(lsp, lst + 1);
    //         for (std::size_t i = 0; i < lsp; ++i) {
    //             for (std::size_t j = 0; j <= lst; ++j) {
    //                 data_resampled(i, j) =
    //                     data_intp(psi_sample[i],
    //                               (static_cast<val_type>(j) + .5) *
    //                               theta_delta_);
    //             }
    //         }
    //         // construct the interpolation function for output
    //         return intp::InterpolationFunction<val_type, 2, ORDER_OUT,
    //         val_type>(
    //             {false, true}, std::move(data_resampled),
    //             intp::util::get_range(psi_sample),
    //             std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 *
    //             theta_delta_));
    // #endif
    //     }
    //
    //     intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>
    //     create_1d_spline_(const std::vector<val_type>& data,
    //                       const std::vector<val_type>& psi_sample) const {
    //         // interpolate the even-spaced data
    //         intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>
    //         data_intp(
    //             std::make_pair(psi_delta(),
    //                            psi_delta() * static_cast<val_type>(lsp - 1)),
    //             intp::util::get_range(data), false);
    //
    //         if (psi_sample.empty()) { return data_intp; }
    //
    //         // resample on the interpolated function
    //         std::vector<val_type> data_resampled;
    //         data_resampled.reserve(psi_sample.size());
    //         for (const auto& psi : psi_sample) {
    //             data_resampled.push_back(data_intp(psi));
    //         }
    //
    //         return intp::InterpolationFunction1D<ORDER_OUT, val_type,
    //         val_type>{
    //             intp::util::get_range(psi_sample),
    //             intp::util::get_range(data_resampled)};
    //     }
};

// template <typename T, typename F1, typename F2>
// auto construct_magnetic_equilibrium(const GFileRawData<T>& gfile_data,
//                                     std::size_t radial_grid_num,
//                                     std::size_t poloidal_grid_num,
//                                     std::size_t radial_sample_num,
//                                     bool use_si,
//                                     F1 field_1d_func,
//                                     F2 field_2d_func) {
//     return MagneticEquilibrium<T, F1, F2>(gfile_data, radial_grid_num,
//                                           poloidal_grid_num,
//                                           radial_sample_num, use_si,
//                                           field_1d_func, field_2d_func);
// }

#endif  // MEQ_MAGNETIC_EQUILIBRIUM_H
