#ifndef MEQ_GFILE_RAW_DATA_H
#define MEQ_GFILE_RAW_DATA_H

#include <fstream>
#include <sstream>  // stringstream
#include <string>
#include <vector>

#include "BSplineInterpolation/src/include/Interpolation.hpp"
#include "Vec.h"
#include "util.h"

template <typename T = double>
struct GFileRawData {
    using val_type = T;

    std::array<char, 48> metadata;
    std::string extra_metadata;
    std::string code_name;
    std::string date;
    std::string shot_id;
    std::string time;

    // horizontal sample number
    unsigned int nw;
    // vertical sample number
    unsigned int nh;
    // horizontal and vertical dimension
    Vec<2, double> dim;
    // center position of horizontal range
    val_type r_center;
    // left position of horizontal range
    val_type r_left;
    // center position of vertical range
    val_type z_mid;
    // magnetic axis coordinate
    Vec<2, val_type> magnetic_axis;
    // magnetic flux at magnetic axis
    val_type flux_magnetic_axis;
    // magnetic flux at Last Closed Flux Surface
    val_type flux_LCFS;
    // magnetic field strength at center position
    val_type b_center;
    // current
    val_type current;
    // magnetic flux at separatrix
    val_type flux_sep;
    // separatrix coordinate
    Vec<2, val_type> sep;

    // poloidal current density profile
    std::vector<val_type> f_pol;
    // pressure profile
    std::vector<val_type> pressure;
    // f*f^{\prime}
    std::vector<val_type> f_f_prime;
    // p^{\prime}
    std::vector<val_type> p_prime;

    // magnetic flux
    intp::Mesh<val_type, 2> flux;

    // safety factor profile
    std::vector<val_type> safety_factor;

    // boundary point number
    unsigned int boundary_num;
    // limiter pointer number
    unsigned int limiter_num;

    // boundary point coordinates
    std::vector<Vec<2, val_type>> boundary;
    // limiter point coordinates
    std::vector<Vec<2, val_type>> limiter;

    // extra data

    std::vector<val_type> geometric_poloidal_angles;

    // constructors

    GFileRawData() : flux{1, 1}, complete_(false) {}
    GFileRawData(const GFileRawData&) = delete;
    GFileRawData(GFileRawData&&) = default;

    // properties

    bool is_complete() const noexcept { return complete_; }

    // post-process

    void rearrange_boundary() {
        geometric_poloidal_angles.reserve(boundary_num);

        size_t middle = 0;
        for (size_t i = 0; i < boundary_num; ++i) {
            geometric_poloidal_angles.push_back(
                util::arctan(boundary[i] - magnetic_axis));

            if (i > 0 && std::abs(geometric_poloidal_angles[i] -
                                  geometric_poloidal_angles[i - 1]) > M_PI) {
                // middle is the index of the first angle crossing \theta=0
                middle = i;
            }
        }

        std::rotate(boundary.begin(),
                    boundary.begin() + static_cast<std::ptrdiff_t>(middle),
                    boundary.end());
        std::rotate(geometric_poloidal_angles.begin(),
                    geometric_poloidal_angles.begin() +
                        static_cast<std::ptrdiff_t>(middle),
                    geometric_poloidal_angles.end());

        {
            auto last = std::unique(geometric_poloidal_angles.begin(),
                                    geometric_poloidal_angles.end());
            geometric_poloidal_angles.erase(last,
                                            geometric_poloidal_angles.end());
        }
        {
            auto last = std::unique(boundary.begin(), boundary.end());
            boundary.erase(last, boundary.end());
        }

        // make sure boundary points are sorted counterclockwise
        if (geometric_poloidal_angles[0] > geometric_poloidal_angles[1]) {
            std::reverse(boundary.begin(), boundary.end());
            std::reverse(geometric_poloidal_angles.begin(),
                         geometric_poloidal_angles.end());
        }
    }

    // operators

    template <typename U>
    friend std::ifstream& operator>>(std::ifstream& is, GFileRawData<U>& g) {
        using val_type = GFileRawData<U>::val_type;
        std::string s_dum;
        std::string line;
        std::stringstream ss_line;
        val_type d_dum;
        // 1st line, extarcted with a fixed width style to avoid peculiar layout
        // causing problems
        std::getline(is, line);
        ss_line.str(line);
        ss_line.read(g.metadata.begin(), 48);
        char str_tmp[5]{};
        ss_line.read(str_tmp, 4);
        ss_line.read(str_tmp, 4);
        g.nw = static_cast<unsigned>(atoi(str_tmp));
        ss_line.read(str_tmp, 4);
        g.nh = static_cast<unsigned>(atoi(str_tmp));
        if (!ss_line.eof()) {
            // There is extra metadata
            g.extra_metadata =
                std::string(std::istreambuf_iterator<char>(ss_line), {});
        }
        if (ss_line.fail()) {
            std::cout << "Data corruption at 1st line.\n";
            return is;
        }

        // 2nd line
        std::getline(is, line);
        ss_line.clear();
        ss_line.str(line);
        ss_line >> g.dim.x() >> g.dim.y() >> g.r_center >> g.r_left >> g.z_mid;
        if (ss_line.fail()) {
            std::cout << "Data corruption at 2nd line.\n";
            return is;
        }
        // 3rd line
        std::getline(is, line);
        ss_line.clear();
        ss_line.str(line);
        ss_line >> g.magnetic_axis.x() >> g.magnetic_axis.y() >>
            g.flux_magnetic_axis >> g.flux_LCFS >> g.b_center;
        if (ss_line.fail()) {
            std::cout << "Data corruption at 3rd line.\n";
            return is;
        }
        // 4th line
        std::getline(is, line);
        ss_line.clear();
        ss_line.str(line);
        ss_line >> g.current >> d_dum >> d_dum >> d_dum >> d_dum;
        if (ss_line.fail()) {
            std::cout << "Data corruption at 4th line.\n";
            return is;
        }
        // 5th line
        std::getline(is, line);
        ss_line.clear();
        ss_line.str(line);
        ss_line >> d_dum >> d_dum >> g.flux_sep >> g.sep.x() >> g.sep.y();
        if (ss_line.fail()) {
            std::cout << "Data corruption at 5th line.\n";
            return is;
        }

        auto read_vec = [&is](unsigned int count, std::vector<val_type>& arr) {
            arr.reserve(count);
            for (unsigned i = 0; i < count; ++i) {
                val_type d;
                is >> d;
                arr.emplace_back(d);
            }
        };

        // poloidal current
        read_vec(g.nw, g.f_pol);
        if (is.fail()) {
            std::cout << "Data corruption at poloidal current.\n";
            return is;
        }
        // pressure
        read_vec(g.nw, g.pressure);
        if (is.fail()) {
            std::cout << "Data corruption at pressure.\n";
            return is;
        }
        // f*f^{\prime}
        read_vec(g.nw, g.f_f_prime);
        if (is.fail()) {
            std::cout << "Data corruption at f*f'.\n";
            return is;
        }
        // p^{\prime}
        read_vec(g.nw, g.p_prime);
        if (is.fail()) {
            std::cout << "Data corruption at p'.\n";
            return is;
        }
        // flux
        g.flux.resize({g.nw, g.nh});
        for (unsigned j = 0; j < g.nh; ++j) {
            for (unsigned i = 0; i < g.nw; ++i) { is >> g.flux(i, j); }
        }
        if (is.fail()) {
            std::cout << "Data corruption at \\psi(r,z).\n";
            return is;
        }
        // safety factor
        read_vec(g.nw, g.safety_factor);
        if (is.fail()) {
            std::cout << "Data corruption at q.\n";
            return is;
        }
        // # boundary point and # limiter point
        is >> g.boundary_num >> g.limiter_num;
        if (is.fail()) {
            std::cout
                << "Data corruption at boundary number or limiter number.\n";
            return is;
        }
        // boundary points
        for (unsigned i = 0; i < g.boundary_num; ++i) {
            val_type rr, zz;
            is >> rr >> zz;
            g.boundary.emplace_back(rr, zz);
        }
        if (is.fail()) {
            std::cout << "Data corruption at boundary points.\n";
            return is;
        }

        g.rearrange_boundary();

        // limiter points
        for (unsigned i = 0; i < g.limiter_num; ++i) {
            val_type rr, zz;
            is >> rr >> zz;
            g.limiter.emplace_back(rr, zz);
        }
        if (is.fail()) {
            std::cout << "Data corruption at limiter points.\n";
            return is;
        }

        g.complete_ = true;

        return is;
    }

   private:
    // the completeness of raw data
    bool complete_;
};

#endif  // MEQ_GFILE_RAW_DATA_H
