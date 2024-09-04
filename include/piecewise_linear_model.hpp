// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
#warning Compilation with -fopenmp is recommended
typedef int omp_int_t;
inline omp_int_t omp_get_max_threads() { return 1; }
#endif

#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>

template<typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template<typename X, typename Y>
class OptimalPiecewiseLinearModel {
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Slope {
        SX dx{};
        SY dy{};

        bool operator<(const Slope &p) const {
            return dy * p.dx < dx * p.dy;
        }

        bool operator>(const Slope &p) const {
            return dy * p.dx > dx * p.dy;
        }

        bool operator==(const Slope &p) const {
            return dy * p.dx == dx * p.dy;
        }

        bool operator!=(const Slope &p) const {
            return dy * p.dx != dx * p.dy;
        }

        explicit operator long double() const {
            return dy / (long double) dx;
        }
    };

    struct StoredPoint {
        X x;
        Y y;
    };

    struct Point {
        X x{};
        Y y{};

        Slope operator-(const Point &p) const {
            return {SX(x) - p.x, SY(y) - p.y};
        }
    };

    const Y epsilon;
    std::vector<Point> lower;
    std::vector<Point> upper;
    X first_x = 0;
    X last_x = 0;
    size_t lower_start = 0;
    size_t upper_start = 0;
    size_t points_in_hull = 0;
    Point rectangle[4];

    auto cross(const Point &O, const Point &A, const Point &B) const {
        auto OA = A - O;
        auto OB = B - O;
        return (OA.dx * OB.dy) - (OA.dy * OB.dx);
    }

public:

    class CanonicalSegment;

    explicit OptimalPiecewiseLinearModel(Y epsilon) : epsilon(epsilon), lower(epsilon), upper(epsilon) {
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");

        upper.reserve(1u << 16);
        lower.reserve(1u << 16);
    }

    bool add_point(const X &x, const Y &y) {
        if (points_in_hull > 0 && x <= last_x)
            throw std::logic_error("Points must be increasing by x.");

        last_x = x;
        auto max_y = std::numeric_limits<Y>::max();
        auto min_y = std::numeric_limits<Y>::lowest();
        Point p1{x, y >= max_y - epsilon ? max_y : y + epsilon};
        Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

        if (points_in_hull == 0) {
            first_x = x;
            rectangle[0] = p1;
            rectangle[1] = p2;
            upper.clear();
            lower.clear();
            upper.push_back(p1);
            lower.push_back(p2);
            upper_start = lower_start = 0;
            ++points_in_hull;
            return true;
        }

        if (points_in_hull == 1) {
            rectangle[2] = p2;
            rectangle[3] = p1;
            upper.push_back(p1);
            lower.push_back(p2);
            ++points_in_hull;
            return true;
        }

        if (epsilon == 0) {
            auto p1_on_line1 = p1 - rectangle[0] == rectangle[2] - rectangle[0];
            points_in_hull = p1_on_line1 ? points_in_hull + 1 : 0;
            return p1_on_line1;
        }

        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];
        bool outside_line1 = p1 - rectangle[2] < slope1;
        bool outside_line2 = p2 - rectangle[3] > slope2;

        if (outside_line1 || outside_line2) {
            points_in_hull = 0;
            return false;
        }

        if (p1 - rectangle[1] < slope2) {
            // Find extreme slope
            auto min = lower[lower_start] - p1;
            auto min_i = lower_start;
            for (auto i = lower_start + 1; i < lower.size(); i++) {
                auto val = (lower[i] - p1);
                if (val > min)
                    break;
                min = val;
                min_i = i;
            }

            rectangle[1] = lower[min_i];
            rectangle[3] = p1;
            lower_start = min_i;

            // Hull update
            auto end = upper.size();
            for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end);
            upper.resize(end);
            upper.push_back(p1);
        }

        if (p2 - rectangle[0] > slope1) {
            // Find extreme slope
            auto max = upper[upper_start] - p2;
            auto max_i = upper_start;
            for (auto i = upper_start + 1; i < upper.size(); i++) {
                auto val = (upper[i] - p2);
                if (val < max)
                    break;
                max = val;
                max_i = i;
            }

            rectangle[0] = upper[max_i];
            rectangle[2] = p2;
            upper_start = max_i;

            // Hull update
            auto end = lower.size();
            for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end);
            lower.resize(end);
            lower.push_back(p2);
        }

        ++points_in_hull;
        return true;
    }

    CanonicalSegment get_segment() const {
        if (points_in_hull == 1)
            return CanonicalSegment(rectangle[0], rectangle[1], first_x);
        return CanonicalSegment(rectangle, first_x);
    }

    void reset() {
        points_in_hull = 0;
        lower.clear();
        upper.clear();
    }
};

template<typename X, typename Y>
class OptimalPiecewiseLinearModel<X, Y>::CanonicalSegment {
    friend class OptimalPiecewiseLinearModel;

    Point rectangle[4];
    X first;

    CanonicalSegment(const Point &p0, const Point &p1, X first) : rectangle{p0, p1, p0, p1}, first(first) {};

    CanonicalSegment(const Point (&rectangle)[4], X first)
        : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, first(first) {};

    bool one_point() const {
        return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
            && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
    }

public:

    CanonicalSegment() = default;

    explicit CanonicalSegment(X first) : CanonicalSegment({first, 0}, {first, 0}, first) {};

    X get_first_x() const { return first; }

    CanonicalSegment copy(X x) const {
        auto c(*this);
        c.first = x;
        return c;
    }

    std::tuple<int64_t, uint8_t, SY> get_fixed_point_segment(X origin, X max_input) const {
        if (one_point())
            return {0, 0, (rectangle[0].y + rectangle[1].y) / 2};

        auto &p1 = rectangle[1];
        auto max_slope = rectangle[3] - rectangle[1];

        auto is_slope_integral = max_slope.dy % max_slope.dx == 0;
        auto slope_exponent = is_slope_integral ? 0 : (uint8_t) std::ceil(std::log2(max_input)) + 1;
        auto slope_significand = (max_slope.dy << slope_exponent) / max_slope.dx;

        auto intercept_n = max_slope.dy * (SX(origin) - p1.x);
        auto intercept_d = max_slope.dx;
        auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
        auto intercept = (intercept_n + rounding_term) / intercept_d + p1.y;

        return {slope_significand, slope_exponent, intercept};
    }
};

template<typename Fin, typename Fout>
size_t make_segmentation(size_t n, size_t epsilon, Fin in, Fout out) {
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t c = 0;
    size_t start = 0;
    auto p = in(0);

    OptimalPiecewiseLinearModel<X, Y> opt(epsilon);
    opt.add_point(p.first, p.second);

    for (size_t i = 1; i < n; ++i) {
        auto next_p = in(i);
        if (i != start && next_p.first == p.first)
            continue;
        p = next_p;
        if (!opt.add_point(p.first, p.second)) {
            out(opt.get_segment());
            start = i;
            --i;
            ++c;
        }
    }

    out(opt.get_segment());
    return ++c;
}

template<typename Fin, typename Fout>
size_t make_segmentation_par(size_t n, size_t epsilon, Fin in, Fout out) {
    auto parallelism = std::min<size_t>(omp_get_max_threads(), 20);
    auto chunk_size = n / parallelism;
    auto c = 0ull;

    if (parallelism == 1 || n < 1ull << 15)
        return make_segmentation(n, epsilon, in, out);

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    using canonical_segment = typename OptimalPiecewiseLinearModel<X, Y>::CanonicalSegment;
    std::vector<std::vector<canonical_segment>> results(parallelism);

    #pragma omp parallel for reduction(+:c) num_threads(parallelism)
    for (auto i = 0ull; i < parallelism; ++i) {
        auto first = i * chunk_size;
        auto last = i == parallelism - 1 ? n : first + chunk_size;
        if (first > 0) {
            for (; first < last; ++first)
                if (in(first).first != in(first - 1).first)
                    break;
            if (first == last)
                continue;
        }

        auto in_fun = [in, first](auto j) { return in(first + j); };
        auto out_fun = [&results, i](const auto &cs) { results[i].emplace_back(cs); };
        results[i].reserve(chunk_size / (epsilon > 0 ? epsilon * epsilon : 16));
        c += make_segmentation(last - first, epsilon, in_fun, out_fun);
    }

    for (auto &v : results)
        for (auto &cs : v)
            out(cs);

    return c;
}