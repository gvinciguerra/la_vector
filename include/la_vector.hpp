// This file is part of la_vector <https://github.com/gvinciguerra/la_vector>.
// Copyright (c) 2020 Giorgio Vinciguerra.
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

#include <limits>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include "piecewise_linear_model.hpp"

/** Computes (bits_per_correction > 0 ? 2^(bits_per_correction-1) - 1 : 0) without the conditional operator. */
#define BPC_TO_EPSILON(bits_per_correction) (((1ul << (bits_per_correction)) + 1) / 2 - 1)

/** Computes the smallest integral power of two that is not smaller than x. */
#define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))

/** Computes the smallest integral value not less than x / y, where x and y must be positive integers. */
#define CEIL_UINT_DIV(x, y) ((x) / (y) + ((x) % (y) != 0))

/** Computes the number of bits needed to store x, that is, 0 if x is 0, 1 + floor(log2(x)) otherwise. */
#define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

template<class ForwardIt, class T, class Compare = std::less<T>>
ForwardIt upper_bound_branchless(ForwardIt first, ForwardIt last, const T &value, Compare comp = Compare()) {
    auto n = std::distance(first, last);
    while (n > 1) {
        auto half = n / 2;
        __builtin_prefetch(&*first + half / 2, 0, 0);
        __builtin_prefetch(&*first + half + half / 2, 0, 0);
        first = !comp(value, *std::next(first, half)) ? first + half : first;
        n -= half;
    }
    return std::next(first, !comp(value, *first));
}

template<typename K, uint8_t t_bpc = 0>
class la_vector {
    static_assert(std::is_integral_v<K>);
    static_assert(std::is_unsigned_v<K>);
    static_assert(t_bpc < sizeof(K) * CHAR_BIT);

    class segment;
    class constant_bpc;
    class variable_bpc;
    class la_iterator;
    friend class la_iterator;

    using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;
    using signed_position_type = typename std::make_signed_t<position_type>;
    using larger_signed_key_type = typename std::conditional_t<sizeof(K) <= 4, int64_t, __int128>;

public:
    using size_type = size_t;
    using iterator = class la_iterator;
    int32_t extra_segments = 0;

private:
    static constexpr bool auto_bpc = t_bpc < 2;
    static constexpr size_t cache_line_bits = 64 * CHAR_BIT;
    static constexpr size_t extraction_density = auto_bpc ? 1 : BIT_CEIL(4 * cache_line_bits / t_bpc);

    // If auto_bpc, each segment uses a different bit-size for the corrections, stored in segment::bpc. It also stores
    // segment::corrections_offset, which is the cumulative sum of the preceding segments' bpc values and is used as a
    // (bit) pointer into the corrections array.
    // If !auto_bpc, the corrections array is split into two parts. The first part contains the corrections at positions
    // multiples of extraction_density. The second part contains the remaining corrections.

    K front;                          ///< The first element in this container.
    K back;                           ///< The last element in this container.
    size_t n;                         ///< The number of elements in this container.
    size_t val_step;                  ///< The chunk size of val_top_level, in terms of universe values 0,...,back.
    size_t pos_step;                  ///< The chunk size of pos_top_level, in terms of positions 0,...,n.
    size_t top_level_size;            ///< The size of the two *_top_level structures.
    size_t total_bits_corrections;    ///< The number of bits needed to store the corrections.
    std::vector<segment> segments;    ///< The linear models that, together with the corrections, compress the data.
    sdsl::int_vector<64> corrections; ///< The corrections for each compressed element.
    sdsl::int_vector<> val_top_level; ///< Used to speed up segment_for_value, contains positions of segments.
    sdsl::int_vector<> pos_top_level; ///< Used to speed up segment_for_position, contains positions of segments.

    using floating_type = double;
    using canonical_segment = typename OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;
    using base_segment_type = typename std::conditional_t<auto_bpc, variable_bpc, constant_bpc>;

    struct canonical_segment_bpc : canonical_segment {
        uint8_t bpc;
        canonical_segment_bpc() = default;
        canonical_segment_bpc(const canonical_segment &cs, uint8_t bpc) : canonical_segment(cs), bpc(bpc) {};
    };

    /**
     * Returns an iterator to the segment responsible for decompressing the element at the given position.
     * @param i position of the element
     * @return an iterator to the segment responsible for the given position
     */
    typename decltype(segments)::const_iterator segment_for_position(size_t i) const {
        auto k = i / pos_step;
        auto first = segments.begin() + (i < pos_step ? 0 : pos_top_level[k - 1]);
        auto last = segments.begin() + pos_top_level[k];
        auto cmp = [](size_t x, const segment &s) { return x < s.first; };
        return std::prev(upper_bound_branchless(first, last, i, cmp));
    }

    /**
     * Returns an iterator to the segment responsible for decompressing an element that is not less than the given
     * value.
     * @param value value of the element
     * @return an iterator to the segment responsible for the given value
     */
    typename decltype(segments)::const_iterator segment_for_value(K value) const {
        auto k = value / val_step;
        auto first = segments.begin() + (value < val_step ? 0 : val_top_level[k - 1]);
        auto last = segments.begin() + val_top_level[k];
        auto cmp = [](K x, const segment &s) { return x < s.first_key(); };
        auto it = upper_bound_branchless(first, last, value, cmp);
        return it == segments.begin() ? it : std::prev(it);
    }

    template<typename RandomIt>
    static std::vector<canonical_segment> standard_segmentation(RandomIt begin, RandomIt end, uint8_t bpc) {
        auto n = std::distance(begin, end);
        auto eps = BPC_TO_EPSILON(bpc);
        std::vector<canonical_segment> out;
        out.reserve(eps > 0 ? n / (eps * eps) : n / 8);
        auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
        auto out_fun = [&out](auto, auto, auto cs) { out.push_back(cs); };
        make_segmentation_par(n, eps, in_fun, out_fun);
        return out;
    }

    template<typename RandomIt>
    static std::vector<canonical_segment_bpc> space_efficient_segmentation(RandomIt begin, RandomIt end) {
        const auto n = size_t(std::distance(begin, end));
        const auto data = begin;
        auto max_bpc = (1 + uint8_t(std::log2(data[n - 1]))) / 2;

        std::vector<std::vector<canonical_segment>> segmentations(max_bpc + 1);

        #pragma omp parallel for default(none) firstprivate(n) shared(begin, end, segmentations, max_bpc)
        for (uint8_t bpc = 0; bpc <= max_bpc; ++bpc) {
            if (bpc == 1)
                continue;
            segmentations[bpc] = standard_segmentation(begin, end, bpc);
            segmentations[bpc].emplace_back(n);
            segmentations[bpc].shrink_to_fit();
        }

        // Ignore redundant segmentations
        size_t new_size = segmentations.size();
        for (; segmentations[new_size - 2].size() == 2; --new_size);
        segmentations.resize(new_size);
        max_bpc = new_size - 1;

        // Find the shortest path
        std::vector<size_t> distance(n + 1, std::numeric_limits<size_t>::max());
        std::vector<size_t> cursors(max_bpc + 1);
        std::vector<std::unique_ptr<canonical_segment_bpc>> predecessor(n + 1);
        size_t i = 0;
        distance[0] = 0;

        while (i < n) {
            size_t next_i = -1;

            // For each j adjacent to i, do a relaxation
            for (uint8_t bpc = 0; bpc <= max_bpc; bpc += 1 + (bpc == 0)) {
                // Move cursors[bpc] so that the segment at segmentations[bpc][cursors[bpc]] covers i
                for (; cursors[bpc] + 1 < segmentations[bpc].size()
                           && segmentations[bpc][cursors[bpc] + 1].get_first_x() <= i; ++cursors[bpc]);

                auto &s = segmentations[bpc][cursors[bpc]];
                auto j = segmentations[bpc][cursors[bpc] + 1].get_first_x();
                auto weight_ij = bpc * (j - i) + CHAR_BIT * sizeof(segment);
                assert(s.get_first_x() <= i);

                // Relax edge (i, j)
                if (distance[j] > distance[i] + weight_ij) {
                    distance[j] = distance[i] + weight_ij;
                    predecessor[j] = std::make_unique<canonical_segment_bpc>(s.copy(i), bpc);
                }

                next_i = std::min<size_t>(next_i, j);
            }

            i = next_i;
        }

        // Traverse the predecessor map to build the result
        if (!predecessor[n])
            throw std::runtime_error("Cannot reach target vertex");

        size_t result_size = 0;
        for (size_t current = n; current != 0; current = predecessor[current]->get_first_x())
            ++result_size;

        std::vector<canonical_segment_bpc> out(result_size);
        for (size_t current = n, c = 0; current != 0; current = predecessor[current]->get_first_x())
            out[result_size - ++c] = *predecessor[current];

        return out;
    }

    template<typename F>
    void fill_top_level(sdsl::int_vector<> &iv, size_t size, uint64_t step, F op) {
        if (step == 0)
            return;

        iv = sdsl::int_vector<>(size, segments.size(), BIT_WIDTH(segments.size()));
        for (auto i = 0, k = 0; i < size - 1; ++i) {
            while (k < segments.size() && op(segments[k]) < (i + 1) * step)
                ++k;
            iv[i] = k;
        }
    }

    template<class RandomIt>
    void push_segment(const canonical_segment &cs, uint8_t bpc, position_type corrections_offset,
                      RandomIt data, size_t i, size_t j) {
        try {
            segment s(cs, bpc, corrections_offset, data, n, i, j, corrections.data());
            segments.push_back(s);
        }
        catch (const std::overflow_error &) {
            ++extra_segments;
            auto half = (i + j) / 2;
            push_segment(cs, bpc, corrections_offset, data, i, half);
            push_segment(cs.copy(half), bpc, corrections_offset + half * bpc, data, half, j);
        }
    }

public:
    la_vector() = default;

    explicit la_vector(std::vector<K> &data) : la_vector(data.begin(), data.end()) {};

    template<class RandomIt>
    la_vector(RandomIt begin, RandomIt end)
        : front(*begin),
          back(*std::prev(end)),
          n(std::distance(begin, end)),
          total_bits_corrections(),
          segments() {
        if (n == 0)
            return;

        using value_type = typename std::conditional_t<auto_bpc, canonical_segment_bpc, canonical_segment>;
        std::vector<value_type> canonical_segments;
        if constexpr (auto_bpc)
            canonical_segments = space_efficient_segmentation(begin, end);
        else
            canonical_segments = standard_segmentation(begin, end, t_bpc);
        segments.reserve(canonical_segments.size() + 1);

        // Allocate corrections array
        if constexpr (auto_bpc) {
            total_bits_corrections = 0;
            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
                auto i = it->get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
                total_bits_corrections += it->bpc * (j - i);
            }
        } else
            total_bits_corrections = n * t_bpc;

        auto corrections_words = CEIL_UINT_DIV(total_bits_corrections, 64);
        corrections = decltype(corrections)(corrections_words + 1, 0);

        // Store segments and fill the corrections array
        size_t corrections_offset = 0;
        for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
            auto i = it->get_first_x();
            auto j = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
            uint8_t bpc = t_bpc;
            if constexpr (auto_bpc)
                bpc = it->bpc;
            push_segment(*it, bpc, corrections_offset, begin, i, j);
            corrections_offset += bpc * (j - i);
        }

        segments.reserve(segments.size() + 1);
        segments[segments.size()] = segment(n); // extra segment to avoid bound checking in decode() and lower_bound()

        // Fill the top-level structures
        top_level_size = std::min(1u << 16, BIT_CEIL(segments.size()));
        auto u = *std::prev(end);
        val_step = CEIL_UINT_DIV(u, top_level_size);
        pos_step = CEIL_UINT_DIV(n, top_level_size);
        fill_top_level(val_top_level, top_level_size, val_step, [begin](const segment &s) { return begin[s.first]; });
        fill_top_level(pos_top_level, top_level_size, pos_step, [](const segment &s) { return s.first; });
    }

    /**
     * Returns the element at the specified position. No bounds checking is performed.
     * @param i position of the element to return
     * @return the element at the specified position
     */
    K operator[](size_t i) const {
        assert(i < n);
        return segment_for_position(i)->decompress(corrections.data(), n, i);
    }

    /**
     * Returns an iterator pointing to the first element that is not less than the given value.
     * @param value value to compare the elements to
     * @return an iterator pointing to the first element that is not less than value
     */
    la_iterator lower_bound(K value) const {
        if (value > back)
            return end();
        if (value <= front)
            return begin();

        auto it = segment_for_value(value);
        auto &s = *it;
        auto &t = *std::next(it);
        auto[pos, bound] = s.approximate_position(value);
        pos = std::clamp(pos, s.first, t.first - 1);

        auto lo = pos <= bound + s.first ? s.first : pos - bound;
        auto hi = std::min(pos + bound + 1, t.first);

        if (!auto_bpc) {
            // Binary search on the samples
            auto sample_lo = CEIL_UINT_DIV(lo, extraction_density);
            auto sample_hi = (hi - 1) / extraction_density + 1;

            while (sample_lo < sample_hi) {
                size_t mid = sample_lo + (sample_hi - sample_lo) / 2;
                if (s.decompress(corrections.data(), n, mid * extraction_density) < value) {
                    sample_lo = mid + 1;
                    lo = mid * extraction_density;
                } else {
                    sample_hi = mid;
                    hi = mid * extraction_density;
                }
            }

            // Binary search on the compressed data
            while (lo < hi) {
                auto mid = lo + (hi - lo) / 2;
                if (s.decompress(corrections.data(), n, mid) < value)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            if (lo == t.first)
                return la_iterator(this, t.first, std::next(it));
            return la_iterator(this, lo, it);
        }

        auto val = s.decompress(corrections.data(), n, pos);
        auto search_forward = val <= value;
        constexpr auto linear_threshold = 2 * cache_line_bits / (auto_bpc ? 4 : t_bpc);

        if (hi - lo <= linear_threshold) {
            if (search_forward)
                while (pos < t.first && val < value)
                    val = s.decompress(corrections.data(), n, ++pos);
            else
                while (pos > s.first && s.decompress(corrections.data(), n, pos - 1) >= value)
                    --pos;

            if (pos == t.first)
                return la_iterator(this, t.first, std::next(it));
            return la_iterator(this, pos, it);
        }

        lo = search_forward ? pos : lo;
        hi = search_forward ? hi : pos;
        auto val_at_lo = search_forward ? val : s.decompress(corrections.data(), n, lo);
        auto val_at_hi = search_forward ? s.decompress(corrections.data(), n, hi - 1) : val;
        auto count = hi - lo;

        if (hi == t.first and value > val_at_hi)
            return la_iterator(this, t.first, std::next(it));

        while (count > linear_threshold) {
            auto x = larger_signed_key_type(value) - val_at_lo;
            auto dx = val_at_hi - val_at_lo;
            auto dy = count - 1;
            auto step = x * dy / dx;
            auto p = lo + step;
            auto val_at_p = s.decompress(corrections.data(), n, p);

            if (value > val_at_p) {
                lo = p + 1;
                count -= step + 1;
                val_at_lo = s.decompress(corrections.data(), n, lo);
                if (val_at_lo >= value) {
                    if (lo == t.first)
                        return la_iterator(this, t.first, std::next(it));
                    return la_iterator(this, lo, it);
                }
            } else {
                hi = p;
                count = step;
                val_at_hi = val_at_p;
            }
        }

        for (; lo < hi && s.decompress(corrections.data(), n, lo) < value; ++lo);

        if (lo == t.first)
            return la_iterator(this, t.first, std::next(it));
        return la_iterator(this, lo, it);
    }

    /**
     * Returns the number of elements in the container that are less than or equal to @p value.
     * @param value value to compare elements to
     * @return the number of elements that are less than or equal to @p value
     */
    size_t rank(K value) const {
        return std::distance(begin(), lower_bound(value));
    }

    /**
     * Returns the i-th smallest element in the container.
     * @param i rank of the element, must be between 1 and @ref size()
     * @return the i-th smallest element
     */
    K select(size_t i) const {
        assert(i > 0 && i <= n);
        return operator[](i - 1);
    }

    /**
     * Decodes all the elements of this container into the memory beginning at @p out. The caller is responsible for
     * allocating enough memory for @p out, that is, at least @ref size() * sizeof(K) bytes.
     * @param out the beginning of the destination of the decoded elements
     */
    void decode(K *out) const {
        for (auto it = segments.begin(); it != segments.end(); ++it) {
            auto &s = *it;
            auto covered = std::next(it)->first - s.first;
            auto significand = s.slope_significand;
            auto exponent = s.slope_exponent;
            auto intercept = s.intercept - BPC_TO_EPSILON(s.bpc);

            #pragma omp simd
            for (auto j = 0; j < covered; ++j)
                out[j] = ((j * significand) >> exponent) + intercept;

            for (auto j = 0; j < covered; ++j)
                out[j] += s.get_correction(corrections.data(), n, j + s.first);

            out += covered;
        }
    }

    /**
     * Decodes all the elements of this container into a vector.
     * @return a vector with the decoded elements
     */
    std::vector<K> decode() const {
        std::vector<K> out(n);
        decode(out.data());
        return out;
    }

    /**
     * Returns an iterator to the first element.
     * @return an iterator to the first element
     */
    iterator begin() const {
        return la_iterator(this, 0, segments.begin());
    }

    /**
     * Returns an iterator to the element following the last element.
     * @return an iterator to the element following the last element
     */
    iterator end() const {
        return la_iterator(this, n, segments.end());
    }

    /**
     * Returns the number of bytes used by this container to encode its elements.
     * @return the size in bytes of this container
     */
    size_t size_in_bytes() const {
        return total_bits_corrections / CHAR_BIT + segments_count() * sizeof(segment)
            + (val_top_level.bit_size() + pos_top_level.bit_size()) / CHAR_BIT;
    }

    /**
     * Returns the number of segments (linear models) used in this container to encode its elements.
     * @return the number of segments in this container
     */
    size_t segments_count() const {
        return segments.size();
    }

    /**
     * Returns the average number of bits per element, that is, @ref size_in_bytes() * 8. / @ref size().
     * @return the average number of bits per element
     */
    double bits_per_element() const {
        return size_in_bytes() * CHAR_BIT / (double) n;
    }

    /**
     * Returns the number of elements in this container.
     * @return the number of elements in this container
     */
    size_t size() const {
        return n;
    }

    /**
     * Serializes the vector to a stream.
     * @param out output stream
     * @param v parent node in the structure tree
     * @param name name of the structure tree node
     * @return the number of bytes written to out
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += sdsl::write_member(n, out, child, "size");
        written_bytes += sdsl::write_member(front, out, child, "front");
        written_bytes += sdsl::write_member(back, out, child, "back");
        written_bytes += sdsl::write_member(val_step, out, child, "val_step");
        written_bytes += sdsl::write_member(pos_step, out, child, "pos_step");
        written_bytes += sdsl::write_member(top_level_size, out, child, "top_level_size");
        written_bytes += sdsl::write_member(total_bits_corrections, out, child, "total_bits_corrections");
        written_bytes += corrections.serialize(out, child, "corrections");
        written_bytes += val_top_level.serialize(out, child, "val_top_level");
        written_bytes += pos_top_level.serialize(out, child, "pos_top_level");
        written_bytes += sdsl::write_member(segments.size(), out, child, "segments_size");
        for (auto &s: segments) {
            out.write((char *) &s, sizeof(segment));
            written_bytes += sizeof(segment);
        }
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /**
     * Loads the vector from a stream.
     * @param in input stream
     */
    void load(std::istream &in) {
        sdsl::read_member(n, in);
        sdsl::read_member(front, in);
        sdsl::read_member(back, in);
        sdsl::read_member(val_step, in);
        sdsl::read_member(pos_step, in);
        sdsl::read_member(top_level_size, in);
        sdsl::read_member(total_bits_corrections, in);
        corrections.load(in);
        val_top_level.load(in);
        pos_top_level.load(in);
        size_t segments_size;
        sdsl::read_member(segments_size, in);
        segments = std::vector<segment>();
        segments.reserve(segments_size + 1);
        for (size_t i = 0; i < segments_size; ++i) {
            segment s;
            in.read((char *) &s, sizeof(segment));
            segments.push_back(s);
        }
        segments[segments.size()] = segment(n);
    }
};

#pragma pack(push, 1)

template<typename K, uint8_t t_bpc>
struct la_vector<K, t_bpc>::constant_bpc {
    static constexpr uint8_t bpc = t_bpc;
    uint32_t first_correction: t_bpc;
    constant_bpc() = default;
    constant_bpc(uint8_t, position_type) : first_correction(0) {};
};

template<typename K, uint8_t t_bpc>
struct la_vector<K, t_bpc>::variable_bpc {
    uint8_t bpc;
    position_type corrections_offset;
    uint32_t first_correction: 16;
    variable_bpc() = default;
    variable_bpc(uint8_t bpc, position_type offset) : bpc(bpc), corrections_offset(offset), first_correction(0) {};
};

template<typename K, uint8_t t_bpc>
struct la_vector<K, t_bpc>::segment : base_segment_type {
    static constexpr auto exponent_bits = 5;
    static constexpr auto significand_bits = sizeof(K) <= 4 ? 32 - exponent_bits : 56 - exponent_bits;
    position_type first;
    signed_position_type intercept;
    uint8_t slope_exponent: exponent_bits;
    uint64_t slope_significand: significand_bits;

    segment() = default;

    template<typename RandomIt>
    segment(const canonical_segment &cs, uint8_t bpc, position_type corrections_offset,
            RandomIt data, size_t n, size_t i, size_t j, uint64_t *corrections)
        : base_segment_type(bpc, corrections_offset) {
        auto epsilon = BPC_TO_EPSILON(bpc);
        auto[cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(cs.get_first_x(), j - i + 1);

        if (BIT_WIDTH(cs_exponent) > exponent_bits || BIT_WIDTH(cs_significand) > significand_bits)
            throw std::overflow_error("Bit fields' sizes are not large enough");

        first = cs.get_first_x();
        intercept = cs_intercept;
        slope_exponent = cs_exponent;
        slope_significand = cs_significand;

        for (auto k = i; k < j; k++) {
            auto error = data[k] - approximate(k);
            auto correction = uint64_t(error + epsilon);
            set_correction(corrections, n, k, correction);
        }
    }

    explicit segment(position_type first)
        : base_segment_type(0, 0),
          first(first),
          slope_significand(0),
          slope_exponent(0),
          intercept(std::numeric_limits<K>::max()) {}

    size_t get_correction_bit_offset(size_t n, size_t i) const {
        if constexpr (auto_bpc)
            return this->corrections_offset + (i - first) * this->bpc;
        if (i % extraction_density == 0)
            return this->bpc * (i / extraction_density);
        return this->bpc * (i + n / extraction_density - i / extraction_density);
    }

    void set_correction(uint64_t *corrections, size_t n, size_t i, uint64_t value) {
        if (this->bpc == 0)
            return;
        if (BIT_WIDTH(value) > this->bpc)
            throw std::overflow_error("Segment correction too large");
        if (i == first)
            this->first_correction = value;

        auto idx = get_correction_bit_offset(n, i);
        sdsl::bits::write_int(corrections + (idx >> 6), value, idx & 0x3F, this->bpc);
    }

    K get_correction(const uint64_t *corrections, size_t n, size_t i) const {
        auto idx = get_correction_bit_offset(n, i);
        return sdsl::bits::read_int(corrections + (idx >> 6u), idx & 0x3F, this->bpc);
    }

    larger_signed_key_type approximate(size_t i) const {
        return ((larger_signed_key_type(slope_significand) * (i - first)) >> slope_exponent) + intercept;
    }

    std::pair<position_type, position_type> approximate_position(const K &value) const {
        auto significand = larger_signed_key_type(slope_significand == 0 ? 1 : slope_significand);
        auto p = ((larger_signed_key_type(value) - intercept) << slope_exponent) / significand + first;
        auto epsilon = larger_signed_key_type(BPC_TO_EPSILON(this->bpc));
        auto bound = 1 + (epsilon << slope_exponent) / larger_signed_key_type(significand);
        return {std::max<larger_signed_key_type>(0, p), bound};
    }

    K first_key() const {
        auto epsilon = BPC_TO_EPSILON(this->bpc);
        auto correction = this->first_correction;
        return intercept + correction - epsilon;
    }

    K decompress(const uint64_t *corrections, size_t n, size_t i) const {
        auto epsilon = BPC_TO_EPSILON(this->bpc);
        auto correction = get_correction(corrections, n, i);
        return approximate(i) + correction - epsilon;
    }
};

#pragma pack(pop)

template<typename K, uint8_t t_bpc>
class la_vector<K, t_bpc>::la_iterator {
    using parent_type = const la_vector<K, t_bpc>;
    using segments_iterator = typename decltype(parent_type::segments)::const_iterator;

    parent_type *p;
    size_t i;
    segments_iterator s_it;

    template<bool Forward = true>
    void move_segment_cursor() {
        bool is_segment_cursor_invalid = s_it == p->segments.end();
        if (is_segment_cursor_invalid) {
            s_it = p->segment_for_position(i);
            return;
        }

        if constexpr (Forward) {
            if (i >= std::next(s_it)->first)
                ++s_it;
        } else {
            if (s_it > p->segments.begin() && i < s_it->first)
                --s_it;
        }
    }

public:
    using value_type = const K;
    using pointer = const K *;
    using reference = const K;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;

    la_iterator() : p(nullptr), i(0), s_it(p->segments.begin()) {}
    la_iterator(parent_type *p, size_t i, segments_iterator s_it) : p(p), i(i), s_it(s_it) {}
    la_iterator(parent_type *p, size_t i) : p(p), i(i), s_it(p->segments.end()) { move_segment_cursor<>(); }

    reference operator*() const { return s_it->decompress(p->corrections.data(), p->n, i); }
    reference operator[](difference_type m) const { return (*p)[i + m]; }

    la_iterator &operator++() {
        ++i;
        move_segment_cursor<true>();
        return *this;
    }

    la_iterator &operator--() {
        --i;
        move_segment_cursor<false>();
        return *this;
    }

    la_iterator operator++(int) {
        la_iterator r(*this);
        ++i;
        move_segment_cursor<true>();
        return r;
    }

    la_iterator operator--(int) {
        la_iterator r(*this);
        --i;
        move_segment_cursor<false>();
        return r;
    }

    la_iterator &operator+=(difference_type d) {
        i += d;
        if (i < s_it->first || i >= std::next(s_it)->first) {
            // Reposition the segment cursor
            s_it = p->segments.end();
            move_segment_cursor<>();
        }
        return *this;
    }

    la_iterator &operator-=(difference_type d) {
        *this += -d;
        return *this;
    }

    la_iterator operator+(difference_type d) const { return la_iterator(p, i + d); }
    la_iterator operator-(difference_type d) const { return la_iterator(p, i - d); }

    difference_type operator-(const la_iterator &r) const { return i - r.i; }

    bool operator<(const la_iterator &r) const { return i < r.i; }
    bool operator<=(const la_iterator &r) const { return i <= r.i; }
    bool operator>(const la_iterator &r) const { return i > r.i; }
    bool operator>=(const la_iterator &r) const { return i >= r.i; }
    bool operator!=(const la_iterator &r) const { return i != r.i; }
    bool operator==(const la_iterator &r) const { return i == r.i; }
};
