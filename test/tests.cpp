#include <cstdio>
#include <string>
#include <unordered_set>
#include "catch.hpp"
#include "la_vector.hpp"
#include "sdsl/sd_vector.hpp"

#define TEMPLATE_ARGS ((typename T, uint8_t B), T, B), \
    (uint32_t, 0), (uint32_t, 4), (uint32_t, 8), \
    (uint64_t, 0), (uint64_t, 4), (uint64_t, 8)

template<typename T>
std::vector<T> generate_unique(size_t n, int seed = 42) {
    std::uniform_int_distribution<T> distribution(0, 10000000);
    std::vector<T> data(n);
    std::mt19937 generator(seed);
    size_t i = 0;
    size_t u = distribution.max() - distribution.min();
    for (auto k = 0; k < u && i < n; ++k)
        if (generator() % (u - k) < n - i)
            data[i++] = k + distribution.min();
    return data;
}

TEMPLATE_TEST_CASE_SIG("Decode", "", TEMPLATE_ARGS) {
    auto data = generate_unique<T>(1000000);
    REQUIRE(data == la_vector<T, B>(data.begin(), data.end()).decode());
}

TEMPLATE_TEST_CASE_SIG("Iterator", "", TEMPLATE_ARGS) {
    auto data = generate_unique<T>(100000);
    size_t i = 0;
    la_vector<T, B> v(data.begin(), data.end());
    for (auto it = v.begin(); it != v.end(); ++it)
        REQUIRE(data[i++] == *it);
}

TEMPLATE_TEST_CASE_SIG("Select", "", TEMPLATE_ARGS) {
    auto data = generate_unique<T>(100000);
    la_vector<T, B> v(data.begin(), data.end());
    for (auto i = 1; i <= 100000; ++i)
        REQUIRE(data[i - 1] == v.select(i));
}

TEMPLATE_TEST_CASE_SIG("Rank", "", TEMPLATE_ARGS) {
    auto data = generate_unique<T>(100000);
    auto queries = generate_unique<T>(1000, 8);

    la_vector<T, B> v(data.begin(), data.end());
    sdsl::sd_vector<> ef(data.begin(), data.end());
    sdsl::sd_vector<>::rank_1_type ef_rank;
    sdsl::util::init_support(ef_rank, &ef);

    for (auto q : queries)
        REQUIRE(ef_rank(q) == v.rank(q));
}

TEMPLATE_TEST_CASE_SIG("Serialize", "", TEMPLATE_ARGS) {
    auto data = generate_unique<T>(100000);
    std::string tmp = std::tmpnam(nullptr);

    {
        la_vector<T, B> v(data.begin(), data.end());
        sdsl::store_to_file(v, tmp);
    }

    la_vector<T, B> v;
    sdsl::load_from_file(v, tmp);
    REQUIRE(data == v.decode());
    std::remove(tmp.c_str());
}