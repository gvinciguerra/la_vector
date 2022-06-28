#include <iostream>
#include "la_vector.hpp"

int main() {
    // Generate some random data
    std::vector<uint32_t> data(10000000);
    std::generate(data.begin(), data.end(), std::rand);
    std::sort(data.begin(), data.end());

    // Construct the la_vector by allocating 6 bits per correction (bpc)
    constexpr auto bpc = 6;
    la_vector<uint32_t, bpc> v1(data);
    std::cout << "Elements: " << v1.size() << std::endl
              << "Bytes: " << v1.size_in_bytes() << std::endl;

    // Construct the la_vector by minimising the space
    la_vector_opt<uint32_t> v2(data);
    std::cout << "Bytes: " << v2.size_in_bytes() << std::endl;

    //
    // Use la_vector as a compressed bitvector
    //

    v1.rank(500);  // number of elements <= 500
    v1.select(10); // 10th smallest element

    //
    // Use la_vector as a compressed C++ container
    //

    std::cout << "First 10 elements: ";
    for (auto it = v1.begin(); it < v1.begin() + 10; ++it)
        std::cout << *it << ", ";

    v1[14];               // element at position 14
    *v1.lower_bound(125); // iterator to the first element >= 125
    v1.decode();          // decompress v1 into a new std::vector

    return 0;
}
