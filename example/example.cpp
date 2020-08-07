#include "utils/hash/murmur2_hash.hpp"

#include "include/cuckoo_dysect.h"
#include <unordered_map>


int main(int, char**)
{
    size_t capacity = 100000;

    using key_type    = size_t;
    using mapped_type = size_t;

    // Our hash tables can be used very similar to std::unordered_map
    std::unordered_map<key_type, mapped_type>    standard_map(capacity);

    // There is just one additional parameter that controls the space efficiency
    dysect::cuckoo_dysect<key_type, mapped_type> dysect_standard(capacity, 1.1);
    // This will instantiate a hash table, that will never be larger than:
    //     1.1*sizeof(std::pair<key_type,mapped_type>) + O(1)

    // If you want to have some control on template parameters, following
    dysect::cuckoo_dysect<
        key_type, mapped_type,
        utils_tm::hash_tm::murmur2_hash, // hash function, the result has to have 64 significant bits
        dysect::cuckoo_config<8,   // bucket size
                              3,   // number of hash functions
                              256, // number of subtables
                              dysect::cuckoo_displacement::bfs // displacement algorithm
                              >
        > dysect_special(capacity,  // obvious
                         1.05,      // memory constraint factor
                         10000);    // probing depth until an insertion is declared unsuccessful


    // INSERT
    auto [it, b] = dysect_standard.insert(5, 8);

    // b is a bool indicating that the element has been inserted
    if (b) std::cout << "key 5 - inserted - mapped data 8" << std::endl;
    else   std::cout << "key 5 - insert unsuccessful (already used?)" << std::endl;

    // it can be used to iterate over the table, or to access the element
    if (it == dysect_standard.end())
        std::cout << "key 5 - insert - key not present"
                  << " (found no free space - probing depth?)" << std::endl;
    else
        std::cout << "key 5 - insert - mapped data is " << it->second
                  << " (interesting if it was contained before)" << std::endl;



    // FIND
    auto it2 = dysect_standard.find(7);

    // the iterator is similar to before
    if (it2 == dysect_standard.end())
        std::cout << "key 7 was not found!" << std::endl;
    else
        std::cout << "key 7 was found contained " << it2->second << std::endl;



    // ACCESSOR
    dysect_standard[42] = 666;
    // this inserts 666 with key 42 or replaces the previous value with that key
}
