#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template<class Parent>
class dstrat_nbfs
{
private:
    using key_type       = typename Parent::key_type;
    using mapped_type    = typename Parent::mapped_type;
    using value_intern   = std::pair<key_type,mapped_type>;
    using Parent_t       = typename Parent::This_t;
    using Hashed_t       = typename Parent::Hashed_t;
    using Bucket_t       = typename Parent::Bucket_t;


    using BFSItem        = std::tuple<value_intern*, int, Bucket_t*>;
    using BFSQueue       = std::vector<BFSItem>;

    Parent&      tab;
    const size_t steps;
    static constexpr size_t nh = Parent_t::nh;

public:
    dstrat_nbfs(Parent& parent, size_t steps = 256, size_t = 0)
        : tab(parent), steps(steps+1)
    { /* parameter is for symmetry with "rwalk" therefore unused*/ }

    dstrat_nbfs(Parent& parent, dstrat_nbfs&& rhs)
        : tab(parent), steps(rhs.steps)
    { }

    inline std::pair<int, value_intern*>
    insert(value_intern t, Hashed_t hash)
    {
        BFSQueue  bq;
        Bucket_t* b[nh];

        tab.getBuckets(hash, b);

        value_intern t_copy = t;

        for (size_t i = 0; i < nh; ++i)
        {
            bq.push_back(BFSItem(&t_copy, -1, b[i]));
        }

        for (size_t i = 0; i < steps; ++i)
        {
            if (expand(bq, i))
            {
                value_intern* pos = rollBackDisplacements(bq);
                return std::make_pair((pos) ? bq.size()-nh : -1, pos);
            }
        }

        return std::make_pair(-1, nullptr);
    }

private:
    inline bool expand(BFSQueue& q, size_t index)
    {
        Bucket_t* b = std::get<2>(q[index]);

        for (size_t i = 0; i < tab.bs && q.size() < steps; ++i)
        {
            value_intern* current = &(b->elements[i]);
            key_type      k       = current->first;

            auto hash = tab.hasher(k);

            Bucket_t* ptr[nh];
            tab.getBuckets(hash, ptr);
            for (size_t ti = 0; ti < nh; ++ti)
            {
                if (ptr[ti] != b) // POTENTIAL BUG!!! continous bucket problem
                {
                    //if (overwatch) std::cout << i << " " << ti << std::endl;
                    q.emplace_back(current, index, ptr[ti]);
                    if (ptr[ti]->space()) return true;
                }
            }
        }
        return false;
    }

    inline value_intern* rollBackDisplacements(BFSQueue& bq)
    {
        value_intern* k1;
        int           prev1;
        Bucket_t*     b1;
        std::tie(k1, prev1, b1) = bq[bq.size()-1];

        value_intern* t = b1->probePtr(key_type()).second;
        (*t) = *k1;

        value_intern* k2;
        int           prev2;
        Bucket_t*     b2;

        value_intern* k0;
        while (prev1 >= 0)
        {
            std::tie(k2,prev2,b2) = bq[prev1];
            (*k1) = *k2;

            k0 = k1;
            k1 = k2; prev1 = prev2; b1 = b2;
        }

        return k0;
    }

};
