#pragma once


#include <iostream>
#include <vector>
#include <tuple>

template<class Parent>
class dstrat_bfs
{
public:
    using Key            = typename Parent::Key;
    using Data           = typename Parent::Data;
    using Parent_t       = typename Parent::This_t;
    using Hashed_t       = typename Parent::Hashed_t;
    using Bucket_t       = typename Parent::Bucket_t;

    using BFSQueue       = std::vector<std::tuple<Key, int, Bucket_t*> >;

    Parent&      tab;
    const size_t steps;
    static constexpr size_t nh = Parent_t::nh;

    dstrat_bfs(Parent& parent, size_t steps = 256, size_t = 0)
        : tab(parent), steps(steps+1)
    {
        /* parameter is for symmetry with "rwalk" therefore unused*/
    }

    dstrat_bfs(Parent& parent, dstrat_bfs&& rhs)
        : tab(parent), steps(rhs.steps)
    { }

    inline bool expand(BFSQueue& q, size_t index)
    {
        //bool overwatch = false;
        //if (std::get<0>(q[0]) == 480139307921249140) overwatch = true;
        //if (overwatch) std::cout << "!" << std::endl;
        Bucket_t* b = std::get<2>(q[index]);

        for (size_t i = 0; i < tab.bs && q.size() < steps; ++i)
        {
            Key k = b->get(i).first;

            auto hash = tab.hasher(k);

            Bucket_t* ptr[nh];
            tab.getBuckets(hash, ptr);
            for (size_t ti = 0; ti < nh; ++ti)
            {
                if (ptr[ti] != b)
                {
                    //if (overwatch) std::cout << i << " " << ti << std::endl;
                    q.emplace_back(k, index, ptr[ti]);
                    if (ptr[ti]->space()) return true;
                }
            }
        }
        return false;
    }


    inline bool rollBackDisplacements(std::pair<Key,Data> t, BFSQueue& bq)
    {
        Key       k1;
        int       prev1;
        Bucket_t* b1;
        std::tie(k1, prev1, b1) = bq[bq.size()-1];

        Key       k2;
        int       prev2;
        Bucket_t* b2;
        while (prev1 >= 0)
        {
            std::tie(k2,prev2,b2) = bq[prev1];

            auto pop = b2->pop(k1);
            // if (!pop.first)
            // {
            //     std::cout << "serious issue with rollBack " << k1 << " from " << b2
            //               << " pop " << pop.first << " " << pop.second << std::endl;

            //     return false;
            // }
            // if (!
            b1->insert(k1, pop.second); //)
            // {
            //     std::cout << "even more serious issue with rollBack" << std::endl;
            //     return false;
            // }

            k1 = k2; prev1 = prev2; b1 = b2;
        }

        // if (!
        b1->insert(t); //)
        // {   std::cout << "failed final insert" << std::endl; return false; }

        return true;
    }

    inline int insert(std::pair<Key,Data> t, Hashed_t hash)
    {
        BFSQueue  bq;

        Bucket_t* b[nh];
        tab.getBuckets(hash, b);

        for (size_t i = 0; i < nh; ++i)
        {
            bq.push_back(std::tuple<Key, int, Bucket_t*>(t.first, -1, b[i]));
        }

        for (size_t i = 0; i < steps; ++i)
        {
            /*
            if (i >= bq.size())
            {
                std::cout << "WTF" << std::endl;
                for (size_t i = 0; i < bq.size(); ++i)
                {
                    Key k;
                    int p;
                    Bucket_t* b;
                    std::tie(k,p,b) =  bq[i];
                    std::cout << "| " << k << " | " << p << " | " << b << " |" << std::endl;
                }
                for (size_t i = 0; i < nh; ++i)
                {
                    for (size_t j = 0; j < tab.bs; ++j)
                    {
                        Key k = b[i]->get(j).first;
                        Bucket_t* b2[nh];
                        tab.getBuckets(tab.hasher(k), b2);
                        std::cout << "| k " << k;
                        for (size_t z = 0; z < nh; ++z)
                        {
                            std::cout << " | " << b2[z];
                        }
                        std::cout << " |" << std::endl;
                    }
                }
            }
            else
            */
            if (expand(bq, i))
            {
                return (rollBackDisplacements(t, bq)) ? bq.size()-nh : -1;
            }
        }

        return -1;
    }
};
