#pragma once



static constexpr size_t ct_log(size_t k)
{ return (k-1) ? 1+ct_log(k>>1) : 0; }

template<class Parent, size_t NH>
struct ExternLocExtractor;

template<class Parent, size_t NH>
struct ExternTabExtractor;

template<size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
union ExternHashed;

template <class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
class Hasher
{
private:
    /* auxiliary definitions **************************************************/
    using This_t    = Hasher<Key, HFCT, TAB_WIDTH, LOC_WIDTH, NSTUFF, NWORDS>;
    using HashFct_t = HFCT;

public:
    using Hashed_t = ExternHashed<TAB_WIDTH, LOC_WIDTH, NSTUFF, NWORDS>;

    static_assert (sizeof(Hashed_t) == NWORDS*8, "HashSplitter Bigger Than Expected");

    template<size_t NH>
    using TabExtractor = ExternTabExtractor<This_t, NH>;
    template<size_t NH>
    using LocExtractor = ExternLocExtractor<This_t, NH>;

    /* hasher (itself) ********************************************************/
private:
    HashFct_t fct[NWORDS];

public:

    Hasher()
    {
        for (size_t i = 0; i < NWORDS; ++i)
        {
            fct[i] = HashFct_t(i*8768656543548765336ull);
        }
    }

    Hashed_t operator()(Key k) const
    {
        Hashed_t result;
        for (size_t i = 0; i < NWORDS; ++i)
        {
            result.hash[i] = fct[i](k);
        }
        return result;
    }

};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
struct ExternTabExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,NSTUFF,NWORDS>, NSTUFF>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,NSTUFF,NWORDS>;
    using Hashed_t = typename Parent_t::Hashed_t;

    //static_assert(NH == NSTUFF, "wrong sized extractor");
    static inline size_t tab(Hashed_t h, size_t i) { return h.pair[i].tab; }
};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
struct ExternLocExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,NSTUFF,NWORDS>, NSTUFF>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,NSTUFF,NWORDS>;
    using Hashed_t = typename Parent_t::Hashed_t;

    //static_assert(NH == NSTUFF, "wrong sized extractor");
    static inline size_t loc(Hashed_t h, size_t i) { return h.pair[i].loc; }
};



template <class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH>
class Hasher<Key, HFCT, TAB_WIDTH, LOC_WIDTH, 2, 1>
{
private:
    /* auxiliary definitions **************************************************/
    using This_t    = Hasher<Key, HFCT, TAB_WIDTH, LOC_WIDTH, 2, 1>;
    using HashFct_t = HFCT;

public:

    using Hashed_t = ExternHashed<TAB_WIDTH, LOC_WIDTH, 2, 1>;
    //static_assert (sizeof(Hashed_t) == 16, "HashSplitter Bigger Than Expected");

    template<size_t NH>
    using TabExtractor = ExternTabExtractor<This_t, NH>;
    template<size_t NH>
    using LocExtractor = ExternLocExtractor<This_t, NH>;

    /* hasher (itself) ********************************************************/
private:
    HashFct_t fct;

public:
    Hashed_t operator()(Key k) const
    {
        Hashed_t result;
        result.hash[0] = fct(k);
        return result;
    }
};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NH>
struct ExternTabExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>, NH>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>;
    using Hashed_t = typename Parent_t::Hashed_t;

    static constexpr size_t tab_bitmask = (1ull << TAB_WIDTH) - 1;

    static inline size_t tab(Hashed_t h, size_t i)
    { return (h.pair[0].tab + i*h.pair[1].tab) & tab_bitmask; }
};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH>
struct ExternTabExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>, 2>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>;
    using Hashed_t = typename Parent_t::Hashed_t;

    static constexpr size_t tab_bitmask = (1ull << TAB_WIDTH) - 1;

    static inline size_t tab(Hashed_t h, size_t i) { return h.pair[i].tab; }
};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NH>
struct ExternLocExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>, NH>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>;
    using Hashed_t = typename Parent_t::Hashed_t;

    static constexpr size_t loc_bitmask = (1ull << LOC_WIDTH) - 1;

    static inline size_t loc(Hashed_t h, size_t i)
    { return (h.pair[0].loc + i*h.pair[1].loc) & loc_bitmask; }
};

template<class Key, class HFCT, size_t TAB_WIDTH, size_t LOC_WIDTH>
struct ExternLocExtractor<Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>, 2>
{
    using Parent_t = Hasher<Key,HFCT,TAB_WIDTH,LOC_WIDTH,2,1>;
    using Hashed_t = typename Parent_t::Hashed_t;

    static constexpr size_t loc_bitmask = (1ull << LOC_WIDTH) - 1;

    static inline size_t loc(Hashed_t h, size_t i) { return h.pair[i].loc; }
};

template<size_t TAB_WIDTH, size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
union ExternHashed
{
    uint64_t hash[NWORDS];
    struct Pair
    {
        uint64_t tab : TAB_WIDTH;
        uint64_t loc : LOC_WIDTH;
    };
    Pair pair[NSTUFF];
};

template<size_t LOC_WIDTH, size_t NSTUFF, size_t NWORDS>
union ExternHashed<0, LOC_WIDTH, NSTUFF, NWORDS>
{
    uint64_t hash[NWORDS];
    struct Pair
    {
        uint64_t loc : LOC_WIDTH;
    };
    Pair pair[NSTUFF];
};
