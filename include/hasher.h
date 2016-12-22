#pragma once



static constexpr size_t ct_log(size_t k)
{ return (k-1) ? 1+ct_log(k>>1) : 0; }

template<size_t t_w, bool dpair>
struct Split;

template <class Hashed, size_t tab_width, bool dpair, bool lcomb>
class ExternExtractor;




/* MAIN CLASS *****************************************************************/
template <class Key, class HFct, size_t TAB_WIDTH, size_t NH,
          bool DPAIR = true, bool LCOMB = true>
class Hasher
{
private:
    using This_t    = Hasher<Key, HFct, TAB_WIDTH, NH, DPAIR, LCOMB>;
    using HashFct_t = HFct;
    static constexpr size_t tab_w   = TAB_WIDTH;
    static constexpr size_t pair_w  = (DPAIR) ? 32 : 64;
    static constexpr size_t loc_w   = pair_w - tab_w;
    static constexpr size_t n_hpair = (LCOMB) ? 2 : NH;
    static constexpr size_t n_hval  = NH;
    static constexpr size_t n_hfct  = (DPAIR) ? (n_hpair>>1)+(1&n_hpair)
                                              :  n_hpair;
    static_assert( tab_w < pair_w,
                   "TAB_WIDTH has to be smaller than PAIR_WIDTH.");

    using Split_t = Split<tab_w, DPAIR>;

public:
    union Hashed_t
    {
        uint64_t hash [n_hfct];
        Split_t  split[n_hfct];
    };

    static_assert (sizeof(Hashed_t) == n_hfct*8,
                   "HashSplitter Bigger Than Expected");

    using Extractor_t = ExternExtractor<Hashed_t, tab_w, DPAIR, LCOMB>;

    /* hasher (itself) ********************************************************/
private:
    HashFct_t fct[n_hfct];

public:

    Hasher()
    {
        for (size_t i = 0; i < n_hfct; ++i)
        {
            fct[i] = HashFct_t(  2345745572344267838ull +
                               i*8768656543548765336ull);
        }
    }

    Hashed_t operator()(Key k) const
    {
        Hashed_t result;
        for (size_t i = 0; i < n_hfct; ++i)
        {
            result.hash[i] = fct[i](k);
        }
        return result;
    }
};




/* Hash Splitter Specializations **********************************************/
template<size_t t_w>
struct Split<t_w, true>
{
    uint64_t tab0 : t_w;
    uint64_t loc0 : 32-t_w;
    uint64_t tab1 : t_w;
    uint64_t loc1 : 32-t_w;
};

template<size_t t_w>
struct Split<t_w, false>
{
    uint64_t tab  : t_w;
    uint64_t loc  : 64-t_w;
};

template<>
struct Split<0, true>
{
    static constexpr uint64_t tab0 = 0;
    static constexpr uint64_t tab1 = 0;
    uint64_t loc0 : 32;
    uint64_t loc1 : 32;
};

template<>
struct Split<0, false>
{
    static constexpr uint64_t tab = 0;
    uint64_t loc  : 64;
};

/* Extractor Specializations without LCOMB ************************************/
template <class Hashed, size_t tab_width>
class ExternExtractor<Hashed, tab_width, false, false>
{
public:
    inline static size_t tab(const Hashed& h, size_t i)
    { return h.split[i].tab; }
    inline static size_t loc(const Hashed& h, size_t i)
    { return h.split[i].loc; }
};

template <class Hashed, size_t tab_width>
class ExternExtractor<Hashed, tab_width, true, false>
{
public:
    inline static size_t tab(const Hashed& h, size_t i)
    { return (i & 1) ? h.split[i>>1].tab1
                     : h.split[i>>1].tab0; }
    inline static size_t loc(const Hashed& h, size_t i)
    { return (i & 1) ? h.split[i>>1].loc1
                     : h.split[i>>1].loc0; }
};

/* Extractor Specializations with LCOMB ***************************************/
template <class Hashed, size_t tab_width>
class ExternExtractor<Hashed, tab_width, false, true>
{
private:
    static_assert((tab_width != 0) && (tab_width != 64),
                  "Illegal TAB_WIDTH value 0 or 64.");
    static constexpr size_t tab_mask = (1ull<<tab_width     )-1;
    static constexpr size_t loc_mask = (1ull<<(64-tab_width))-1;
public:
    inline static size_t tab(const Hashed& h, size_t i)
    { return (h.split[0].tab + i*h.split[1].tab) & tab_mask; }
    inline static size_t loc(const Hashed& h, size_t i)
    { return (h.split[0].loc + i*h.split[1].loc) & tab_mask; }
};

template <class Hashed, size_t tab_width>
class ExternExtractor<Hashed, tab_width, true, true>
{
private:
    static constexpr size_t tab_mask = (1ull<<tab_width     )-1;
    static constexpr size_t loc_mask = (1ull<<(32-tab_width))-1;
public:
    inline static size_t tab(const Hashed& h, size_t i)
    { return (h.split[0].tab0 + i*h.split[0].tab1) & tab_mask; }
    inline static size_t loc(const Hashed& h, size_t i)
    { return (h.split[0].loc0 + i*h.split[0].loc1) & loc_mask; }
};
