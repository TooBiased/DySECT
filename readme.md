# DySECT

## Motivation In many circumstances hash tables can only be space
efficient when they can also adapt to the number of inserted elements.
Otherwise programmers have to correctly estimate the final table size
to create densely filled hash tables.  Guessing conservatively will
create sparser tables and guessing optimistically will create slower
tables (too full) or it might even lose elements.

There has been a lot of research in the area of space efficient hash
tables.  But dynamically growing these space efficient tables has not
received the same attention.  The conventional wisdom still seems to
be full table migration (into a newly allocated larger table).  This
technique violates space constraints even when growing in small steps
because both the source and the target table are allocated at the same
time.

## Contents

### DySECT
This is our clean and simple data-structure presented at [ESA 2017]
(http://drops.dagstuhl.de/opus/volltexte/2017/7848/pdf/LIPIcs-ESA-2017-58.pdf).
It is based on a number of subtables that can grow independently.
Elements can be moved between subtatbles -- using techniques similar
to cuckoo hashing -- such that the free space generated from growing
one subtable can be used efficiently.

### Common Hashing Techniques
We also offer a number implementations of common hashing techniques
(cuckoo hashing, linear probing, robin hood hashing, and hopscotch
hashing).  For every one of them we implement a cache efficient
migration technique that is based on the premise that elements are
basically sorted by their hash value (hash value in [0,1) is used as
scaling factor).

These cannot truely be space efficient, because during a growth step,
both the new and the old table are allocated.  But outside of the
migration they are.  This necessitates small growing factors and thus
is relatively inefficient.

### Inplace migration (through overallocation)
This technique allows us to increase the size of the hash table in
place.  Therefore removing the necessity of having both an old and a
new table during the migration.  Instead the whole table is reordered
in place.  To make this work, we (ab)use the way virtual memory works.
Instead of allocating a piece of memory that has the same size as the
initial size of the hash table, we allocate a large chunk of memory
(maximum final size).  This memory will be purely virtual until it is
accessed.  Therefore, only the beginning part (where we build the hash
table) is actually mapped to physical memory.  Whenever the table is
grown, we access more of the virtual memory, therefore, resizing the
table in place.

## Implementation
The goal for our implementation was to stay as close
to possible to the interface established by `std::unordered_map`.  We
are still working on achieving this goal, so the actual interface
might still go through some minor changes.  *Take a look at the
`example` folder!*

```cpp
// our dysect data-structure
dysect::cuckoo_dysect

// common hashing techniques
dysect::cuckoo_standard
dysect::prob_linear
dysect::prob_robin
dysect::prob_hopscotch

// in place variants
dysect::cuckoo_dysect_inplace // uses virtual memory trick for subtable migration
dysect::cuckoo_standard_inplace
dysect::prob_linear_inplace
dysect::prob_robin_inplace
dysect::prob_hopscotch_inplace

// multitable variants of common techniques
dysect::cuckoo_independent_2lvl
dysect::multitable_linear
dysect::multitable_robin

// experimental stuff
dysect::cuckoo_deamortized
dysect::cuckoo_overlap
dysect::cuckoo_overlap_inplace
dysect::prob_linear_doubling

```

## Installation / Usage
Our implementations are all header only, so including the correct file
should be enough.

Try our test files by:

```sh
mkdir build
cd build
cmake ..
make
```

This will create a multitude of folders with different tests, each
built with many of our hashing techniques. Use ccmake, to change
parameters like the hash function, and virtual memory size (for in place
variants).
