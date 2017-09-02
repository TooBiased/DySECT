* DySECT

** Motivation
In many circumstances hash tables can only be space efficient when
they can also adapt to the number of inserted elements.  In many cases
programmers have to correctly estimate the final table size, to create
densely filled hash tables.  Guessing too conservatively will create
sparser tables and guessing to optimistically will create slower
tables (too full) or it will make the table grow.

There has been a lot of research in the area of space efficient hash
tables.  But dynamically growing these space efficient tables has not
received the same attention.  The conventional wisdom still seems to
be full table migration (into a newly allocated larger table).  This
technique violates space constraints even when growing in small steps
since both tables are allocated at the same time.

** Contents
We present to different solutions to this problem.

*** Inplace implementations of common hashing techniques
We have a multitude of implementations of common hashing techniques
(linear probing, robin hood hashing, cuckoo hashing, hopscotch hashing
...).  For each of these techniques, we implement a variant with a fast
cache efficient migration.  For each of these tables, we also implement
a variant that uses memory overallocation to increase the memory
inplace and do an inplace migration.

*** DySECT
The memory overallocation tables still have some problems due to the
linear work during each full table migration.  Additionally,
overallocation is not possible to use in many projects.  Therefore, we
need a space efficient table that works without it.

** Implementation
The goal for our implementation was to stay as close to possible to
the interface established by ~std::unordered_map~.  We are still
working on achieving this goal, so the actual interface might still go
through some minor changes.

** Installation / Usage
Our implementations are all header only, so including the correct file
should be enough.

To try our test files, you should do the following

#+BEGIN_SRC comment
mkdir build
cd build
cmake ..
make
#+END_SRC

This will create a multitude of folders with different tests, each
built with all possible hash functions.
