#ifndef COUNTING_WAIT_H
#define COUNTING_WAIT_H


#include <stdlib.h>
#include <atomic>
#include <sys/time.h>
#include <memory>
#include <iostream>

#include <linux/futex.h>
#include <unistd.h>
#include <sys/syscall.h>

static long sys_futex(void *addr1, int op, int val1, struct timespec *timeout, void *addr2, int val3)
{
    return syscall(SYS_futex, addr1, op, val1, timeout, addr2, val3);
}

class alignas(64) counting_wait
{
public:
    inline counting_wait(int start = 0) : counter(start)
    {
        if (sizeof(std::atomic_int) != 4)
            std::cout << "std::atomic_int has wrong size:"
                      << sizeof(std::atomic_int)<< std::endl;
    }

    inline bool inc_if(int exp)
    {
        auto temp = exp;
        return counter.compare_exchange_strong(temp, exp+1,
                                               std::memory_order_acq_rel,
                                               std::memory_order_acquire);
    }

    inline bool wait_if(int exp)
    {
        //while (counter.load(std::memory_order_acquire) < l_epoch) ;  //temporary should soon be removed
        auto ecode = sys_futex(&counter, FUTEX_WAIT, exp, NULL, NULL, 0);
        return !ecode;
    }

    inline uint wake(uint n_threads = 9999)
    {
        return sys_futex(&counter, FUTEX_WAKE, n_threads, NULL, NULL, 0);
    }

private:
    std::atomic_int counter;
};


#endif // COUNTING_WAIT_H
