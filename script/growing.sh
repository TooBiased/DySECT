#!/bin/bash

pre=8000000
n=10000000
steps=8192

binfolder="bin"
outfolder="out"

for load in 0.8 0.85 0.9 0.925 0.95 0.9625 0.975 0.98
do
    for tab in $binfolder/time_multi_*
    do
        for B in 12 8 4
        do
            for H in 4 3 2
            do
                ./$tab -n $n -pre $pre -cap 50000 -steps $steps -load $load -tl 256 -bs $B -nh $H -bfs \
                       -out ${outfolder}/${tab}
                ./$tab -n $n -pre $pre            -steps $steps -load $load -tl 256 -bs $B -nh $H -bfs \
                       -out ${outfolder}/${tab}
            done
        done
    done

    # for tab in $binfolder/time_hop_*
    # do
    # ./$tab -n $n -pre $pre -cap 50000 -steps $steps -load $load -ns 64 \
        #        -out ${outfolder}/${tab}
    # ./$tab -n $n -pre $pre            -steps $steps -load $load -ns 64 \
        #        -out ${outfolder}/${tab}
    # done

    for tab in $binfolder/time_triv_*
    do
        ./$tab -n $n -pre $pre -cap 50000 -steps $steps -load $load \
               -out ${outfolder}/${tab}
        ./$tab -n $n -pre $pre            -steps $steps -load $load \
               -out ${outfolder}/${tab}
    done
done
