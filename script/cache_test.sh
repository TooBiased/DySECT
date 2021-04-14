#!/bin/bash

n=7000000
cap=10000000
steps=8192

binfolder="bin"
outfolder="out"

for tab in $binfolder/displ_multi_*
do
    for B in 8
    do
        for H in 3
        do
            ./$tab -n $n -cap $cap -steps $steps -tl 256 -bs $B -nh $H -bfs \
                   -out ${outfolder}/${tab}
        done
    done
done

for tab in $binfolder/displ_hop_*
do
    ./$tab -n $n -cap $cap -steps $steps -ns 64 \
           -out ${outfolder}/${tab}
done

for tab in $binfolder/displ_triv_*
do
    ./$tab -n $n -cap $cap -steps $steps \
           -out ${outfolder}/${tab}
done
