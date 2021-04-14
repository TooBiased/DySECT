#!/bin/bash

n=1000000
steps=8192

binfolder="bin"
outfolder="out"

for tab in $binfolder/sing_multi_*
do
    for B in 12 8 4
    do
        for H in 4 3 2
        do
            pname="$(basename -- $tab)"
            ./$tab -n $n -cap 50000 -steps $steps -tl 256 -bs $B -nh $H -bfs \
                   -out ${outfolder}/${pname} &
        done
    done
done

for job in `jobs -p`
do
    wait $job
done
