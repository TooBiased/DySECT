#!/bin/bash

steps=8192

binfolder="bin"
outfolder="out"

for neps in 16777216 20971520 25165824 29360128
do
    for tab in $binfolder/eps_multi_*
    do
        for B in 12 8 4
        do
            for H in 4 3 2
            do
                pname="$(basename -- $tab)"
                ./$tab -n $neps -steps $steps -tl 256 -bs $B -nh $H -bfs -load 0.8 \
                      -out ${outfolder}/${pname}
            done
        done
    done

    # for tab in $binfolder/eps_hop_*
    # do
    #     pname="$(basename -- $tab)"
    #     ./$tab -n $neps -steps $steps -ns 64 -load 0.8 \
    #                   -out ${outfolder}/${pname}
    # done

    for tab in $binfolder/eps_triv_*
    do
        pname="$(basename -- $tab)"
        ./$tab -n $neps -steps $steps -load 0.8 \
                      -out ${outfolder}/${pname}
    done
done
