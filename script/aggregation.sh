#!/bin/bash

steps=8192

binfolder="bin"
outfolder="out"
infile="input"

for load in 0.8 0.85 0.9 0.925 0.95 0.9625 0.975 0.98
do
    for tab in $binfolder/crawl_multi_*
    do
        for B in 8
        do
            for H in 3
            do
                pname="$(basename -- $tab)"
                ./$tab -cap 50000 -steps $steps -load $load -tl 256 -bs $B -nh $H -bfs \
                       -in $ifile -out ${outfolder}/${pname}
            done
        done
    done

    # for tab in $binfolder/crawl_hop_*
    # do
            # pname="$(basename -- $tab)"
            # ./$tab -cap 50000 -steps $steps -load $load -ns 64 \
            #            -in $ifile -out ${outfolder}/${pname}
    # done

    for tab in $binfolder/crawl_triv_*
    do
        pname="$(basename -- $tab)"
        ./$tab -cap 50000 -steps $steps -load $load \
                       -in $ifile -out ${outfolder}/${pname}
    done
done
