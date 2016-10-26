#!/bin/bash

cp ../build/in .

n=1000000
steps=1024



for dis in bfs rwalk rwalkcyc
do
   for tl in 64 128 256
   do
        for bs in 4 8
        do
            for cap in 1 $n
            do
                for alpha in 1.10 1.12 1.15
                do
                    echo "Parameters: -n $n -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} -out out/${dis}_TL${tl}BS${bs}_al${alpha}_cap${cap}"
                    ./in -n $n -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} -out out/${dis}_TL${tl}BS${bs}_al${alpha}_cap${cap}
                done
            done
        done
    done
done

