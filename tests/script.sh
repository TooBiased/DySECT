#!/bin/bash

cp ../build/step_all .
cp ../build/in/* .
cp ../build/time/* .

n=100000000
nt=30000000
nh=50000000
ns=70000000
steps=512

echo "step_all"
./step_all -bfs -n 50000000 -cap 10000000 -out 11_11/steps_all

echo "grows"
for dis in bfs rwalk rwalkcyc
do
   for tl in 64 256 2048
   do
        for bs in 4 6 8 16
        do
            for cap in 1 $nh $n
            do
                for alpha in 1.10 1.15 1.2 1.25
                do
                    echo "TL ${tl} BS ${bs} al ${alpha}"
                    ./in_grows -n ${nt} -pre ${ns} -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                               -out 11_11/${dis}_TL${tl}_BS${bs}_al${alpha}_n${nt}_cap${cap}
                    ./time_grows -n $n -pre 0 -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                               -out 11_11/grows_${dis}_TL${tl}_BS${bs}_al${alpha}_n${n}_cap${cap}
                done
            done
        done
    done
done

for tab in cuckoo hom2lvl
do
    echo "type: $tab"
    for dis in bfs
    do
        for tl in 256
        do
            for bs in 8
            do
                for cap in $n
                do
                    for alpha in 1.1 1.15 1.2
                    do
                        ./in_${tab}   -n ${nt} -pre ${ns} -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                                      -out 11_11/${tab}_${dis}_TL${tl}_BS${bs}_al${alpha}_n${nt}_cap${cap}
                        ./time_${tab} -n $n -pre 0 -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                                      -out 11_11/${tab}_${dis}_TL${tl}_BS${bs}_al${alpha}_n${n}_cap${cap}
                    done
                done
            done
        done
    done
done

echo "linprob"
./in_linprob   -n ${nt} -pre ${ns} -cap $n   -steps $steps -bfs \
               -out 11_11/linprob_n${nt}_cap${n}
./time_linprob -n $n    -pre 0     -cap $n   -steps $steps -bfs \
               -out 11_11/linprob_n${n}_cap${n}
