#!/bin/bash

cp ../build/step_all .
cp ../build/in/* .
cp ../build/time/* .

n=10000000
nt=3000000
nh=5000000
ns=7000000
steps=512

echo "Step_all"
# ./step_all -bfs -n 50000000 -cap 10000000 -out 17_11/steps_all

echo "Grows"
for dis in bfs rwalk rwalkcyc
do
   for tl in 64 128 256 512 1024 2048
   do
        for bs in 4 6 8 12 16
        do
            for cap in 1 $nh $n
            do
                for alpha in 1.05 1.1 1.15 1.2 1.25
                do
                    ./time_grows -n $n    -pre 0     -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                                 -out 17_11/grows_${dis}
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
                    for alpha in 1.05 1.1 1.15 1.2 1.25
                    do
                        ./time_${tab} -n $n -pre 0 -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                                      -out 17_11/${tab}_${dis}
                    done
                done
            done
        done
    done
done



echo "Linprob"
for cap in 1 $nh $n
do
    ./time_linprob -n $n -pre 0 -cap $cap -steps $steps -bfs \
               -out 17_11/linprob
done

echo "Spaceprob"
for alpha in 1.1 1.15 1.2 1.25
do
    ./time_spaceprob -n $n  -pre 0  -cap $n  -alpha $alpha  -steps $steps -bfs \
                   -out 17_11/spaceprob
done

echo "Hopscotch"
for cap in 1 $nh $n
do
    for alpha in 1.1 1.15 1.2 1.25
    do
        ./time_hopscotch -n $n  -pre 0  -cap $cap  -alpha $alpha  -steps $steps -bfs \
                   -out 17_11/hopscotch
    done
done


### VARIABLE NR. OF ELEMENTS ############################################

echo "Variable Nr"
for t in 1000000 2000000 4000000 6000000 8000000 10000000 12500000 15000000 17500000 20000000 #1000000 # 10000000 25000000 50000000 75000000 100000000 1250000000 150000000
do
    for alpha in 1.1 1.15 1.2 1.25
    do
        ./time_grows     -bfs -n $t -pre 0 -cap 1 -steps $steps -alpha $alpha \
                         -out 17_11/grows_incr_bfs
        ./time_hopscotch -bfs -n $t -pre 0 -cap 1 -steps $steps -alpha $alpha \
                         -out 17_11/hopscotch_incr_bfs
    done
    ./time_linprob -bfs -n $t -pre 0 -cap 1 -steps $steps \
                   -out 17_11/linprob_incr_bfs
done

### DISPLACEMENT MEASUREMENTS ##################################################

echo "Displacement Measurements"
echo "Grows"
for dis in bfs rwalk rwalkcyc
do
   for tl in 64 128 256 512 1024 2048
   do
        for bs in 4 6 8 12 16
        do
            for cap in 1 $nh $n
            do
                for alpha in 1.05 1.1 1.15 1.2 1.25
                do
                    echo "TL ${tl} BS ${bs} al ${alpha}"
                    ./in_grows   -n ${nt} -pre ${ns} -cap $cap -steps $steps -alpha $alpha -tl $tl -bs $bs -${dis} \
                                 -out 17_11/grows_${dis}_TL${tl}_BS${bs}_al${alpha}_n${nt}_cap${cap}
                done
            done
        done
   done
done
