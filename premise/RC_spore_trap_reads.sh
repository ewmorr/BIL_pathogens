#!/bin/bash

cd FEDRR_all_2024

for i in reads/*_ST_*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}
    
    mv $dir/$r1File $dir/${r1File}.tmp && mv $dir/$r2File $dir/$r1File && mv $dir/${r1File}.tmp $dir/$r2File
)
done
