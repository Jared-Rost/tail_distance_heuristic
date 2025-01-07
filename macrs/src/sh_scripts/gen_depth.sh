#!/bin/sh

cd $(dirname $0)


#number of iterations per value
num_iterations=100;

    for j in $(seq 1 $num_iterations);
    do 
        combined_seq=$(./macrs gen -e 150 20 50);
        read -r seq1 seq2 <<< $combined_seq;
        echo "$seq1" > tests/depth/test$j-a.txt;
        echo "$seq2" > tests/depth/test$j-b.txt;
    done 

