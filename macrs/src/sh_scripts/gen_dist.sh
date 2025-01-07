#!/bin/sh

#chmod 755 run_leaves.sh

cd $(dirname $0)

#set value to increase each test by
#increment=10;

#number of different values
#num_values=8;

#array of test values
test_values=(10 20 50 100 200)

#number of iterations per value
num_iterations=100;


for i in ${test_values[@]};
do 
    for j in $(seq 1 $num_iterations);
    do 
        combined_seq=$(./macrs gen -e 200 20 $i);
        read -r seq1 seq2 <<< $combined_seq;
        echo "$seq1" > tests/distance/dist_$i/test$j-a.txt;
        echo "$seq2" > tests/distance/dist_$i/test$j-b.txt;
    done 
done

