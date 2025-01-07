#!/bin/sh

cd $(dirname $0)

#set value to increase each test by
#increment=10;

#number of different values
#num_values=10;

#array of test values
test_values=(50 70 100 200 300)

#number of iterations per value
num_iterations=100;


for i in ${test_values[@]};
do 
    for j in $(seq 1 $num_iterations);
    do 
        num_ret=$((i / 10));
        combined_seq=$(./macrs gen -e $i $num_ret 50);
        read -r seq1 seq2 <<< $combined_seq;
        echo "$seq1" > tests/leaves/leaves_$i/test$j-a.txt;
        echo "$seq2" > tests/leaves/leaves_$i/test$j-b.txt;
    done 
done

