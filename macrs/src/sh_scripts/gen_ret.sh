#!/bin/sh

cd $(dirname $0)

#set value to increase each test by
#increment=3;

#number of different values
#num_values=8;

#array of test values
test_values=(0 10 20 30 50)

#number of iterations per value
num_iterations=100;


for i in ${test_values[@]};
do 
    for j in $(seq 1 $num_iterations);
    do 
        combined_seq=$(./macrs gen -e 200 $i 50);
        read -r seq1 seq2 <<< $combined_seq;
        echo "$seq1" > tests/ret/ret_$i/test$j-a.txt;
        echo "$seq2" > tests/ret/ret_$i/test$j-b.txt;
    done 
done

