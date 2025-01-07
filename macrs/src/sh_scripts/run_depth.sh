#!/bin/sh

cd $(dirname $0)

#variable with output file name
output_file="results/depth_results.txt";

#set value to increase each test by
increment=2;

#number of different values
num_values=5;

#number of iterations per value
num_iterations=100;

# create file
printf "Depth,Returned Distance\n" > $output_file

for i in $(seq 0 $num_values);
do 
    for j in $(seq 1 $num_iterations);
    do 
        #come up with new value based on increment amount
        test_value=$((i * increment));
        #run the test
        result=$(./macrs heuristic -l tests/depth/test$j-a.txt tests/depth/test$j-b.txt $test_value);
        #add result to csv 
        printf "$test_value" >> "$output_file";
        printf "," >> "$output_file";
        printf "$result" >> "$output_file";

        #check whether we need to print the last new line character (don't print if at last iteration)
        if [[ "$i" != "$num_values" ]] || [[ "$j" != "$num_iterations" ]];
        then
            printf "\n" >> "$output_file";
        fi
    done
done