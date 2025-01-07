#!/bin/sh

cd $(dirname $0)

#variable with output file name
output_file="results/leaves_results.txt";

#set value to increase each test by
#increment=10;

#variable for depth
depth=4;

#number of different values
#num_values=10;

#array of test values
test_values=(50 70 100 200 300)

#number of iterations per value
num_iterations=100;

# create file
printf "Number leaves,Returned Distance\n" > $output_file

for i in ${test_values[@]};
do 
    for j in $(seq 1 $num_iterations);
    do 
        #run the test
        result=$(./macrs heuristic -l tests/leaves/leaves_$i/test$j-a.txt tests/leaves/leaves_$i/test$j-b.txt $depth);
        #add result to csv 
        printf "$i" >> "$output_file";
        printf "," >> "$output_file";
        printf "$result" >> "$output_file";

        #check whether we need to print the last new line character (don't print if at last iteration)
        if [[ "$i" != "300" ]] || [[ "$j" != "$num_iterations" ]];
        then
            printf "\n" >> "$output_file";
        fi
    done
done