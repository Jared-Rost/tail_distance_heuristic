#!/bin/sh

cd $(dirname $0)

#variable with output file name
output_file="results/ret_results.txt";

#set value to increase each test by
#increment=3;

#variable for depth
depth=4;

#number of different values
#num_values=8;

#array of test values
test_values=(0 10 20 30 50)

#number of iterations per value
num_iterations=100;

# create file
printf "Reticulations,Returned Distance\n" > $output_file

for i in ${test_values[@]};
do 
    for j in $(seq 1 $num_iterations);
    do 
        #run the test
        result=$(./macrs heuristic -l tests/ret/ret_$i/test$j-a.txt tests/ret/ret_$i/test$j-a.txt $depth);
        #add result to csv 
        printf "$i" >> "$output_file";
        printf "," >> "$output_file";
        printf "$result" >> "$output_file";

        #check whether we need to print the last new line character (don't print if at last iteration)
        if [[ "$i" != "50" ]] || [[ "$j" != "$num_iterations" ]];
        then
            printf "\n" >> "$output_file";
        fi
    done
done