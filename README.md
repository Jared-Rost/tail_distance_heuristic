## NSERC Research Award Project

As part of the NSERC award at the University of Manitoba I worked on a Bioinfomatics Project where I developed a Heuristic for calculating the distance between Phylogenetic Networks under Olivier Tremblay Savard.

I made the src/network/heuristic.rs file, some .sh scripts to test the code, and made some modifications to the Main file to allow the heuristic version to be run from the command line; all other work was done by Kaari Landry and Olivier Tremblay Savard.

## Heuristic Algorithm

Phylogenetic networks are better suited than trees to represent complex evolutionary relationships.

Over time research has focused on solving problems related to phylogenetic networks like construction, and tree containment.

Cherry-picking sequences (CPSs) are a tool for solving such problems.

The tail distance, proposed by Landry et al., is based on the structure of the networks’ CPSs and seeks to maximize the common suffix.

CPSs are made of operations that reduce a network by removing 1 leaf of a simple cherry or 1 edge of a reticulated cherry until the network is fully reduced.

This heuristic estimates that cherry picking distance by choosing the option that leads to the most common cherries:

1. Makes a list of common cherries vs unique cherries in each network.

2. Finds the unique/common cherry that, if picked, results in the most common cherries in subnetworks. Depth determines how many cycles ahead it looks.

3. Removes a unique cherry from the larger network, if none are available then it removes a common cherry.

### Optimizations

1. Using a hashmap to store all common/unique cherries associated with a pair of networks

2. Using a hashmap to store the cherry that leads to the most common cherries in subnetworks if picked associated with a pair of networks

3. Branch and Bound for the recursive algorithm, ending the recursion early if on current path we find more unique cherries than the best answer found so far (with the idea being a better answer would result in more common cherries and less unique cherries than the previous best answer).

4. Updating cherries instead of recalculating them when picking cherries from the network

## How To Run

Compile the code:

```
cargo build
```

Generate input networks in newick using:

```
./macrs gen -e {num leaves} {num reticulations} {desired cherry picking distance}
```

(or use one of the gen scripts in src/sh_scripts)

Run the code using the generated networks:

```
./macrs heuristic -l {first input network in newick}.txt {second input network in newick}.txt {depth, essentially how far }
```

(or use one of the run scripts in src/sh scripts)
