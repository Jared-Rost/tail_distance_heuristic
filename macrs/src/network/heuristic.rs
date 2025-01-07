/**
* heuristic
*
* COMP3350 SECTION A02
*
* Jared Rost
* August 21, 2024
*
* PURPOSE: Given two input networks, to return two cherry picking sequences that attempt to maximize the common tail by picking common cherries first.
*  
**/
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;

use crate::network::cherry::get_reducible_cherries;
use crate::network::cherry::reduce_cherry;
use crate::network::cherry::reduce_cherry_multi_label;
use crate::network::cherry::update_reducible_cherries;
use crate::network::cherry::Cherry;
use crate::network::cherry::CherryType;
use crate::network::cherry::CherryType::Retic;
use crate::network::cherry::CherryType::Simple;
use crate::network::Network;
use crate::network::Node;
use std::time::{Duration, Instant};

/**
 * Struct CherryInfo
 *
 * Used to store all info related to a cherry (cherry object + labels associated with it) for the purpose of faster comparisons using BtreeSets.  
 * Single labels only.
 */
#[derive(Clone)]
pub struct CherryInfo {
    // labels associated with this cherry
    pub label1: String,
    pub label2: String,
    // cherry object itself
    pub cherry_object: Cherry,
}
impl PartialEq for CherryInfo {
    fn eq(&self, other: &Self) -> bool {
        // since the single label version compares cherry versions ({a,b} and {b,a}) then for them to be considered common the labels must exactly match.
        self.label1 == other.label1
            && self.label2 == other.label2
            && self.cherry_object.cherry_type == other.cherry_object.cherry_type
    }
}
impl Eq for CherryInfo {}
impl PartialOrd for CherryInfo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for CherryInfo {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.label1, &self.label2, self.cherry_object.cherry_type).cmp(&(
            &other.label1,
            &other.label2,
            other.cherry_object.cherry_type,
        ))
    }
}

/**
 * Struct CherryInfoMultiLabel
 *
 * Used to store all info related to a cherry (cherry object + labels associated with it) for the purpose of faster comparisons using BtreeSets.  
 * Multi labels only.
 */
#[derive(Clone)]
pub struct CherryInfoMultiLabel {
    // labels associated with this cherry
    pub label1: BTreeSet<String>,
    pub label2: BTreeSet<String>,
    // cherry object itself
    pub cherry_object: Cherry,
}
impl PartialEq for CherryInfoMultiLabel {
    fn eq(&self, other: &Self) -> bool {
        // since the multi label version only compares 1 version of simple cherries (only {a,b}) then that means we check both possible intersections
        if self.cherry_object.cherry_type == Simple {
            (self.label1.intersection(&other.label1).count() > 0
                && self.label2.intersection(&other.label2).count() > 0
                && self.cherry_object.cherry_type == other.cherry_object.cherry_type)
                || (self.label1.intersection(&other.label2).count() > 0
                    && self.label2.intersection(&other.label1).count() > 0
                    && self.cherry_object.cherry_type == other.cherry_object.cherry_type)
        }
        // since the multi label version compares both versions of reticulated cherries ({a,b} and {b,a}) then that means we check only one intersection
        else {
            self.label1.intersection(&other.label1).count() > 0
                && self.label2.intersection(&other.label2).count() > 0
                && self.cherry_object.cherry_type == other.cherry_object.cherry_type
        }
    }
}
impl Eq for CherryInfoMultiLabel {}
impl PartialOrd for CherryInfoMultiLabel {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for CherryInfoMultiLabel {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.label1, &self.label2, self.cherry_object.cherry_type).cmp(&(
            &other.label1,
            &other.label2,
            other.cherry_object.cherry_type,
        ))
    }
}

/**
 * find_heuristic_sequence
 *
 * Main function for single label CPS, takes two networks and the depth then returns the CPSs for both
 * @param Network<Node<String>> n1 - first network to find CPS for
 * @param Network<Node<String>> n2 - second network to find CPS for
 * @returns Vec<Vec<(String, String, CherryType)>> - Contains CPS for n1 in index 0 and CPS for n2 in index 1
 */
pub fn find_heuristic_sequence(
    mut n1: Network<Node<String>>,
    mut n2: Network<Node<String>>,
    num_iterations: usize,
) -> Vec<Vec<(String, String, CherryType)>> {
    // used for time measurements
    //let mut start = Instant::now();

    // size of networks
    let mut n1_size = n1.get_n();
    let mut n2_size = n2.get_n();

    // special case where we reduce a common cherry in n2 then we want to immediately go reduce the same cherry in n1
    // boolean indicating whether we are in the special case
    let mut reduced_n2_common_cherry: bool = false;
    // initialize variable storing cherry to default value, not intended to be used as default value ever
    let mut n2_common_cherry_target: CherryInfo = CherryInfo {
        label1: "".to_string(),
        label2: "".to_string(),
        cherry_object: Cherry {
            cherry: (usize::MAX, usize::MAX),
            cherry_type: Retic,
        },
    };

    // vector holding return value
    let mut result: Vec<Vec<(String, String, CherryType)>> = Vec::new();

    // create new variables to hold cherry picking sequences
    let mut n1_sequence: Vec<(String, String, CherryType)> = Vec::new();
    let mut n2_sequence: Vec<(String, String, CherryType)> = Vec::new();

    // create new variables to store the current cherries
    let mut n1_reducible_cherries: Vec<Cherry> = get_reducible_cherries(&n1);
    let mut n2_reducible_cherries: Vec<Cherry> = get_reducible_cherries(&n2);

    // make hashmap storing dynamic programming idea #1: storing results of previous recursive searches within the same depth level
    let mut in_cycle_results_hash: HashMap<String, (usize, usize)> = HashMap::new();

    // make hashmap storing dynamic programming idea #2: storing the cherry types for a specific network combination
    let mut cherry_types_hash: HashMap<String, Vec<Vec<Cherry>>> = HashMap::new();

    // create new variables to hold the common and unique cherries
    let mut cherry_types: Vec<Vec<Cherry>> =
        get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);

    // print statements used for debugging
    /*  print_cherrypicking_sequence(&n1, &n1_reducible_cherries);
    print_cherrypicking_sequence(&n2, &n2_reducible_cherries);

    print_cherrypicking_sequence(&n1, &cherry_types[0]);
    print_cherrypicking_sequence(&n2, &cherry_types[1]);
    print_cherrypicking_sequence(&n1, &cherry_types[2]);
    print_cherrypicking_sequence(&n2, &cherry_types[3]); */

    // loop until there is only 1 cherry left in both of them
    while n1_size > 3 || n2_size > 3 {
        if n1_size > n2_size {
            // create variable to store largest found so far
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early
            let selected_cherry; // variable to store the cherry we ultimetely decide on removing
            let mut max_index: usize = 0; // index of the cherry who chose to remove in list of reducible cherries
            let string_representation; // newick representations of the networks used for hashmaps

            // special case to ensure we reduce same common cherry from both
            if reduced_n2_common_cherry {
                // default location to impossible value
                let mut n2_common_cherry_index = usize::MAX;

                // find the location of the selected cherry
                for (i, current_n1_cherry) in n1_reducible_cherries.iter().enumerate() {
                    if n1.get_label(current_n1_cherry.cherry.0) == n2_common_cherry_target.label1
                        && n1.get_label(current_n1_cherry.cherry.1)
                            == n2_common_cherry_target.label2
                        && current_n1_cherry.cherry_type
                            == n2_common_cherry_target.cherry_object.cherry_type
                    {
                        n2_common_cherry_index = i;
                        break;
                    }
                }

                // if we didn't find it then exit - should never happen
                if n2_common_cherry_index == usize::MAX {
                    panic!("Auto remove indexes do not line up");
                }

                // add the chosen cherry to n1 sequence
                n1_sequence.push((
                    n1.get_label(n1_reducible_cherries[n2_common_cherry_index].cherry.0),
                    n1.get_label(n1_reducible_cherries[n2_common_cherry_index].cherry.1),
                    n1_reducible_cherries[n2_common_cherry_index]
                        .cherry_type
                        .clone(),
                ));

                // get the chosen cherry
                selected_cherry = n1_reducible_cherries.swap_remove(n2_common_cherry_index);

                reduced_n2_common_cherry = false;
            }
            // check to see if there are unique cherries available
            else if cherry_types[2].len() > 0 {
                // look ahead at all possible options
                for (i, current_n1_cherry) in cherry_types[2].iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries(
                        n1.clone(), // same networks, the reducing of the cherry is done in the curirsive calls
                        n2.clone(),
                        Some(current_n1_cherry.clone()), // we want to reduce the current n1 cherry so we wrap it in an option
                        None::<Cherry>, // we do not want to reduce any n2 cherries so we
                        n1_reducible_cherries.clone(), // clone reducible cherries so they don't need to be recalculated
                        n2_reducible_cherries.clone(), // clone reducible cherries so they don't need to be recalculated
                        num_iterations,                // depth, given by user
                        (0, 0), // we don't have any common or unique cherries yet
                        most_common_cherries.1, // send it the # of unique cherries in our best branch found so far
                        &mut in_cycle_results_hash, // hash maps for reducing time
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // check to see if this is better than our current best
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n1_sequence.push((
                    n1.get_label(cherry_types[2][max_index].cherry.0),
                    n1.get_label(cherry_types[2][max_index].cherry.1),
                    cherry_types[2][max_index].cherry_type.clone(),
                ));

                // get the chosen cherry
                selected_cherry = cherry_types[2].swap_remove(max_index);
            }
            // otherwise look at common cherries too
            else {
                // look ahead at all possible options
                for (i, current_n1_cherry) in n1_reducible_cherries.iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries(
                        n1.clone(),
                        n2.clone(),
                        Some(current_n1_cherry.clone()),
                        None::<Cherry>,
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n1_sequence.push((
                    n1.get_label(n1_reducible_cherries[max_index].cherry.0),
                    n1.get_label(n1_reducible_cherries[max_index].cherry.1),
                    n1_reducible_cherries[max_index].cherry_type.clone(),
                ));

                // get the chosen cherry
                selected_cherry = n1_reducible_cherries.swap_remove(max_index);
            }

            // reduce the network by the chosen cherry
            reduce_cherry(&mut n1, &selected_cherry);

            // get new reducable cherries for n1
            n1_reducible_cherries =
                update_reducible_cherries(&n1, &selected_cherry, true, n1_reducible_cherries);

            // get hash key
            string_representation = get_hash_key(&n1, &n2);

            // calculate the new common cherries and unique cherries (if not already calculated)
            if cherry_types_hash.contains_key(&string_representation) {
                cherry_types = cherry_types_hash[&string_representation].clone();
            } else {
                cherry_types =
                    get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
                cherry_types_hash.insert(string_representation, cherry_types.clone());
            }

            // recalculate the size of the network now that a cherry has been removed
            n1_size = n1.get_n();
        }
        // otherwise we remove a cherry from n2
        else {
            // create variable to store largest found so far
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early
            let mut max_index: usize = 0;
            let selected_cherry;
            let string_representation;

            // check to see if there are unique cherries available
            if cherry_types[3].len() > 0 {
                // look ahead at all possible options
                for (i, current_n2_cherry) in cherry_types[3].iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries(
                        n1.clone(),
                        n2.clone(),
                        None::<Cherry>,
                        Some(current_n2_cherry.clone()),
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n2_sequence.push((
                    n2.get_label(cherry_types[3][max_index].cherry.0),
                    n2.get_label(cherry_types[3][max_index].cherry.1),
                    cherry_types[3][max_index].cherry_type.clone(),
                ));

                // get the selected cherry
                selected_cherry = cherry_types[3].swap_remove(max_index);
            }
            // otherwise use common cherries
            else {
                // look ahead at all possible options
                for (i, current_n2_cherry) in n2_reducible_cherries.iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries(
                        n1.clone(),
                        n2.clone(),
                        None::<Cherry>,
                        Some(current_n2_cherry.clone()),
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n2_sequence.push((
                    n2.get_label(n2_reducible_cherries[max_index].cherry.0),
                    n2.get_label(n2_reducible_cherries[max_index].cherry.1),
                    n2_reducible_cherries[max_index].cherry_type.clone(),
                ));

                // get the selected cherry
                selected_cherry = n2_reducible_cherries.swap_remove(max_index);

                reduced_n2_common_cherry = true;
                n2_common_cherry_target = CherryInfo {
                    label1: n2.get_label(selected_cherry.cherry.0),
                    label2: n2.get_label(selected_cherry.cherry.1),
                    cherry_object: selected_cherry.clone(),
                };
            }

            // reduce the network by the chosen cherry
            reduce_cherry(&mut n2, &selected_cherry);

            // get new reducable cherries for n2
            n2_reducible_cherries =
                update_reducible_cherries(&n2, &selected_cherry, true, n2_reducible_cherries);

            // get hash key
            string_representation = get_hash_key(&n1, &n2);

            // calculate the new common cherries and unique cherries (if not already calculated)
            if cherry_types_hash.contains_key(&string_representation) {
                cherry_types = cherry_types_hash[&string_representation].clone();
            } else {
                cherry_types =
                    get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
                cherry_types_hash.insert(string_representation, cherry_types.clone());
            }

            // recalculate the size of the network now that a cherry has been removed
            n2_size = n2.get_n();
        }

        // at the end of every cycle clear the hashcap because it doesn't work when the depth is different
        in_cycle_results_hash.clear();
    }

    // manually add the last cherry
    n1_sequence.push((
        n1.get_label(n1_reducible_cherries[0].cherry.0),
        n1.get_label(n1_reducible_cherries[0].cherry.1),
        n1_reducible_cherries[0].cherry_type.clone(),
    ));
    n2_sequence.push((
        n2.get_label(n2_reducible_cherries[0].cherry.0),
        n2.get_label(n2_reducible_cherries[0].cherry.1),
        n2_reducible_cherries[0].cherry_type.clone(),
    ));

    // used for measuring time
    //println!("{:?} time", start.elapsed());

    // add the two vectors and return it
    result.push(n1_sequence);
    result.push(n2_sequence);
    result
}
/**
 * find_heuristic_sequence_multi_label
 *
 * Main function for multi label CPS, takes two networks and the depth then returns the multi label CPSs for both
 * @param Network<Node<String>> n1 - first network to find CPS for
 * @param Network<Node<String>> n2 - second network to find CPS for
 * @returns Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> - Contains CPS for n1 in index 0 and CPS for n2 in index 1
 */
pub fn find_heuristic_sequence_multi_label(
    mut n1: Network<Node<String>>,
    mut n2: Network<Node<String>>,
    num_iterations: usize,
) -> Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> {
    // used for measuring time
    //let mut start = Instant::now();

    // size of networks
    let mut n1_size = n1.get_n();
    let mut n2_size = n2.get_n();

    // vector holding return value
    let mut result: Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> = Vec::new();

    // create new variables to hold cherry picking sequences
    let mut n1_sequence: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> = Vec::new();
    let mut n2_sequence: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> = Vec::new();

    // create new variables to store the current cherries
    let mut n1_reducible_cherries: Vec<Cherry> = get_reducible_cherries_multi_label(&n1);
    let mut n2_reducible_cherries: Vec<Cherry> = get_reducible_cherries_multi_label(&n2);

    // make hashmap storing dynamic programming idea #1: storing results of previous recursive searches within the same depth level
    let mut in_cycle_results_hash: HashMap<String, (usize, usize)> = HashMap::new();

    // make hashmap storing dynamic programming idea #2: storing the cherry types for a specific network combination
    let mut cherry_types_hash: HashMap<String, Vec<Vec<Cherry>>> = HashMap::new();

    // create new variables to hold the common and unique cherries
    let mut cherry_types: Vec<Vec<Cherry>> =
        get_multi_label_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);

    // used for debugging
    /* print_cherrypicking_sequence(&n1, &n1_reducible_cherries);
       print_cherrypicking_sequence(&n2, &n2_reducible_cherries);

       print_cherrypicking_sequence(&n1, &cherry_types[0]);
       print_cherrypicking_sequence(&n2, &cherry_types[1]);
       print_cherrypicking_sequence(&n1, &cherry_types[2]);
       print_cherrypicking_sequence(&n2, &cherry_types[3]);
    */

    // loop until there is only 1 cherry left in both of them
    while n1_size > 3 || n2_size > 3 {
        // if n1 is bigger than n2 reduce a n1 cherry
        if n1_size > n2_size {
            // create variable to store largest found so far
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early
            let selected_cherry; // variable to store the cherry we ultimetely decide on removing
            let mut max_index: usize = 0; // index of the cherry who chose to remove in list of reducible cherries
            let string_representation; // newick representations of the networks used for hashmaps

            // if there are unique cherries available reduce one of them
            if cherry_types[2].len() > 0 {
                // look ahead at all possible options
                for (i, current_n1_cherry) in cherry_types[2].iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries_multi_label(
                        n1.clone(), // same networks, the reducing of the cherry is done in the curirsive calls
                        n2.clone(),
                        Some(current_n1_cherry.clone()), // we want to reduce the current n1 cherry so we wrap it in an option
                        None::<Cherry>, // we do not want to reduce any n2 cherries so we
                        n1_reducible_cherries.clone(), // clone reducible cherries so they don't need to be recalculated
                        n2_reducible_cherries.clone(), // clone reducible cherries so they don't need to be recalculated
                        num_iterations,                // depth, given by user
                        (0, 0), // we don't have any common or unique cherries yet
                        most_common_cherries.1, // send it the # of unique cherries in our best branch found so far
                        &mut in_cycle_results_hash, // hash maps for reducing time
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // check to see if this is better than our current best
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n1_sequence.push((
                    n1.get_labels(cherry_types[2][max_index].cherry.0),
                    n1.get_labels(cherry_types[2][max_index].cherry.1),
                    cherry_types[2][max_index].cherry_type.clone(),
                ));

                // get the chosen cherry
                selected_cherry = cherry_types[2].swap_remove(max_index);
            }
            // otherwise look at common cherries instead
            else {
                // look ahead at all possible options
                for (i, current_n1_cherry) in n1_reducible_cherries.iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries_multi_label(
                        n1.clone(),
                        n2.clone(),
                        Some(current_n1_cherry.clone()),
                        None::<Cherry>,
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n1_sequence.push((
                    n1.get_labels(n1_reducible_cherries[max_index].cherry.0),
                    n1.get_labels(n1_reducible_cherries[max_index].cherry.1),
                    n1_reducible_cherries[max_index].cherry_type.clone(),
                ));

                // get the chosen cherry
                selected_cherry = n1_reducible_cherries.swap_remove(max_index);
            }

            // reduce the network by the chosen cherry
            reduce_cherry_multi_label(&mut n1, &selected_cherry);

            // get new reducable cherries for n1
            n1_reducible_cherries = update_reducible_cherries_multi_label(
                &n1,
                &selected_cherry,
                true,
                n1_reducible_cherries,
            );

            // get hash key
            string_representation = get_hash_key(&n1, &n2);

            // calculate the new common cherries and unique cherries (if not already calculated)
            if cherry_types_hash.contains_key(&string_representation) {
                cherry_types = cherry_types_hash[&string_representation].clone();
            } else {
                cherry_types = get_multi_label_cherry_types(
                    &n1,
                    &n2,
                    &n1_reducible_cherries,
                    &n2_reducible_cherries,
                );
                cherry_types_hash.insert(string_representation, cherry_types.clone());
            }

            // recalculate the size of the network now that a cherry has been removed
            n1_size = n1.get_n();
        } else {
            // create variable to store largest found so far
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early
            let mut max_index: usize = 0;
            let selected_cherry;
            let string_representation;

            // check to see if there are unique cherries available
            if cherry_types[3].len() > 0 {
                // look ahead at all possible options
                for (i, current_n2_cherry) in cherry_types[3].iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries_multi_label(
                        n1.clone(),
                        n2.clone(),
                        None::<Cherry>,
                        Some(current_n2_cherry.clone()),
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // add the chosen cherry
                n2_sequence.push((
                    n2.get_labels(cherry_types[3][max_index].cherry.0),
                    n2.get_labels(cherry_types[3][max_index].cherry.1),
                    cherry_types[3][max_index].cherry_type.clone(),
                ));

                // get the selected cherry
                selected_cherry = cherry_types[3].swap_remove(max_index);
            }
            // otherwise use common cherries
            else {
                // look ahead at all possible options
                for (i, current_n2_cherry) in n2_reducible_cherries.iter().enumerate() {
                    let temp_num_cherry_types = predict_cherries_multi_label(
                        n1.clone(),
                        n2.clone(),
                        None::<Cherry>,
                        Some(current_n2_cherry.clone()),
                        n1_reducible_cherries.clone(),
                        n2_reducible_cherries.clone(),
                        num_iterations,
                        (0, 0),
                        most_common_cherries.1,
                        &mut in_cycle_results_hash,
                        &mut cherry_types_hash,
                    );

                    // first check to see this is the first iteration
                    if i == 0 {
                        // if so then this is the best
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                    // otherwise check to see if this version has more common cherries
                    else if temp_num_cherry_types.0 > most_common_cherries.0 {
                        most_common_cherries = temp_num_cherry_types;
                        max_index = i;
                    }
                }

                // push the selected cherry
                n2_sequence.push((
                    n2.get_labels(n2_reducible_cherries[max_index].cherry.0),
                    n2.get_labels(n2_reducible_cherries[max_index].cherry.1),
                    n2_reducible_cherries[max_index].cherry_type.clone(),
                ));

                // get the selected cherry
                selected_cherry = n2_reducible_cherries.swap_remove(max_index);
            }

            // reduce the network by the chosen cherry
            reduce_cherry_multi_label(&mut n2, &selected_cherry);

            // get new reducable cherries for n2
            n2_reducible_cherries = update_reducible_cherries_multi_label(
                &n2,
                &selected_cherry,
                true,
                n2_reducible_cherries,
            );

            // get hash key
            string_representation = get_hash_key(&n1, &n2);

            // calculate the new common cherries and unique cherries (if not already calculated)
            if cherry_types_hash.contains_key(&string_representation) {
                cherry_types = cherry_types_hash[&string_representation].clone();
            } else {
                cherry_types = get_multi_label_cherry_types(
                    &n1,
                    &n2,
                    &n1_reducible_cherries,
                    &n2_reducible_cherries,
                );
                cherry_types_hash.insert(string_representation, cherry_types.clone());
            }

            // recalculate the size of the network now that a cherry has been removed
            n2_size = n2.get_n();
        }

        // at the end of every cycle clear the hashcap because it doesn't work when the depth is different
        in_cycle_results_hash.clear();
    }

    // manually add the last cherry
    n1_sequence.push((
        n1.get_labels(n1_reducible_cherries[0].cherry.0),
        n1.get_labels(n1_reducible_cherries[0].cherry.1),
        n1_reducible_cherries[0].cherry_type.clone(),
    ));
    n2_sequence.push((
        n2.get_labels(n2_reducible_cherries[0].cherry.0),
        n2.get_labels(n2_reducible_cherries[0].cherry.1),
        n2_reducible_cherries[0].cherry_type.clone(),
    ));

    // used to measure time
    //println!("{:?} time", start.elapsed());

    // combine the two vectors together into 1 vector and return it
    result.push(n1_sequence);
    result.push(n2_sequence);
    result
}
/**
 * get_cherry_types
 *
 * Given two networks and the list of available cherries for each, returns a lists of common cherries and unique cherries in each (single label version)
 * @param &Network<Node<String>> n1 - input network to get labels from
 * @param &Network<Node<String>> n2 - input network to get labels from
 * @returns Vec<Vec<Cherry>> - Contains n1 common cherries in [0], n2 common cherries in [1], n1 unique cherries in [2], and n2 unique cherries in [3]
 */
fn get_cherry_types(
    n1: &Network<Node<String>>,
    n2: &Network<Node<String>>,
    n1_cherries: &Vec<Cherry>,
    n2_cherries: &Vec<Cherry>,
) -> Vec<Vec<Cherry>> {
    // let mut start = Instant::now();

    // declare variable to store the common cherries
    let mut types_of_cherries: Vec<Vec<Cherry>> = Vec::new();
    let mut n1_common_cherries_tree: BTreeSet<CherryInfo> = BTreeSet::new();
    let mut n2_common_cherries_tree: BTreeSet<CherryInfo> = BTreeSet::new();
    let mut n1_unique_cherries_tree: BTreeSet<CherryInfo> = BTreeSet::new();

    let n1_common_cherries: Vec<Cherry>;
    let n2_common_cherries: Vec<Cherry>;
    let n1_unique_cherries: Vec<Cherry>;
    let n2_unique_cherries: Vec<Cherry>;

    // step 0, make new btreesets
    let mut s1: BTreeSet<CherryInfo> = BTreeSet::new();
    let mut s2: BTreeSet<CherryInfo> = BTreeSet::new();

    // step 1, get all labels associated with each cherry and put them into a new object
    // step 2, put that object into btreeset

    for current_n1_cherry in n1_cherries.into_iter() {
        s1.insert(CherryInfo {
            label1: n1.get_label(current_n1_cherry.cherry.0),
            label2: n1.get_label(current_n1_cherry.cherry.1),
            cherry_object: current_n1_cherry.clone(),
        });
    }

    for current_n2_cherry in n2_cherries.into_iter() {
        s2.insert(CherryInfo {
            label1: n2.get_label(current_n2_cherry.cherry.0),
            label2: n2.get_label(current_n2_cherry.cherry.1),
            cherry_object: current_n2_cherry.clone(),
        });
    }

    // step 3, use various set sorting methods to organize it appropriately

    for x in s1.into_iter() {
        if s2.contains(&x) {
            n2_common_cherries_tree.insert(s2.take(&x).unwrap());
            n1_common_cherries_tree.insert(x);
        } else {
            n1_unique_cherries_tree.insert(x);
        }
    }

    n1_common_cherries = n1_common_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n2_common_cherries = n2_common_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n1_unique_cherries = n1_unique_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n2_unique_cherries = s2.into_iter().map(|p| p.cherry_object).collect();

    types_of_cherries.push(n1_common_cherries);
    types_of_cherries.push(n2_common_cherries);
    types_of_cherries.push(n1_unique_cherries);
    types_of_cherries.push(n2_unique_cherries);

    //  println!("{:?} time", start.elapsed());

    types_of_cherries
}
/**
 * get_multi_label_cherry_types
 *
 * Given two networks and the list of available cherries for each, returns a lists of common cherries and unique cherries in each (multi label version)
 * @param &Network<Node<String>> n1 - input network to get labels from
 * @param &Network<Node<String>> n2 - input network to get labels from
 * @returns Vec<Vec<Cherry>> - Contains n1 common cherries in [0], n2 common cherries in [1], n1 unique cherries in [2], and n2 unique cherries in [3]
 */
fn get_multi_label_cherry_types(
    n1: &Network<Node<String>>,
    n2: &Network<Node<String>>,
    n1_cherries: &Vec<Cherry>,
    n2_cherries: &Vec<Cherry>,
) -> Vec<Vec<Cherry>> {
    // let mut start = Instant::now();

    // declare variable to store the common cherries
    let mut types_of_cherries: Vec<Vec<Cherry>> = Vec::new();
    let mut n1_common_cherries_tree: BTreeSet<CherryInfoMultiLabel> = BTreeSet::new();
    let mut n2_common_cherries_tree: BTreeSet<CherryInfoMultiLabel> = BTreeSet::new();
    let mut n1_unique_cherries_tree: BTreeSet<CherryInfoMultiLabel> = BTreeSet::new();

    let n1_common_cherries: Vec<Cherry>;
    let n2_common_cherries: Vec<Cherry>;
    let n1_unique_cherries: Vec<Cherry>;
    let n2_unique_cherries: Vec<Cherry>;

    // step 0, make new btreesets
    let mut s1: BTreeSet<CherryInfoMultiLabel> = BTreeSet::new();
    let mut s2: BTreeSet<CherryInfoMultiLabel> = BTreeSet::new();

    // step 1, get all labels associated with each cherry and put them into a new object
    // step 2, put that object into btreeset

    for current_n1_cherry in n1_cherries.into_iter() {
        s1.insert(CherryInfoMultiLabel {
            label1: n1.get_labels(current_n1_cherry.cherry.0),
            label2: n1.get_labels(current_n1_cherry.cherry.1),
            cherry_object: current_n1_cherry.clone(),
        });
    }

    for current_n2_cherry in n2_cherries.into_iter() {
        s2.insert(CherryInfoMultiLabel {
            label1: n2.get_labels(current_n2_cherry.cherry.0),
            label2: n2.get_labels(current_n2_cherry.cherry.1),
            cherry_object: current_n2_cherry.clone(),
        });
    }

    // step 3, use various set sorting methods to organize it appropriately

    for x in s1.into_iter() {
        if s2.contains(&x) {
            n2_common_cherries_tree.insert(s2.take(&x).unwrap());
            n1_common_cherries_tree.insert(x);
        } else {
            n1_unique_cherries_tree.insert(x);
        }
    }

    n1_common_cherries = n1_common_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n2_common_cherries = n2_common_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n1_unique_cherries = n1_unique_cherries_tree
        .into_iter()
        .map(|p| p.cherry_object)
        .collect();
    n2_unique_cherries = s2.into_iter().map(|p| p.cherry_object).collect();

    types_of_cherries.push(n1_common_cherries);
    types_of_cherries.push(n2_common_cherries);
    types_of_cherries.push(n1_unique_cherries);
    types_of_cherries.push(n2_unique_cherries);

    // println!("{:?} time", start.elapsed());

    types_of_cherries
}
/**
 * get_reducible_cherries_multi_label
 *
 * Given a network it retrieves a list of reducible cherries from the network then removes half the simple cherries (so if {a,b} and {b,a} are both simple cherries it removes one)
 * @param &Network<Node<String>> input_network - input network to get cherries from
 * @returns Vec<Vec<Cherry>> - list of trimmed cherries
 */
fn get_reducible_cherries_multi_label(input_network: &Network<Node<String>>) -> Vec<Cherry> {
    // get the non trimmed list of cherries
    let nontrimmed_cherries: Vec<Cherry> = get_reducible_cherries(&input_network);
    let mut trimmed_cherries: Vec<Cherry> = Vec::new();

    // iterate through the list and remove all duplicates
    for current_cherry in nontrimmed_cherries.iter() {
        // by default we add the cherry unless we find the opposite cherry in the list of added cherries already
        let mut add_cherry: bool = true;

        // reticulated cherries must both be added since the label ends up associated with different nodes in each cherry
        if current_cherry.cherry_type == Simple {
            // see if we have added the opposite to the vector already
            for current_trimmed_cherry in trimmed_cherries.iter() {
                if current_cherry.cherry.0 == current_trimmed_cherry.cherry.1
                    && current_cherry.cherry.1 == current_trimmed_cherry.cherry.0
                    && current_cherry.cherry_type == current_trimmed_cherry.cherry_type
                {
                    add_cherry = false;
                    break;
                }
            }
        }

        // if not then we add return
        if add_cherry {
            trimmed_cherries.push(current_cherry.clone());
        }
    }

    trimmed_cherries
}
/**
 * update_reducible_cherries_multi_label
 *
 * Given a network it retrieves a list of updated reducible cherries from the network then removes half the simple cherries (so if {a,b} and {b,a} are both simple cherries it removes one)
 * @param &Network<Node<String>> input_network - input network to get cherries from
 * @param &Cherry reduced_cherry - cherry that was just reduced, needed for network function
 * @param bool is_reduction - needed for network function
 * @param Vec<Cherry> cur_reducible_cherries - current list of cherries, needed for network function
 * @returns Vec<Vec<Cherry>> - list of trimmed cherries
 */
fn update_reducible_cherries_multi_label(
    n1: &Network<Node<String>>,
    reduced_cherry: &Cherry,
    is_reduction: bool,
    cur_reducible_cherries: Vec<Cherry>,
) -> Vec<Cherry> {
    // get the non trimmed list of cherries
    let nontrimmed_cherries: Vec<Cherry> =
        update_reducible_cherries(&n1, &reduced_cherry, is_reduction, cur_reducible_cherries);
    let mut trimmed_cherries: Vec<Cherry> = Vec::new();

    // iterate through the list and remove all duplicates
    // because the list returns the cherries ordered we only need to check the one most recently added
    for current_cherry in nontrimmed_cherries.iter() {
        let mut add_cherry: bool = true;

        // see if we have added the opposite to the vector already
        for current_trimmed_cherry in trimmed_cherries.iter() {
            if current_cherry.cherry.0 == current_trimmed_cherry.cherry.1
                && current_cherry.cherry.1 == current_trimmed_cherry.cherry.0
                && current_cherry.cherry_type == current_trimmed_cherry.cherry_type
            {
                add_cherry = false;
                break;
            }
        }

        // if not then we return it
        if add_cherry {
            trimmed_cherries.push(current_cherry.clone());
        }
    }

    trimmed_cherries
}
/**
 * predict_cherries
 *
 * Finds and returns the number of common cherries if you pick the selected cherries in the given networks (single label version)
 * @param &Network<Node<String>> n1 - input network
 * @param &Network<Node<String>> n2 - input network
 * @param Option<Cherry> selected_n1_cherry - cherry to pick from n1, may be None
 * @param Option<Cherry> selected_n2_cherry - cherry to pick from n2, may be None
 * @param Vec<Cherry> n1_reducible_cherries - list of available reducible cherries in n1
 * @param Vec<Cherry> n2_reducible_cherries - list of available reducible cherries in n2
 * @param usize iterations_left - countdown to base case, starts at depth (user input)
 * @param (usize, usize) current_num_cherries - current best number of common/unique cherries found so far
 * @param usize max_unique_cherries - Used for branch and bound, if number of unique cherries exceeds this then end early.
 * @param HashMap<String, (usize, usize)> in_cycle_results_hash - hashmap mapping a certain combination of networks to the max number of common/unique cherries found in it (avoid repeat work)
 * @param HashMap<String, (usize, usize)> cherry_types_hash - hashmap mapping a certain combination of networks to the lists of common / unique cherries available (avoid repeat work)
 * @returns (usize, usize) - number of common and unique cherries found in best case, ie max number of common cherries found + number of unique cherries in that version
 */
fn predict_cherries(
    mut n1: Network<Node<String>>,
    mut n2: Network<Node<String>>,
    selected_n1_cherry: Option<Cherry>,
    selected_n2_cherry: Option<Cherry>,
    mut n1_reducible_cherries: Vec<Cherry>,
    mut n2_reducible_cherries: Vec<Cherry>,
    iterations_left: usize,
    mut current_num_cherries: (usize, usize),
    max_unique_cherries: usize,
    mut in_cycle_results_hash: &mut HashMap<String, (usize, usize)>,
    mut cherry_types_hash: &mut HashMap<String, Vec<Vec<Cherry>>>,
) -> (usize, usize) {
    let num_cherry_types: (usize, usize); // stores max number of common cherries available from this current network at given depth (and number of unique cherries associated with it)
    let string_representation; // stores string representation of networks for hashing

    // case where n1 is already reduced
    if n1.get_n() == 1 {
        // n1 has 0 cherries, so 0 new common cherries and all n2 cherries are unique
        num_cherry_types = (
            current_num_cherries.0,
            current_num_cherries.1 + n2_reducible_cherries.len(),
        );
    }
    // case where n2 is already reduced
    else if n2.get_n() == 1 {
        // n2 has 0 cherries, so 0 new common cherries and all n1 cherries are unique
        num_cherry_types = (
            current_num_cherries.0,
            current_num_cherries.1 + n1_reducible_cherries.len(),
        );
    }
    // case where n1 is already almost reduced so we can no longer reduce it further
    else if n1.get_n() == 3 {
        // check whether the remaining cherry is equal to anything in n2
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types =
                get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // if len > 0 that means that the single cherry n1 has is in common with n2 and the rest are unique
        if cherry_types[0].len() > 0 {
            num_cherry_types = (
                current_num_cherries.0 + 1,
                current_num_cherries.1 + cherry_types[3].len(),
            );
        }
        // otherwise it is not common and everything is unique
        else {
            num_cherry_types = (
                current_num_cherries.0,
                current_num_cherries.1 + cherry_types[3].len(),
            );
        }
    }
    // case where n2 is already almost reduced so we can no longer reduce it further
    else if n2.get_n() == 3 {
        // check whether the remaining cherry is equal to anything in n2
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types =
                get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // if len > 0 that means that the single cherry n2 has is in common with n1 and rest is unique
        if cherry_types[0].len() > 0 {
            num_cherry_types = (
                current_num_cherries.0 + 1,
                current_num_cherries.1 + cherry_types[2].len(),
            );
        }
        // otherwise it is not common and everything is unique
        else {
            num_cherry_types = (
                current_num_cherries.0 + 0,
                current_num_cherries.1 + cherry_types[2].len(),
            );
        }
    }
    // case where there are no more recursion left to perform
    else if iterations_left == 0 {
        // now, we want to simulate removing the cherry(s) then counting how many common and unique cherries there are
        // so, we check if n1 cherry is null and if not we remove that
        // next, we check if n2 cherry is null and if not we remove that
        // then we count all the common and unique cherries and return

        // reduce the network by the selected cherries and recalculate reducible cherries if they exist
        match selected_n1_cherry {
            Some(p) => {
                reduce_cherry(&mut n1, &p);
                n1_reducible_cherries =
                    update_reducible_cherries(&n1, &p, true, n1_reducible_cherries);
            }
            None => {}
        }

        match selected_n2_cherry {
            Some(p) => {
                reduce_cherry(&mut n2, &p);
                n2_reducible_cherries =
                    update_reducible_cherries(&n2, &p, true, n2_reducible_cherries);
            }
            None => {}
        }

        // count the number of common cherries and the number of unique cherries
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types =
                get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // add all current common and unique cherries
        num_cherry_types = (
            current_num_cherries.0 + cherry_types[0].len() + cherry_types[1].len(),
            current_num_cherries.1 + cherry_types[2].len() + cherry_types[3].len(),
        );
    }
    // otherwise, reduce the network by the selected cherry and simulate removing every other cherry and keeping track of the maximum possible number of common cherries from each
    else {
        // reduce the network by the selected cherries then recalculate reducible cherries
        match selected_n1_cherry {
            Some(p) => {
                reduce_cherry(&mut n1, &p);
                n1_reducible_cherries =
                    update_reducible_cherries(&n1, &p, true, n1_reducible_cherries);
            }
            None => {}
        }

        match selected_n2_cherry {
            Some(p) => {
                reduce_cherry(&mut n2, &p);
                n2_reducible_cherries =
                    update_reducible_cherries(&n2, &p, true, n2_reducible_cherries);
            }
            None => {}
        }

        // create variables to store the types of cherries
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types =
                get_cherry_types(&n1, &n2, &n1_reducible_cherries, &n2_reducible_cherries);
            cherry_types_hash.insert(string_representation.clone(), cherry_types.clone());
        }

        // update cherries found so far
        current_num_cherries.0 =
            current_num_cherries.0 + cherry_types[0].len() + cherry_types[1].len();
        current_num_cherries.1 =
            current_num_cherries.1 + cherry_types[2].len() + cherry_types[3].len();

        // check to see if we have seen this specific network combo before
        if in_cycle_results_hash.contains_key(&string_representation) {
            // if so simply skip straight to the base case
            num_cherry_types = in_cycle_results_hash[&string_representation];
        }
        // check to see if current number of unique cherries is more than the number found in our best answer so far
        else if current_num_cherries.1 > max_unique_cherries {
            // if so, just cut the branch here
            num_cherry_types = current_num_cherries;
        }
        // if it does not then proceed as normal
        else {
            // create variables that are used to store the best answer found so far
            let mut temp_num_cherry_types;
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early

            // in this case, we remove 1 cherry from n1
            if n1.get_n() > n2.get_n() {
                // if unique cherries available only use those
                if cherry_types[2].len() > 0 {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n1_cherry) in cherry_types[2].iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries(
                            n1.clone(),
                            n2.clone(),
                            Some(current_n1_cherry.clone()),
                            None::<Cherry>,
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
                // otherwise use common cherries
                else {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n1_cherry) in n1_reducible_cherries.iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries(
                            n1.clone(),
                            n2.clone(),
                            Some(current_n1_cherry.clone()),
                            None::<Cherry>,
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
            }
            // in this case, we remove 1 cherry from n2
            else {
                // if unique cherries available use them
                if cherry_types[3].len() > 0 {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n2_cherry) in cherry_types[3].iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries(
                            n1.clone(),
                            n2.clone(),
                            None::<Cherry>,
                            Some(current_n2_cherry.clone()),
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
                // otherwise use common cherries
                else {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n2_cherry) in n2_reducible_cherries.iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries(
                            n1.clone(),
                            n2.clone(),
                            None::<Cherry>,
                            Some(current_n2_cherry.clone()),
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
            }

            num_cherry_types = most_common_cherries;

            // update the hashmap of networks
            in_cycle_results_hash.insert(string_representation, num_cherry_types);
        }
    }

    num_cherry_types
}

/**
 * predict_cherries_multi_label
 *
 * Finds and returns the number of common cherries if you pick the selected cherries in the given networks (multi label version)
 * @param &Network<Node<String>> n1 - input network
 * @param &Network<Node<String>> n2 - input network
 * @param Option<Cherry> selected_n1_cherry - cherry to pick from n1, may be None
 * @param Option<Cherry> selected_n2_cherry - cherry to pick from n2, may be None
 * @param Vec<Cherry> n1_reducible_cherries - list of available reducible cherries in n1
 * @param Vec<Cherry> n2_reducible_cherries - list of available reducible cherries in n2
 * @param usize iterations_left - countdown to base case, starts at depth (user input)
 * @param (usize, usize) current_num_cherries - current best number of common/unique cherries found so far
 * @param usize max_unique_cherries - Used for branch and bound, if number of unique cherries exceeds this then end early.
 * @param HashMap<String, (usize, usize)> in_cycle_results_hash - hashmap mapping a certain combination of networks to the max number of common/unique cherries found in it (avoid repeat work)
 * @param HashMap<String, (usize, usize)> cherry_types_hash - hashmap mapping a certain combination of networks to the lists of common / unique cherries available (avoid repeat work)
 * @returns (usize, usize) - number of common and unique cherries found in best case, ie max number of common cherries found + number of unique cherries in that version
 */
fn predict_cherries_multi_label(
    mut n1: Network<Node<String>>,
    mut n2: Network<Node<String>>,
    selected_n1_cherry: Option<Cherry>,
    selected_n2_cherry: Option<Cherry>,
    mut n1_reducible_cherries: Vec<Cherry>,
    mut n2_reducible_cherries: Vec<Cherry>,
    iterations_left: usize,
    mut current_num_cherries: (usize, usize),
    max_unique_cherries: usize,
    mut in_cycle_results_hash: &mut HashMap<String, (usize, usize)>,
    mut cherry_types_hash: &mut HashMap<String, Vec<Vec<Cherry>>>,
) -> (usize, usize) {
    let num_cherry_types: (usize, usize); // stores max number of common cherries available from this current network at given depth (and number of unique cherries associated with it)
    let string_representation; // stores string representation of networks for hashing

    // case where n1 is already reduced
    if n1.get_n() == 1 {
        // n1 has 0 cherries, so 0 new common cherries and all n2 cherries are unique
        num_cherry_types = (
            current_num_cherries.0,
            current_num_cherries.1 + n2_reducible_cherries.len(),
        );
    }
    // case where n2 is already reduced
    else if n2.get_n() == 1 {
        // n2 has 0 cherries, so 0 new common cherries and all n1 cherries are unique
        num_cherry_types = (
            current_num_cherries.0,
            current_num_cherries.1 + n1_reducible_cherries.len(),
        );
    }
    // case where n1 is already almost reduced so we can no longer reduce it further
    else if n1.get_n() == 3 {
        // check whether the remaining cherry is equal to anything in n2
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types = get_multi_label_cherry_types(
                &n1,
                &n2,
                &n1_reducible_cherries,
                &n2_reducible_cherries,
            );
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // if len > 0 that means that the single cherry n1 has is in common with n2 and the rest are unique
        if cherry_types[0].len() > 0 {
            num_cherry_types = (
                current_num_cherries.0 + 1,
                current_num_cherries.1 + cherry_types[3].len(),
            );
        }
        // otherwise it is not common and everything is unique
        else {
            num_cherry_types = (
                current_num_cherries.0,
                current_num_cherries.1 + cherry_types[3].len(),
            );
        }
    }
    // case where n2 is already almost reduced so we can no longer reduce it further
    else if n2.get_n() == 3 {
        // check whether the remaining cherry is equal to anything in n2
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types = get_multi_label_cherry_types(
                &n1,
                &n2,
                &n1_reducible_cherries,
                &n2_reducible_cherries,
            );
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // if len > 0 that means that the single cherry n2 has is in common with n1 and rest is unique
        if cherry_types[0].len() > 0 {
            num_cherry_types = (
                current_num_cherries.0 + 1,
                current_num_cherries.1 + cherry_types[2].len(),
            );
        }
        // otherwise it is not common and everything is unique
        else {
            num_cherry_types = (
                current_num_cherries.0 + 0,
                current_num_cherries.1 + cherry_types[2].len(),
            );
        }
    }
    // case where there are no more recursion left to perform
    else if iterations_left == 0 {
        // now, we want to simulate removing the cherry(s) then counting how many common and unique cherries there are
        // so, we check if n1 cherry is null and if not we remove that
        // next, we check if n2 cherry is null and if not we remove that
        // then we count all the common and unique cherries and return

        // reduce the network by the selected cherries and recalculate reducible cherries if they exist
        match selected_n1_cherry {
            Some(p) => {
                reduce_cherry_multi_label(&mut n1, &p);
                n1_reducible_cherries =
                    update_reducible_cherries_multi_label(&n1, &p, true, n1_reducible_cherries);
            }
            None => {}
        }

        match selected_n2_cherry {
            Some(p) => {
                reduce_cherry_multi_label(&mut n2, &p);
                n2_reducible_cherries =
                    update_reducible_cherries_multi_label(&n2, &p, true, n2_reducible_cherries);
            }
            None => {}
        }

        // count the number of common cherries and the number of unique cherries
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types = get_multi_label_cherry_types(
                &n1,
                &n2,
                &n1_reducible_cherries,
                &n2_reducible_cherries,
            );
            cherry_types_hash.insert(string_representation, cherry_types.clone());
        }

        // add all current common and unique cherries
        num_cherry_types = (
            current_num_cherries.0 + cherry_types[0].len() + cherry_types[1].len(),
            current_num_cherries.1 + cherry_types[2].len() + cherry_types[3].len(),
        );
    }
    // otherwise, reduce the network by the selected cherry and simulate removing every other cherry and keeping track of the maximum possible number of common cherries from each
    else {
        // reduce the network by the selected cherries then recalculate reducible cherries
        match selected_n1_cherry {
            Some(p) => {
                reduce_cherry_multi_label(&mut n1, &p);
                n1_reducible_cherries =
                    update_reducible_cherries_multi_label(&n1, &p, true, n1_reducible_cherries);
            }
            None => {}
        }

        match selected_n2_cherry {
            Some(p) => {
                reduce_cherry_multi_label(&mut n2, &p);
                n2_reducible_cherries =
                    update_reducible_cherries_multi_label(&n2, &p, true, n2_reducible_cherries);
            }
            None => {}
        }

        // create variables to store the types of cherries
        let cherry_types;

        // get hash key
        string_representation = get_hash_key(&n1, &n2);

        // calculate the new common cherries and unique cherries (if not already calculated)
        if cherry_types_hash.contains_key(&string_representation) {
            cherry_types = cherry_types_hash[&string_representation].clone();
        } else {
            cherry_types = get_multi_label_cherry_types(
                &n1,
                &n2,
                &n1_reducible_cherries,
                &n2_reducible_cherries,
            );
            cherry_types_hash.insert(string_representation.clone(), cherry_types.clone());
        }

        // update cherries found so far
        current_num_cherries.0 =
            current_num_cherries.0 + cherry_types[0].len() + cherry_types[1].len();
        current_num_cherries.1 =
            current_num_cherries.1 + cherry_types[2].len() + cherry_types[3].len();

        // check to see if we have seen this specific network combo before
        if in_cycle_results_hash.contains_key(&string_representation) {
            // if so simply skip straight to the base case
            num_cherry_types = in_cycle_results_hash[&string_representation];
        }
        // check to see if current number of unique cherries is more than the number found in our best answer so far
        else if current_num_cherries.1 > max_unique_cherries {
            // if so, just cut the branch here
            num_cherry_types = current_num_cherries;
        }
        // if it does not then proceed as normal
        else {
            // create variables that are used to store the best answer found so far
            let mut temp_num_cherry_types;
            let mut most_common_cherries: (usize, usize) = (0, usize::MAX); // starts with 0 for common cherries so anything is bigger and MAX so we never cut the first branch off early

            // in this case, we remove 1 cherry from n1
            if n1.get_n() > n2.get_n() {
                // if unique cherries available only use those
                if cherry_types[2].len() > 0 {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n1_cherry) in cherry_types[2].iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries_multi_label(
                            n1.clone(),
                            n2.clone(),
                            Some(current_n1_cherry.clone()),
                            None::<Cherry>,
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
                // otherwise use common cherries
                else {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n1_cherry) in n1_reducible_cherries.iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries_multi_label(
                            n1.clone(),
                            n2.clone(),
                            Some(current_n1_cherry.clone()),
                            None::<Cherry>,
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
            }
            // in this case, we remove 1 cherry from n2
            else {
                // if unique cherries available use them
                if cherry_types[3].len() > 0 {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n2_cherry) in cherry_types[3].iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries_multi_label(
                            n1.clone(),
                            n2.clone(),
                            None::<Cherry>,
                            Some(current_n2_cherry.clone()),
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
                // otherwise use common cherries
                else {
                    //  try every possible cherry and then store the results in a vector
                    for (i, current_n2_cherry) in n2_reducible_cherries.iter().enumerate() {
                        // call itself recursively with same networks, selected cherry, with depth - 1
                        temp_num_cherry_types = predict_cherries_multi_label(
                            n1.clone(),
                            n2.clone(),
                            None::<Cherry>,
                            Some(current_n2_cherry.clone()),
                            n1_reducible_cherries.clone(),
                            n2_reducible_cherries.clone(),
                            iterations_left - 1,
                            current_num_cherries,
                            most_common_cherries.1,
                            &mut in_cycle_results_hash,
                            &mut cherry_types_hash,
                        );

                        // first check to see this is the first iteration
                        if i == 0 {
                            // if so then this is the best even if it has 0 common cherries
                            most_common_cherries = temp_num_cherry_types;
                        }
                        // otherwise check to see if this version has more common cherries
                        else if temp_num_cherry_types.0 > most_common_cherries.0 {
                            most_common_cherries = temp_num_cherry_types;
                        }
                    }
                }
            }

            num_cherry_types = most_common_cherries;

            // update the hashmap of networks
            in_cycle_results_hash.insert(string_representation, num_cherry_types);
        }
    }

    num_cherry_types
}

/**
 * get_hash_key
 *
 * Given two networks reduced unique hashkey to be used in hashmaps
 * @param &Network<Node<String>> n1 - input network
 * @param &Network<Node<String>> n2 - input network
 * @returns String - unqiue hashkey based on newick format of each network
 */
fn get_hash_key(n1: &Network<Node<String>>, n2: &Network<Node<String>>) -> String {
    return format!("{};{}", n1.to_string(), n2.to_string());
}

/**
 * print_cherrypicking_sequence
 *
 * Given a network and a list of available reducible cherries from that network, print the cherries (currently not used)
 * @param &Network<Node<String>> n1 - input network
 * @param &Vec<Cherry> sequence - list of reducible cherries
 */
fn print_cherrypicking_sequence(n1: &Network<Node<String>>, sequence: &Vec<Cherry>) {
    let mut i = 0;

    print!("{{");

    while i < sequence.len() {
        print!(
            "({0}, {1}, ",
            n1.get_label(sequence[i].cherry.0),
            n1.get_label(sequence[i].cherry.1)
        );

        if sequence[i].cherry_type == Retic {
            print!("R");
        } else {
            print!("S");
        }

        if i != (sequence.len() - 1) {
            print!("), ");
        } else {
            print!(")");
        }

        i = i + 1;
    }

    print!("}}\n");
}

/**
 * print_CherryInfo_sequence
 *
 * Given a BTreeSet of CherryInfo objects it goes through and prints all of them (currently not used)
 * @param &BTreeSet<CherryInfo> sequence - list of CherryInfo objects
 */
fn print_CherryInfo_sequence(sequence: &BTreeSet<CherryInfo>) {
    print!("{{");

    for (i, current) in sequence.iter().enumerate() {
        print!("({0}, {1}, ", current.label1, current.label2);

        if current.cherry_object.cherry_type == Retic {
            print!("R");
        } else {
            print!("S");
        }

        if i != (sequence.len() - 1) {
            print!("), ");
        } else {
            print!(")");
        }
    }

    print!("}}\n");
}
