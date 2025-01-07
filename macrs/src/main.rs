use macrs::network::cherry::CherryType;
use macrs::network::cherry::CherryType::Retic;
use macrs::network::distance::DTree;
use macrs::network::Network;
use macrs::network::Node;
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::time::Instant;

fn main() {
    env::set_var("RUST_BACKTRACE", "1");
    let args: Vec<String> = env::args().collect();
    //     ---------        set up environment variables
    let mut is_debug = false;
    let mut is_multi_label = false;
    let mut is_exact = false;
    let mut is_gen_mode = false;
    let mut is_heuristic_mode = false;
    let mut req_args: Vec<&str> = vec![];
    //       --------------             parse command
    match args[1].as_str() {
        "gen" => is_gen_mode = true,
        "new" => is_gen_mode = false,
        "heuristic" => is_heuristic_mode = true,
        _ => {
            help("error: provide command for mode.");
            return;
        }
    }
    //    -------------          parse option flags
    //short form and long form, in alphabetical order
    for (i, arg) in args.iter().enumerate() {
        if i == 0 || i == 1 {
            continue;
        }
        if arg.starts_with("--") {
            let option = arg.as_str().trim_start_matches('-');
            println!("option: {option}");
        } else if arg.starts_with('-') {
            let options = arg.trim_start_matches('-');
            for c in options.chars() {
                match c {
                    'd' => is_debug = true,
                    'v' => is_debug = true, //legacy "verbose"
                    'e' => is_exact = true,
                    'l' => is_multi_label = true,
                    'h' => {
                        help("printing help...");
                        return;
                    }
                    'm' => {
                        println!("a link to the manual.");
                    } //TODO
                    _ => {
                        print!("{c} ");
                        help("unknown flag encountered");
                        return;
                    }
                }
            }
        } else {
            req_args.push(arg);
        }
    }
    if is_heuristic_mode {
        if req_args.len() != 3 {
            help("Not the required number of arguments for heuristic");
            return;
        }

        let file_result_1 = fs::read_to_string(req_args[0]);
        let newick1 = match file_result_1 {
            Ok(new_string) => new_string,
            Err(_) => {
                help("Could not read file 1.");
                return;
            }
        };

        let file_result_2 = fs::read_to_string(req_args[1]);
        let newick2 = match file_result_2 {
            Ok(new_string) => new_string,
            Err(_) => {
                help("Could not read file 2.");
                return;
            }
        };

        let num_iterations;
        let check_numeric = req_args[2].trim().parse::<usize>();

        match check_numeric {
            Ok(ok) => {
                num_iterations = ok;
            }
            Err(_) => {
                help("heuristic mode required argument is not numeric.");
                return;
            }
        }

        call_heuristic_mode(is_multi_label, newick1, newick2, num_iterations);
    } else if is_gen_mode {
        let mut leaves = 0;
        let mut reticulations = 0;
        let mut distance = 0;
        for (i, req_arg_result) in req_args.into_iter().enumerate() {
            if let Ok(req_arg) = req_arg_result.parse::<usize>() {
                match i {
                    0 => leaves = req_arg,
                    1 => reticulations = req_arg,
                    2 => distance = req_arg,
                    _ => {
                        help("Not the required number of arguments for random generation mode");
                        return;
                    }
                }
            } else {
                help("gen mode required argument is not numeric.");
                return;
            }
        }
        //      ------------   gen   call distance with given args and options
        call_gen_mode(is_debug, is_exact, leaves, reticulations, distance);
    } else {
        /*
        if req_args.len() != 2 {
            help("Not the required number of arguments for extended newick format string parser mode");
            return;
        }
        let file_result_1 = fs::read_to_string(req_args[0]);
        let newick1 = match file_result_1 {
            Ok(new_string) => new_string,
            Err(_) => {
                help("Could not read file 1.");
                return;
            }
        };
        let file_result_2 = fs::read_to_string(req_args[1]);
        let newick2 = match file_result_2 {
            Ok(new_string) => new_string,
            Err(_) => {
                help("Could not read file 2.");
                return;
            }
        };
        // ---------   new     call distance with selected options and arguments
        call_new_mode(is_debug, newick1, newick2);
        */
    }
}
fn help(message: &str) {
    println!("exit message: {}", message);
    println!("
    Cherry distance program calculates and outputs cherry distance to stdout. It runs in two modes: 

    COMMANDS
    gen     random generation mode
    new     extended format newick string parsing mode
    heuristic     Estimation mode
    
    Random generation mode will generate a network with the specified number of leaves and up to the specified number of reticulations (see --exact option below), then randomly modifies the generated network by the given distance. The calculated distance may be exactly 2 less than given (see manual elaboration). 

    Extended format newick string parsing mode takes two required arguments, two files which should each contain an extended Newick format string. The program will parse each into a network, outputting the resulting parsed network each into an output file. 

    Estimation mode takes 3 arguments, two files which should each contain an extended Newick format string and an integer specifying what level of depth to look at. The program will parse each into a network, then estimate the distance between them. 

    USAGE
    1. random reneration mode
    usage: ./dist gen [options] <leaves> <reticulations> <distance> 

    2. extended newick format string parsing mode
    usage: ./dist new [options] <file1> <file2>

    3. Estimation mode
    usage: ./dist heuristic [options] <file1> <file2> <depth(integer)>
    
    OPTIONS
    -d, --debug        verbose printing mode outputs mainly runtime information
    -e, --exact        option only available on random generation mode, specifies that the exact number of reticulations requested should be reached in randomly generated network, note this may increase runtime
    -h, --help         print this help guide   
    -l, --multi label  If running in heuristic mode, then use multi labels instead of single labels 
    -m, --manual    links to version of manual on web

    ");
}
fn call_gen_mode(debug: bool, exact: bool, leaves: usize, reticulations: usize, distance: usize) {
    //let dist;
    let n1 = Network::new_random(leaves, reticulations, exact);
    let (n1, n2) = Network::random_modify(n1, distance);
    println!("{n1}; {n2};");
}
fn call_new_mode(debug: bool, new1: String, new2: String) {
    let n1 = Network::parse_newick(&new1);
    let n2 = Network::parse_newick(&new2);
    //calculate distance
    //let dist = Network::find_cherry_distance_smart_enum(n1, n2, debug);
    //last output
    //println!("Cherry distance of random networks: {dist}");
}
/**
 * call_heuristic_mode
 *
 * Call functions, with appropriate parameters, that find the estimate of the CPS, then calculate the tail distance
 * @param bool is_multi_label - input network
 * @param String new1 - network 1 in newick
 * @param String new2 - network 2 in newick
 * @param usize num_iterations - depth of the look ahead recursive calls
 */
fn call_heuristic_mode(is_multi_label: bool, new1: String, new2: String, num_iterations: usize) {
    let n1: Network<Node<String>> = Network::parse_newick(&new1);
    let n2: Network<Node<String>> = Network::parse_newick(&new2);
    //calculate sequence
    if is_multi_label {
        // get multi label CPS
        let sequences =
            Network::<Node<String>>::find_heuristic_sequence_multi_label(n1, n2, num_iterations);

        // use function to calculate tail size and then use tail distance formula
        print!(
            "{0}",
            sequences[0].len() + sequences[1].len()
                - 2 * find_single_label_cps_distance_from_mutli_label(
                    sequences[0].clone(),
                    sequences[1].clone()
                )
                .pop_last()
                .unwrap()
        );
    } else {
        let sequences = Network::<Node<String>>::find_heuristic_sequence(n1, n2, num_iterations);
        // print_cherrypicking_sequence(&sequences[0]);
        // print_cherrypicking_sequence(&sequences[1]);
        print!(
            "{0}",
            find_cherrypicking_distance(&sequences[0], &sequences[1])
        );
    }
}
/**
 * print_cherrypicking_sequence
 *
 * Given a list of cherries (labels and type), prints the list (single label version)
 * @param &Vec<(String, String, CherryType)> sequence - sequence to print
 */
fn print_cherrypicking_sequence(sequence: &Vec<(String, String, CherryType)>) {
    let mut i = 0;

    print!("{{");

    while i < sequence.len() {
        print!("({0}, {1}, ", sequence[i].0, sequence[i].1);

        if sequence[i].2 == Retic {
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
 * print_cherrypicking_sequence_multi_label
 *
 * Given a list of cherries (labels and type), prints the list (multi label version)
 * @param &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> sequence - sequence to print
 */
fn print_cherrypicking_sequence_multi_label(
    sequence: &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
) {
    let mut i = 0;

    print!("{{");

    while i < sequence.len() {
        print!("(");

        print!("[");

        for current in sequence[i].0.iter() {
            print!("{0} ", current);
        }

        print!("], ");

        print!("[");

        for current in sequence[i].1.iter() {
            print!("{0} ", current);
        }

        print!("], ");

        if sequence[i].2 == Retic {
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
 * find_cherrypicking_distance
 *
 * Given CPSs for two networks returns the tail distance (single label version)
 * @param &Vec<(String, String, CherryType)> sequence_1 - CPS 1
 * @param &Vec<(String, String, CherryType)> sequence_2 - CPS 2
 * @returns usize - tail distance between the two networks the CPSs reduce
 */
fn find_cherrypicking_distance(
    sequence_1: &Vec<(String, String, CherryType)>,
    sequence_2: &Vec<(String, String, CherryType)>,
) -> usize {
    let mut tail_size: usize = 0;
    let result: usize;
    let mut sequence_1_iter = (sequence_1.iter()).rev();
    let mut sequence_2_iter = (sequence_2.iter()).rev();

    // loop until there is no more macrs
    loop {
        let seq_1_cherry;
        let seq_2_cherry;

        // get the next value from seq 1
        match sequence_1_iter.next() {
            Some(x) => {
                seq_1_cherry = x;
            }
            None => break,
        }

        // get the next value from seq 2
        match sequence_2_iter.next() {
            Some(x) => {
                seq_2_cherry = x;
            }
            None => break,
        }

        // check to see if the values match
        if seq_1_cherry.0 == seq_2_cherry.0
            && seq_1_cherry.1 == seq_2_cherry.1
            && seq_1_cherry.2 == seq_2_cherry.2
        {
            tail_size = tail_size + 1;
        } else {
            break;
        }
    }

    // ensure this is the same distance that Kaari uses
    result = sequence_1.len() + sequence_2.len() - 2 * tail_size;
    result
}

/**
 * find_cherrypicking_distance_mutli_label
 *
 * Given CPSs for two networks returns the tail distance (multi label version)
 * @param &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> sequence_1 - CPS 1
 * @param &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> - CPS 2
 * @returns usize - tail distance between the two networks the CPSs reduce
 */
fn find_cherrypicking_distance_mutli_label(
    sequence_1: &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    sequence_2: &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
) -> usize {
    let mut tail_size: usize = 0;
    let result: usize;
    let mut sequence_1_iter = (sequence_1.iter()).rev();
    let mut sequence_2_iter = (sequence_2.iter()).rev();

    // loop until there is no more macrs
    loop {
        let seq_1_cherry;
        let seq_2_cherry;

        // get the next value from seq 1
        match sequence_1_iter.next() {
            Some(x) => {
                seq_1_cherry = x;
            }
            None => break,
        }

        // get the next value from seq 2
        match sequence_2_iter.next() {
            Some(x) => {
                seq_2_cherry = x;
            }
            None => break,
        }

        // check to see if the values match
        if seq_1_cherry.0.intersection(&seq_2_cherry.0).count() > 0
            && seq_1_cherry.1.intersection(&seq_2_cherry.1).count() > 0
            && seq_1_cherry.2 == seq_2_cherry.2
            || seq_1_cherry.0.intersection(&seq_2_cherry.1).count() > 0
                && seq_1_cherry.1.intersection(&seq_2_cherry.0).count() > 0
                && seq_1_cherry.2 == seq_2_cherry.2
        {
            tail_size = tail_size + 1;
        } else {
            break;
        }
    }

    // ensure this is the same distance that Kaari uses
    result = sequence_1.len() + sequence_2.len() - 2 * tail_size;
    result
}

/**
 * find_cherrypicking_distance_mutli_label
 *
 * Given two multi label CPSs for two networks finds all single label CPSs (some of them may still have multi labels in first few cherries since the algorithm ends early if no matching labels can be found between a cherry pair)
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s1 - Multi Label CPS 1
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s2 - Multi Label CPS 2
 * @returns Vec<Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>>> - vector of all single label CPSs
 */
fn find_single_label_cps_from_mutli_label(
    s1: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    s2: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
) -> Vec<Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>>> {
    let mut result: Vec<Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>>> = Vec::new();
    find_single_label_cps_from_mutli_label_recursive(s1, s2, 1, &mut result);
    result
}

/**
 * find_single_label_cps_from_mutli_label_recursive
 *
 * Recursive helper for find_single_label_cps_from_mutli_label
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s1 - Multi Label CPS 1
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s2 - Multi Label CPS 2
 * @param usize end_index - Recursive counter indicating how far along in the CPSs we are. 1 = use last cherry pair, 2 = use second last cherry pair etc
 * @returns Vec<Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>>> - vector of all single label CPSs
 */
fn find_single_label_cps_from_mutli_label_recursive(
    s1: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    s2: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    end_index: usize,
    result: &mut Vec<Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>>>,
) {
    // check to see if we have hit recursive limit
    // the limit being we have already considered all cherries in one CPS
    if end_index > s1.len() || end_index > s2.len() {
        // if so, enter what we got into the results vector and end
        let mut sub_result: Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> = Vec::new();
        sub_result.push(s1);
        sub_result.push(s2);
        result.push(sub_result);
    } else {
        // loop through the labels for the cherry that was chosen
        let s1_cherry_index = s1.len() - end_index;
        let s2_cherry_index = s2.len() - end_index;

        // verify that the cherries are the same type
        if s1[s1_cherry_index].2 != s2[s2_cherry_index].2 {
            // if not, enter what we got into the results vector and end
            let mut sub_result: Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> =
                Vec::new();
            sub_result.push(s1);
            sub_result.push(s2);
            result.push(sub_result);
        }
        // otherwise we continue
        else {
            // variable determining whether or not we found at least one match, if not then we end here as the tail distance cannot be extended
            let mut found_one: bool = false;

            // compare Cherry A label list 1/2 to Cherry B label 1/2
            for label1_current in s1[s1_cherry_index].0.iter() {
                // check to see if s2 label1 has the same corresponding label
                if s2[s2_cherry_index].0.contains(&label1_current as &str) {
                    // if it does, then we have found 1/2 of the needed pairings
                    // next, we do the same thing for label2
                    for label2_current in s1[s1_cherry_index].1.iter() {
                        // check to see if s2 label2 has the same corresponding label
                        if s2[s2_cherry_index].1.contains(&label2_current as &str) {
                            // We found a match!
                            found_one = true;

                            // now we need to make a copy, remove all other labels, and then make a recursive call
                            let mut s1_copy = s1.clone();
                            let mut s2_copy = s2.clone();

                            // call helper
                            remove_all_other_labels(
                                &mut s1_copy,
                                &mut s2_copy,
                                label1_current.clone(),
                                label2_current.clone(),
                            );

                            // make recursive call to move onto the next cherry
                            find_single_label_cps_from_mutli_label_recursive(
                                s1_copy,
                                s2_copy,
                                end_index + 1,
                                result,
                            );
                        }
                    }
                }
            }

            // compare Cherry A label list 1/2 to Cherry B label 2/1
            for label1_current in s1[s1_cherry_index].0.iter() {
                // check to see if s2 label1 has the same corresponding label
                if s2[s2_cherry_index].1.contains(&label1_current as &str) {
                    // if it does, then we have found 1/2 of the needed pairings
                    // next, we do the same thing for label2
                    for label2_current in s1[s1_cherry_index].0.iter() {
                        // check to see if s2 label2 has the same corresponding label
                        if s2[s2_cherry_index].1.contains(&label2_current as &str) {
                            // We found a match!
                            found_one = true;

                            // now we need to make a copy, remove all other labels, and then make a recursive call
                            let mut s1_copy = s1.clone();
                            let mut s2_copy = s2.clone();

                            // call helper
                            remove_all_other_labels(
                                &mut s1_copy,
                                &mut s2_copy,
                                label1_current.clone(),
                                label2_current.clone(),
                            );

                            // make recursive call to move onto the next cherry
                            find_single_label_cps_from_mutli_label_recursive(
                                s1_copy,
                                s2_copy,
                                end_index + 1,
                                result,
                            );
                        }
                    }
                }
            }

            // if we did not find any matches then this is as long as the tail distance is going to get and we add it to the list
            if !found_one {
                let mut sub_result: Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> =
                    Vec::new();
                sub_result.push(s1);
                sub_result.push(s2);
                result.push(sub_result);
            }
        }
    }
}

/**
 * find_single_label_cps_distance_from_mutli_label
 *
 * Given two multi label CPSs for two networks finds all possible tail lengths from all possible single label CPSs
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s1 - Multi Label CPS 1
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s2 - Multi Label CPS 2
 * @returns BTreeSet<usize> - vector of all possible tail lengths
 */
fn find_single_label_cps_distance_from_mutli_label(
    s1: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    s2: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
) -> BTreeSet<usize> {
    let mut result: BTreeSet<usize> = BTreeSet::new();
    find_single_label_cps_distance_from_mutli_label_recursive(s1, s2, 1, &mut result);
    result
}

/**
 * find_single_label_cps_distance_from_mutli_label_recursive
 *
 * Recursive helper for find_single_label_cps_distance_from_mutli_label
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s1 - Multi Label CPS 1
 * @param Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s2 - Multi Label CPS 2
 * @param usize end_index - Recursive counter indicating how far along in the CPSs we are. 1 = use last cherry pair, 2 = use second last cherry pair etc
 * @returns BTreeSet<usize> - vector of all possible tail lengths
 */
fn find_single_label_cps_distance_from_mutli_label_recursive(
    s1: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    s2: Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    end_index: usize,
    result: &mut BTreeSet<usize>,
) {
    // check to see if we have hit recursive limit
    // the limit being we have already considered all cherries in one CPS
    if end_index > s1.len() || end_index > s2.len() {
        // if so, enter the tail size and end
        result.insert(end_index - 1);
    } else {
        // loop through the labels for the cherry that was chosen
        let s1_cherry_index = s1.len() - end_index;
        let s2_cherry_index = s2.len() - end_index;

        // verify that the cherries are the same type
        if s1[s1_cherry_index].2 != s2[s2_cherry_index].2 {
            // if so, enter the tail size and end
            result.insert(end_index - 1);
        }
        // otherwise we continue
        else {
            // variable determining whether or not we found at least one match, if not then we end here as the tail distance cannot be extended
            let mut found_one: bool = false;

            // compare Cherry A label list 1/2 to Cherry B label 1/2
            for label1_current in s1[s1_cherry_index].0.iter() {
                // check to see if s2 label1 has the same corresponding label
                if s2[s2_cherry_index].0.contains(&label1_current as &str) {
                    // if it does, then we have found 1/2 of the needed pairings
                    // next, we do the same thing for label2
                    for label2_current in s1[s1_cherry_index].1.iter() {
                        // check to see if s2 label2 has the same corresponding label
                        if s2[s2_cherry_index].1.contains(&label2_current as &str) {
                            // We found a match!
                            found_one = true;

                            // now we need to make a copy, remove all other labels, and then make a recursive call
                            let mut s1_copy = s1.clone();
                            let mut s2_copy = s2.clone();

                            // call helper
                            remove_all_other_labels(
                                &mut s1_copy,
                                &mut s2_copy,
                                label1_current.clone(),
                                label2_current.clone(),
                            );

                            // make recursive call to move onto the next cherry
                            find_single_label_cps_distance_from_mutli_label_recursive(
                                s1_copy,
                                s2_copy,
                                end_index + 1,
                                result,
                            );
                        }
                    }
                }
            }

            // compare Cherry A label list 1/2 to Cherry B label 2/1
            for label1_current in s1[s1_cherry_index].0.iter() {
                // check to see if s2 label1 has the same corresponding label
                if s2[s2_cherry_index].1.contains(&label1_current as &str) {
                    // if it does, then we have found 1/2 of the needed pairings
                    // next, we do the same thing for label2
                    for label2_current in s1[s1_cherry_index].1.iter() {
                        // check to see if s2 label2 has the same corresponding label
                        if s2[s2_cherry_index].0.contains(&label2_current as &str) {
                            // We found a match!
                            found_one = true;

                            // now we need to make a copy, remove all other labels, and then make a recursive call
                            let mut s1_copy = s1.clone();
                            let mut s2_copy = s2.clone();

                            // call helper
                            remove_all_other_labels(
                                &mut s1_copy,
                                &mut s2_copy,
                                label1_current.clone(),
                                label2_current.clone(),
                            );

                            // make recursive call to move onto the next cherry
                            find_single_label_cps_distance_from_mutli_label_recursive(
                                s1_copy,
                                s2_copy,
                                end_index + 1,
                                result,
                            );
                        }
                    }
                }
            }

            // if we did not find any matches then this is as long as the tail distance is going to get and we add it to the list
            if !found_one {
                // if so, enter the tail size and end
                result.insert(end_index - 1);
            }
        }
    }
}

/**
 * remove_all_other_labels
 *
 * Given two multi label CPSs for two networks and two labels, goes through each label list in each cherry in both networks and checks to see if that label list contains either label, if it does then remove all other labels except for that one.
 * @param &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s1 - Multi Label CPS 1
 * @param &Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)> s2 - Multi Label CPS 2
 * @param String label1 - First label to search for
 * @param String label2 - Second label to search for
 */
fn remove_all_other_labels(
    s1: &mut Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    s2: &mut Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>,
    label1: String,
    label2: String,
) {
    // if a label list contains label1, then remove all other labels
    // if a label list contains label2, then remove all other labels

    // check everything in s1
    for current in s1.iter_mut() {
        if current.0.contains(&label1) {
            current.0.retain(|k| *k == label1);
        } else if current.1.contains(&label1) {
            current.1.retain(|k| *k == label1);
        } else if current.0.contains(&label2) {
            current.0.retain(|k| *k == label2);
        } else if current.1.contains(&label2) {
            current.1.retain(|k| *k == label2);
        }
    }

    // check everything in s2
    for current in s2.iter_mut() {
        if current.0.contains(&label1) {
            current.0.retain(|k| *k == label1);
        } else if current.1.contains(&label1) {
            current.1.retain(|k| *k == label1);
        } else if current.0.contains(&label2) {
            current.0.retain(|k| *k == label2);
        } else if current.1.contains(&label2) {
            current.1.retain(|k| *k == label2);
        }
    }
}
