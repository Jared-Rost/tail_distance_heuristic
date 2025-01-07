//! Random Submodule (of Network )
//!
//! Generates a random network topology that is level-1 and orchard
//! 
//! `random` submodule contains all functionality associated with generating 
//! random trees and for attaching additional edges to form reticulations.
//! 

use rand::{distributions::Uniform, Rng};
use std::cmp;
//
use crate::network::Network;
use crate::network::Node;
use crate::network::gen_label;

// --------------------------- CONSTRUCTOR
pub fn new_random(size: usize, r: usize, exact: bool) -> Network<Node<String>> {
    const MAX_TRY_FOR_EXACT: i32 = 100;
    if size < 2 {panic!("cannot construct network with less than 2 leaves")}
    //NOTE:check r isnt too big
    //size is # leaves
    //v is # nodes
    let v = 2*size-1;
    //generate the random numbers needed
    let mut rands: Vec<usize> = vec![]; 
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(1usize..100usize);//[1,99] (percent)
    for _ in 0..v {
        rands.push(rng.sample(dist));
    } 
    //------------------------build tree and add reticulations
    let mut result: Network<Node<String>> = Network::new(size);
    if exact {
        //build tree
        for i in 0..MAX_TRY_FOR_EXACT {
            if i == MAX_TRY_FOR_EXACT-1 {
                panic!("could not add exact number of reticulations in {} tries", MAX_TRY_FOR_EXACT);
            }
            let mut net = Network::new(size);
            net = new_random_tree(&rands, net, 0, size, 0);
            net.update_tree_clusters();
            //add rets
            result = add_random_rets(net, r);
            if result.get_ret_nodes().len() == r {
                break
            }
        }
    } else {
        //build tree
        let mut net = new_random_tree(&rands, result, 0, size, 0);
        net.update_tree_clusters();
        //add rets
        result = add_random_rets(net, r);
    }
    result
}
fn new_random_tree(
        cur_rands: &[usize],
        mut cur_net: Network<Node<String>>,
        cur_index: usize,
        cur_size: usize,//# of leaves
        par_i: usize,
    ) -> Network<Node<String>> {
    //assumes Network will not be singleton
    if cur_size < 2 {//BASE CASE
        //create leaf and fill cur_node position in Network
        cur_net.add_leaf(cur_index, gen_label());
        cur_net.add_edge(par_i, cur_index);
    } else {
        //create+add internal node at cur_node position in Network
        if cur_index == 0 {
            cur_net.add_node_i(0);
        } else {
            cur_net.add_node_i(cur_index);
            cur_net.add_edge(par_i, cur_index);
        }
        //   ---------   recursive calls for left and right subNetworks
        //determine size (# leaves) of subNetworks
        let left_pc = cur_rands[0];
        let left_size = cmp::max(1, left_pc*cur_size/100);
        let right_size = cur_size - left_size;
        let ( _ , cur_rands) = cur_rands.split_at(1);
        //left and right of rands
        let (left_rands, right_rands) = cur_rands.split_at(2*left_size-1);
        //left and right indices
        let left_index = cur_index + 1;
        let right_index = cur_index + (2*left_size-1) + 1;
        //recursive calls
        cur_net = new_random_tree(left_rands, cur_net, left_index, left_size, cur_index);
        cur_net = new_random_tree(right_rands, cur_net, right_index, right_size, cur_index);
    }
    cur_net
}
fn add_random_rets(net: Network<Node<String>>, r:usize) -> Network<Node<String>> {
    //only guarunteed to return up to r edge pairs, not exactly r
    if !net.is_tree() { panic!("reticulations cannot be added piecemeal."); }
    let potential_bicons:Vec<Vec<(usize,usize)>> = vec![net.get_edges()];
    if potential_bicons[0].len() < 2 { net } else {
        add_random_rets_helper(net,r, potential_bicons)
    }
}
fn add_random_rets_helper(
                    mut net: Network<Node<String>>, 
                    r: usize, 
                    mut potential_bicons: Vec<Vec<(usize,usize)>>, 
    ) -> Network<Node<String>> {
    if r==0 || potential_bicons.is_empty()  {
        return net;
    }
    //set up random distribution
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0usize..99usize);//cant use 100% of length as index
    //n = total across all groups. choose in proportion to group size
    let n:usize = potential_bicons.iter().map(|l| l.len()).sum();
    //set up edge selection
    let e1:(usize,usize);let e2:(usize,usize);
    let mut cur_group; let mut g = 0; //g is index of cur group
    let mut edge_meta_i = rng.sample(dist)*n/100;
    //start selection
    'edge_1: loop {
        for group in &mut potential_bicons {
            for (i, _) in group.iter().enumerate() {
                if edge_meta_i == 0 {   
                    e1 = group.swap_remove(i);
                    break 'edge_1;
                }
                edge_meta_i -= 1;
            }
            g+=1;
        }
    }
    cur_group = potential_bicons.swap_remove(g);
    //choose e2
    let edge_i = rng.sample(dist)*cur_group.len()/100;
    e2 = cur_group.swap_remove( edge_i );
    //choose direction of ret edge
    let (_, in_1) = e1; let (_, in_2) = e2; 
    let mut out_edge = e1; let mut in_edge = e2;
    if net.is_below(in_1, in_2) {
        out_edge = e2; in_edge = e1;
    }
    //add the reticulation
    let ret_node_i = net.add_reticulation(out_edge,in_edge,r); 
    let (_, bc) = in_edge;
    let (_, ac) = out_edge;
    cur_group.push((ret_node_i,bc));
    cur_group.push((ret_node_i-1,ac));
    //update potential bicons
    let bicon = net.get_bicon(ret_node_i);
    cur_group.retain(|&(a,b)| !(bicon.contains(&a) && bicon.contains(&b)));
    //var set up
    let mut new_groups:Vec<Vec<(usize,usize)>> = vec![vec![];bicon.len()];
    let mut outer_roots:Vec<Option<usize>> = vec![None;bicon.len()];
    let mut moved = vec![false;cur_group.len()];
    //find "outer roots" group
    for (b_i,v) in bicon.iter().enumerate() {
        //so inefficient
        for (cg_i,(a,b)) in cur_group.iter().enumerate() {
            if a==v {
                outer_roots[b_i] = Some(*b);
                moved[cg_i] = true;//flag for later removal from cur_group
            } 
        }
    }
    for (cg_i, (a,b)) in cur_group.iter().enumerate() {
        //add edges to a new group
        for (r_i,root) in outer_roots.iter().enumerate() {
            if root.is_some() {
                if net.is_below(*a,root.unwrap()) {
                    new_groups[r_i].push((*a,*b));
                    moved[cg_i] = true;
                } else if *b == root.unwrap() {
                    //no reticulation is an outer root so
                    //this is the rooting edge
                    new_groups[r_i].push((*a,*b));
                    moved[cg_i] = true;
                }
            }
        }
    }
    //all other edges are thier own group
    let mut moved_iter = moved.iter();
    cur_group.retain(|_| !(moved_iter.next().unwrap()) );
    new_groups.push(cur_group);
    //add new groups to potential bicons if large enough
    for new_group in new_groups {
        if new_group.len() > 1 {
            //conjecture: 2 connected edges can have ret edge between them
            potential_bicons.push(new_group);
        }
    }
    //update network, down here due to borrowing rules
    net.add_nontriv_bicon(bicon,r);
    //call again for next reticulation
    add_random_rets_helper(net, r-1, potential_bicons)
}
