//! Distance Submodule (of Network)
//!
//! The reconstruction, deconstruction, and tail distance are rearrangement-
//! based network distances that are equal to eachother. Together will be
//! referred to as "cherry distance"
//!
//! `distance` submodule contains all functionality required to calculcate
//! the cherry distance on input networks.
//!
//!
use std::cmp;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;
use std::iter::zip;
use std::time::Instant;

use crate::network::Adj;
use crate::network::Network;
use crate::network::Node;

#[derive(Clone)]
pub struct DTree {
    //TODO pub for testing
    tree: Network<Node<usize>>, //nodes are labelled with corresponding r val
    id_i_map: BTreeMap<usize, usize>, //(k,v) = (ret id in net, node i in tree)
}
// --------------- find reduced set
//for a given network, find the set of ret-trimmed subnetworks
//i.e. the reduced set
fn ƒ<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    net: Network<Node<T>>,
    debug: bool,
) -> Vec<Network<Node<T>>> {
    //returns set of all reticulation-trimmed subnetworks
    let mut result = vec![];
    let ret_edge_pairs = net.get_retic_edge_pairs();
    let f_set = find_f_set(ret_edge_pairs);
    for mut f in f_set {
        if let Some(rt_subnet) = make_rt_subnet(&net, &mut f, debug) {
            result.push(rt_subnet);
        }
    }
    result
}
// --------------- rt subnet maker
//given an F set and network, makes the ret-trimmed subnet
fn make_rt_subnet<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    net: &Network<Node<T>>,
    f: &mut Vec<(usize, usize)>,
    debug: bool,
) -> Option<Network<Node<T>>> {
    //assumes f is disjoint
    let result;
    let mut possible = true;
    let mut new_net = net.clone();
    let f = bu_top_sort(net, f);
    let mut start = Instant::now();
    let mut elapsed;
    for (u, v) in &f {
        if new_net.is_bicon_root_of(*u, *v) {
            possible = false;
            break;
        }
        let mut rets_below: Vec<usize> = new_net
            .get_rets_below(*u)
            .union(&new_net.get_rets_below(*v))
            .cloned()
            .collect();
        rets_below.retain(|x| x != v); //we know v is a ret (below itself)
        elapsed = start.elapsed();
        //if debug {
        //  println!("{:?} all rets below",elapsed);
        //}
        start = Instant::now();
        if rets_below.is_empty() {
            //remove retic edge
            new_net.remove_edge(*u, *v);
            //update u's cluster
            let v_cluster = &new_net.get_cluster(*v);
            new_net.bubble_cluster_up_remove(*u, &v_cluster);
            new_net.bubble_cluster_up(*v, &v_cluster);
            //name parents
            let u_p = new_net.get_parents_i(*u)[0];
            let v_p = new_net.get_parents_i(*v)[0];
            //save old clusters
            let u_cluster = new_net.get_cluster(*u);
            let v_cluster = new_net.get_cluster(*v);
            elapsed = start.elapsed();
            //if debug {
            //  println!("{:?} up to remove subnets",elapsed);
            //}
            start = Instant::now();
            //remove subnets
            new_net.remove_subnet(*u);
            new_net.remove_subnet(*v);
            //update other network structs
            new_net.remove_addn_edge(*v);
            new_net.remove_nontriv_bicon(*v);
            //add new leaves with such clusters
            let u_i = new_net.get_open_id();
            new_net.add_multi_leaf(u_i, u_cluster);
            let v_i = new_net.get_open_id();
            new_net.add_multi_leaf(v_i, v_cluster);
            //attach leaves
            new_net.add_edge(u_p, u_i);
            new_net.add_edge(v_p, v_i);
            //update clusters
        } else {
            // if debug {
            //   println!("rets below empty");
            //}
            possible = false;
            break;
        }
    }
    elapsed = start.elapsed();
    //if debug {
    //  println!("{:?} to make network changes",elapsed);
    //}
    start = Instant::now();
    //--shrink reindex
    //if possible {
    //  let new_new_net = shrink_reindex(&new_net);
    //result = Some(new_new_net);
    //elapsed = start.elapsed();
    //println!("{:?} to shrink/reindex",elapsed);
    //} else {
    //  result = None;
    //}
    //return without it
    //--dont reindex
    if possible {
        result = Some(new_net);
    } else {
        result = None;
    }
    //--
    result
}
fn bu_top_sort<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    net: &Network<Node<T>>,
    f: &mut Vec<(usize, usize)>,
) -> Vec<(usize, usize)> {
    let mut result = vec![];
    while f.len() > 0 {
        let mut index = 0;
        let (a, b) = f.pop().unwrap();
        for (r_a, _) in &result {
            if net.is_below(*r_a, b) || net.is_below(*r_a, a) {
                break;
            } else {
                index += 1;
            }
        }
        result.insert(index, (a, b))
    }
    result.reverse();
    result
}
///this method shrinks the network structure, assiming many "None" values are
/// introduced when removing subnetworks  
/// used after trimming networks
fn shrink_reindex<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    net: &Network<Node<T>>,
) -> Network<Node<T>> {
    //creat reference mapping of old id to new id
    let mut node_map = BTreeMap::new(); //(k,v) = (old_i, new_i)
    let mut i = 0;
    for entry in &net.nodes {
        if let Some(node) = entry {
            //root is still 0
            node_map.insert(node.get_id(), i);
            i += 1;
        }
    }
    //duplicate nodes
    let mut new_nodes = vec![];
    for (old_i, new_i) in &node_map {
        //iterates in order (k/v in same order)
        let old_node = net.get_node_ref(*old_i);
        let new_cluster = old_node.cluster.clone();
        let new_labels = old_node.labels.clone();
        let mut new_adjs = vec![];
        for adj in &old_node.adjs {
            new_adjs.push(Adj {
                is_to: adj.is_to,
                i: node_map.get(&adj.i).unwrap().clone(),
            });
        }
        let new_node = Node::new(*new_i, new_cluster, new_labels, new_adjs);
        new_nodes.push(Some(new_node)); //pushes in order of i!
    }
    //duplicate other data
    let mut new_addn_edges = vec![];
    for (old_a, old_b) in &net.addn_edges {
        new_addn_edges.push((
            node_map.get(old_a).unwrap().clone(),
            node_map.get(old_b).unwrap().clone(),
        ));
    }
    let mut new_nontriv_bicons: HashMap<usize, usize> = HashMap::new();
    for (old_v, bicon_i) in &net.nontriv_bicons {
        new_nontriv_bicons.insert(node_map.get(old_v).unwrap().clone(), *bicon_i);
    }
    //create and return the new network
    Network {
        nodes: new_nodes,
        addn_edges: new_addn_edges,
        nontriv_bicons: new_nontriv_bicons,
    }
}
// ------------------- find f setss
fn find_f_set(s: Vec<Vec<(usize, usize)>>) -> Vec<Vec<(usize, usize)>> {
    //assume that s is a set of pairs
    let mut results = vec![];
    for set in one_of_each(s) {
        let powerset = powerset(set);
        for subset in powerset {
            results.push(subset);
        }
    }
    results.sort();
    results.dedup();
    results
}
fn powerset(s: Vec<(usize, usize)>) -> Vec<Vec<(usize, usize)>> {
    (0..2usize.pow(s.len() as u32))
        .map(|i| {
            s.iter()
                .enumerate()
                .filter(|&(t, _)| (i >> t) % 2 == 1)
                .map(|(_, element)| *element)
                .collect()
        })
        .collect()
}
fn one_of_each(matched_pairs: Vec<Vec<(usize, usize)>>) -> Vec<Vec<(usize, usize)>> {
    //generate all binary numbers of length n
    let n = matched_pairs.len();
    let m = 2_i32.pow(n as u32);
    let mut bin_strings = vec![];
    let mut result = vec![];
    if m > 1 {
        for x in 0..m {
            bin_strings.push(format!("{:0n$b}", x));
        }
    } else {
        result.push(vec![]);
    }
    //for each binary string, use chars to choose matched_pair
    for s in bin_strings {
        let mut temp = vec![];
        for (i, c) in s.chars().enumerate() {
            if c == '0' {
                temp.push(matched_pairs[i][0]);
            } else {
                temp.push(matched_pairs[i][1]);
            }
        }
        result.push(temp);
    }
    result
}
// ------------------- find simple distance DP DP
fn find_simple_macrs_size<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    n1: &Network<Node<T>>,
    n2: &Network<Node<T>>,
) -> Option<usize> {
    //retuns Some<dist> with distance being the SIMPLE distance between inputs
    //returns None if simple distance cannot be found
    //lets think of M as a function rather than a structure!
    //we do not need to fill all of the table M,
    //table M returns size of macrs, use that calc distance
    //WANT size of simple network, not simple dist!
    if let Some(macrs_leaves) = table_m(n1, n2, 0, 0) {
        let macrs_size = 2 * macrs_leaves + 2 * n1.addn_edges.len() - 1;
        Some(macrs_size)
        //Some( n1.get_n() + n2.get_n() - 2*macrs_size //2)//gives dist not size
    } else {
        None
    }
}
fn table_m<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
    n1: &Network<Node<T>>,
    n2: &Network<Node<T>>,
    u: usize,
    v: usize,
) -> Option<usize> {
    let is_rets_below_union_empty = n1
        .get_rets_below(u)
        .union(&n2.get_rets_below(v))
        .cloned()
        .collect::<Vec<usize>>()
        .is_empty();
    if n1.is_trivial(u) && n2.is_trivial(v) {
        if n1.is_leaf_node(u) || n2.is_leaf_node(v) {
            if (!n1
                .get_cluster(u)
                .intersection(&n2.get_cluster(v))
                .cloned()
                .collect::<Vec<T>>()
                .is_empty())
                && is_rets_below_union_empty
            {
                Some(1)
            } else {
                None
            }
        } else {
            let u1 = n1.get_children_i(u)[0];
            let u1c = n1.get_cluster(u1);
            let u2 = n1.get_children_i(u)[1];
            let u2c = n1.get_cluster(u2);
            let v1 = n2.get_children_i(v)[0];
            let v1c = n2.get_cluster(v1);
            let v2 = n2.get_children_i(v)[1];
            let v2c = n2.get_cluster(v2);
            let x11: Vec<T> = u1c.intersection(&v1c).cloned().collect();
            let x12: Vec<T> = u1c.intersection(&v2c).cloned().collect();
            let x21: Vec<T> = u2c.intersection(&v1c).cloned().collect();
            let x22: Vec<T> = u2c.intersection(&v2c).cloned().collect();
            let m1;
            let m2;
            // -- m1
            if is_rets_below_union_empty
                && ((!x11.is_empty() && x22.is_empty()) || (x11.is_empty() && !x22.is_empty()))
            {
                //condition 1 and 2
                m1 = Some(1);
            } else if !x11.is_empty() && !x22.is_empty() {
                //condition 3: recursive "joining" step
                if let Some(l) = table_m(n1, n2, u1, v1) {
                    if let Some(r) = table_m(n1, n2, u2, v2) {
                        m1 = Some(l + r);
                    } else {
                        m1 = None; //if reticulations below, no macrs possible
                    }
                } else {
                    m1 = None; //if reticulations below, no macrs possible
                }
            } else {
                //there are rets below or both are empty
                m1 = None; //the "otherwise" condition
            }
            // -- m2
            if is_rets_below_union_empty
                && ((!x12.is_empty() && x21.is_empty()) || (x12.is_empty() && !x21.is_empty()))
            {
                //condition 1 and 2
                m2 = Some(1);
            } else if !x12.is_empty() && !x21.is_empty() {
                //condition 3: recursive "joining" step
                if let Some(l) = table_m(n1, n2, u1, v2) {
                    if let Some(r) = table_m(n1, n2, u2, v1) {
                        m2 = Some(l + r);
                    } else {
                        m2 = None; //if reticulations below, no macrs possible
                    }
                } else {
                    m2 = None; //if reticulations below, no macrs possible
                }
            } else {
                //there are rets below or both are empty
                m2 = None; //the "otherwise" condition
            }
            // max of m1/m2
            if let Some(m1_val) = m1 {
                if let Some(m2_val) = m2 {
                    Some(cmp::max(m1_val, m2_val))
                } else {
                    //m1 is some but not m2
                    m1
                }
            } else {
                if let Some(_) = m2 {
                    //m2 is some but not m1
                    m2
                } else {
                    //neither are Some
                    None
                }
            }
        }
    } else if (n1.is_trivial(u) && !n2.is_trivial(v)) || (!n1.is_trivial(u) && n2.is_trivial(v)) {
        None
    } else {
        //both u, v are non-trivial
        //conjecture: u and v are the roots of the nontrivial component
        let r1 = n1.get_ret_of_bicon(u);
        let r2 = n2.get_ret_of_bicon(v);
        //paths do not include u , v, or r
        let (pi_1l, pi_1r) = n1.get_component_paths(u, r1);
        let (pi_2l, pi_2r) = n2.get_component_paths(v, r2);
        let mut m1 = None;
        let mut m2 = None;
        if let Some(ret_pair) = table_m(n1, n2, n1.get_children_i(r1)[0], n2.get_children_i(r2)[0])
        {
            // try for m1
            if pi_1l.len() == pi_2l.len() && pi_1r.len() == pi_2r.len() {
                let mut m1_possible = true;
                let mut l_total = ret_pair;
                let mut r_total = 0;
                'outer: for (one, two) in zip(&pi_1l, &pi_2l) {
                    let mut h_one = n1.get_children_i(*one)[0];
                    if pi_1l.contains(&h_one) || h_one == r1 {
                        h_one = n1.get_children_i(*one)[1];
                    }
                    let mut h_two = n2.get_children_i(*two)[0];
                    if pi_2l.contains(&h_two) || h_two == r2 {
                        h_two = n2.get_children_i(*two)[1];
                    }
                    if let Some(l_match) = table_m(n1, n2, h_one, h_two) {
                        l_total += l_match;
                    } else {
                        m1_possible = false;
                        break 'outer;
                    }
                }
                if m1_possible {
                    'outer: for (one, two) in zip(&pi_1r, &pi_2r) {
                        let mut h_one = n1.get_children_i(*one)[0];
                        if pi_1r.contains(&h_one) || h_one == r1 {
                            h_one = n1.get_children_i(*one)[1];
                        }
                        let mut h_two = n2.get_children_i(*two)[0];
                        if pi_2r.contains(&h_two) || h_two == r2 {
                            h_two = n2.get_children_i(*two)[1];
                        }
                        if let Some(r_match) = table_m(n1, n2, h_one, h_two) {
                            r_total += r_match;
                        } else {
                            m1_possible = false;
                            break 'outer;
                        }
                    }
                }
                if m1_possible {
                    m1 = Some(l_total + r_total);
                }
            }
            //try for m2
            if pi_1l.len() == pi_2r.len() && pi_1r.len() == pi_2l.len() {
                let mut m2_possible = true;
                let mut l_total = ret_pair;
                let mut r_total = 0;
                'outer: for (one, two) in zip(&pi_1l, &pi_2r) {
                    let mut h_one = n1.get_children_i(*one)[0];
                    if pi_1l.contains(&h_one) || h_one == r1 {
                        h_one = n1.get_children_i(*one)[1];
                    }
                    let mut h_two = n2.get_children_i(*two)[0];
                    if pi_2r.contains(&h_two) || h_two == r2 {
                        h_two = n2.get_children_i(*two)[1];
                    }
                    if let Some(l_match) = table_m(n1, n2, h_one, h_two) {
                        l_total += l_match;
                    } else {
                        m2_possible = false;
                        break 'outer;
                    }
                }
                if m2_possible {
                    'outer: for (one, two) in zip(&pi_1r, &pi_2l) {
                        let mut h_one = n1.get_children_i(*one)[0];
                        if pi_1r.contains(&h_one) || h_one == r1 {
                            h_one = n1.get_children_i(*one)[1];
                        }
                        let mut h_two = n2.get_children_i(*two)[0];
                        if pi_2l.contains(&h_two) || h_two == r2 {
                            h_two = n2.get_children_i(*two)[1];
                        }
                        if let Some(r_match) = table_m(n1, n2, h_one, h_two) {
                            r_total += r_match;
                        } else {
                            m2_possible = false;
                            break 'outer;
                        }
                    }
                }
                if m2_possible {
                    m2 = Some(l_total + r_total);
                }
            }
        }
        //note: ret pair was added onto l total
        //dealing with m1 and m2 results
        if let Some(m1_val) = m1 {
            if let Some(m2_val) = m2 {
                //both m1 and m2 were some
                Some(cmp::max(m1_val, m2_val))
            } else {
                //m1 some not m2
                m1
            }
        } else if let Some(_) = m2 {
            //m2 some, not m1
            m2
        } else {
            //neither m1 or m2 were some
            None
        }
    }
}

// -------------------- DTree specific implementation
impl DTree {
    //------------- make dependency trees
    pub fn dep_tree_maker<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
        n: &Network<Node<T>>,
    ) -> Self {
        let mut new_id_i_map: BTreeMap<usize, usize> = BTreeMap::new();
        let mut root_label = BTreeSet::new();
        root_label.insert(0);
        new_id_i_map.insert(0, 0);
        let mut new_nodes = vec![Some(Node {
            i: 0,
            cluster: BTreeSet::new(),
            labels: root_label,
            adjs: vec![],
        })];
        for (i, r) in n.get_ret_nodes().iter().enumerate() {
            new_id_i_map.insert(*r, i + 1);
            let mut new_label: BTreeSet<usize> = BTreeSet::new();
            new_label.insert(*r);
            new_nodes.push(Some(Node {
                i: i + 1,
                cluster: BTreeSet::new(),
                labels: new_label,
                adjs: vec![],
            }));
        }
        let mut new_tree = Network {
            nodes: new_nodes,
            addn_edges: vec![],
            nontriv_bicons: HashMap::new(),
        };
        //add edges
        for i in 0..new_tree.nodes.len() {
            let binding = new_tree.get_labels(i);
            let node_id = binding.first().unwrap();
            let children = Self::get_first_rets_below(n, *node_id);
            for child in children {
                let child_i = new_id_i_map.get(&child).unwrap();
                new_tree.add_edge(i, *child_i);
            }
        }
        DTree {
            tree: new_tree,
            id_i_map: new_id_i_map,
        }
    }
    fn get_first_rets_below<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
        n: &Network<Node<T>>,
        cur: usize,
    ) -> Vec<usize> {
        //get one or get the set?
        //TODO the sort and dedup may be able to be solved by not traversing
        //addn_edges: will that leave out rets?
        //conjecture: if ret can be found by addn_edge only, then not level-1
        let mut next_curs = n.get_children_i(cur);
        let mut result = vec![];
        while next_curs.len() > 0 {
            let next_cur = next_curs.pop().unwrap();
            if n.is_ret_node(next_cur) {
                //push ret to result
                result.push(next_cur);
            } else {
                //keep adding children
                for next_child in n.get_children_i(next_cur) {
                    next_curs.push(next_child);
                }
            }
        }
        result.sort();
        result.dedup();
        result
    }
    //-------------------------smart enum
    // ---------------enum
    pub fn make_pairs<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
        t1: &DTree,
        n1: &Network<Node<T>>,
        t2: &DTree,
        n2: &Network<Node<T>>,
    ) -> Vec<(DTree, DTree)> {
        //special use of labels and map in self.
        //self labelled with corresponding index in t1.
        //self map is (k,v)=(i in self, index of t2)
        //---
        //make base t'
        let mut root_label = BTreeSet::new();
        root_label.insert(0);
        let root = Node {
            i: 0,
            cluster: BTreeSet::new(),
            labels: root_label,
            adjs: vec![],
        };
        let t_tree = Network {
            nodes: vec![Some(root)],
            addn_edges: vec![],
            nontriv_bicons: HashMap::new(),
        };
        let mut t_map: BTreeMap<usize, usize> = BTreeMap::new();
        t_map.insert(0, 0);
        let t = DTree {
            tree: t_tree,
            id_i_map: t_map,
        };
        //
        let min_v = cmp::min(t1.tree.get_n(), t2.tree.get_n());
        let mut ts = vec![t.clone()];
        let mut prev_lvl = vec![t.clone()];
        //a level is all trees of the same size
        for i in 1..min_v {
            let mut next_lvl = vec![];
            for t in prev_lvl {
                //for each unique t
                for v_node in &t.tree.nodes {
                    let v = v_node.as_ref().unwrap();
                    let v_chs = t.tree.get_children_i(v.i);
                    //for each node in v
                    let vt1 = t.tree.get_label(v.i);
                    let vt2 = t.id_i_map.get(&v.i).unwrap();
                    let mut chs_vt1 = t1.tree.get_children_i(vt1);
                    //filter out so we look only at children of v not in t
                    chs_vt1.retain(|ch| v_chs.iter().all(|v_ch| *ch != t.tree.get_label(*v_ch))); //the value is not a label of children
                    let mut chs_vt2 = t2.tree.get_children_i(*vt2);
                    chs_vt2.retain(|ch| t.id_i_map.values().all(|v| ch != v)); //the value is not in the map as a value
                    if chs_vt1.len() > 0 && chs_vt2.len() > 0 {
                        for ch_vt1 in chs_vt1 {
                            for ch_vt2 in &chs_vt2 {
                                let mut new_t = t.clone();
                                let new_i = new_t.tree.get_open_id();
                                new_t.tree.add_leaf(new_i, ch_vt1);
                                new_t.tree.add_edge(v.i, new_i);
                                new_t.id_i_map.insert(new_i, *ch_vt2);
                                //next_lvl.push(new_t);
                                if next_lvl.iter().all(|t: &DTree| !t.is_cong(&new_t)) {
                                    //only push if this new tree is unique
                                    next_lvl.push(new_t);
                                }
                            }
                        }
                    }
                }
            }
            prev_lvl = next_lvl.clone();
            ts.append(&mut next_lvl);
        }
        let mut results = vec![];
        // -------------- make the pairs
        for t in ts {
            //make a pair from the t'
            let mut tree1: Network<Node<usize>> = Network::new(1);
            let mut id_i_map_1: BTreeMap<usize, usize> = BTreeMap::new();
            let mut tree2: Network<Node<usize>> = Network::new(1);
            let mut id_i_map_2: BTreeMap<usize, usize> = BTreeMap::new();
            //add root
            tree1.add_leaf(0, 0);
            id_i_map_1.insert(0, 0);
            tree2.add_leaf(0, 0);
            id_i_map_2.insert(0, 0);
            for v_node in &t.tree.nodes {
                let v = v_node.as_ref().unwrap();
                let mut v_1 = t.tree.get_label(v.i);
                let mut v_2 = t.id_i_map.get(&v.i).unwrap();
                //for each vertex, add all its children
                let v_chs = t.tree.get_children_i(v.i);
                for v_ch in v_chs {
                    //tree1
                    //let new_id_1 = tree1.get_open_id();
                    let v_ch_1 = t.tree.get_label(v_ch);
                    let t1_ret = t1.tree.get_label(v_ch_1);
                    tree1.add_leaf(v_ch_1, t1_ret);
                    tree1.add_edge(v_1, v_ch_1);
                    id_i_map_1.insert(t1_ret, v_ch_1);
                    //tree2
                    //let new_id_2 = tree2.get_open_id();
                    let v_ch_2 = t.id_i_map.get(&v_ch).unwrap();
                    let t2_ret = t2.tree.get_label(*v_ch_2);
                    tree2.add_leaf(*v_ch_2, t2_ret);
                    tree2.add_edge(*v_2, *v_ch_2);
                    id_i_map_2.insert(t2_ret, *v_ch_2);
                }
            }
            let mut new_t1 = DTree {
                tree: tree1,
                id_i_map: id_i_map_1,
            };
            let mut new_t2 = DTree {
                tree: tree2,
                id_i_map: id_i_map_2,
            };
            results.push((new_t1, new_t2));
        }
        results
    }
    // ----------------- is iso
    fn is_cong(&self, n: &Self) -> bool {
        //just has to check all the same values and labels , NOT shape!
        self.id_i_map
            .values()
            .all(|s_v| n.id_i_map.values().any(|n_v| *s_v == *n_v))
            && self
                .tree
                .find_cluster()
                .iter()
                .all(|s_l| n.tree.find_cluster().iter().any(|n_l| *s_l == *n_l))
    }
    //----------------embeddings
    //pub fn all_embeddings_of(&self, gt:DTree) -> Vec<DTree> {
    //  let mut wips:Vec<DTree> = vec![];
    //build base DTree
    //let mut new_id_i_map:HashMap<usize,usize> = HashMap::new();
    //let mut root_label = BTreeSet::new();
    //root_label.insert(0);
    //new_id_i_map.insert(0,0);
    //let mut new_nodes = vec![ Some( Node{
    //                          i: 0,
    //                        cluster: BTreeSet::new(),
    //                      labels: root_label,
    //                    adjs: vec![],
    //          } )
    //];
    ///let mut new_tree = Network{
    // nodes: new_nodes,
    //addn_edges: vec![],
    //nontriv_bicons: HashMap::new(),
    //};
    //let new_dtree = DTree{id_i_map:new_id_i_map, tree:new_tree};
    //TODO this will map all the subtrees too along the way???
    //can I use this to my advantage?
    //wips.push(new_dtree);
    //for _ in 0..gt.tree.nodes.len() {
    //  let mut next_wips = vec![];
    //}
    //}
    //pub fn all_embeddings_of_helper(&self, s_cur:usize, gt: DTree, gt_cur:usize, cur_tree: DTree) -> Vec<Option<DTree>> {
    //  let gt_children = gt.tree.get_children_i(gt_cur);
    //let self_children = self.tree.get_children_i(s_cur);
    //let gt_child_count = gt_children.len();
    //let self_child_count = self_children.len();
    //if self_child_count < gt_child_count {
    //what is none?
    //  None//cannot make tree with these curs matching
    //}
    //for self_child in self_children {
    //  for gt_child in gt_children {
    //    let label = self.get_id_from_i(self_child);
    //  let new_cur_tree = cur_tree.clone();
    //tree.add_leaf(self_child, label);
    //new_net
    // }
    // }
    //}
    //pub fn all_noniso_embeddings_of(&self, dt: &DTree) -> Vec<(DTree,DTree)> {
    //pairs of t'1 (self), t'2(dt)
    //without actual isomorphisms, but rather all "same shapes"
    //  let mut results = vec![];
    //preprocess to determine "terminal vertices", all children mutual iso

    //}
    // --------------- generic build for testing
    pub fn build_from_level_seq(seq: Vec<usize>) -> Self {
        let mut new_nodes = vec![];
        let mut new_id_i_map = BTreeMap::new();
        //for each digit in seq create a node
        for (index, _) in seq.iter().enumerate() {
            let mut new_label = BTreeSet::new();
            new_label.insert(index); //generic label equal to index
            new_nodes.push(Some(Node {
                i: index,
                cluster: BTreeSet::new(),
                labels: new_label,
                adjs: vec![],
            }));
            new_id_i_map.insert(index, index);
        }
        let mut new_tree = Network {
            nodes: new_nodes,
            addn_edges: vec![],
            nontriv_bicons: HashMap::new(),
        };
        //create and update adjs
        let mut parent = 0;
        let mut prev_i = 0;
        for (index, i) in seq.iter().enumerate() {
            if *i > prev_i {
                parent = index - 1;
            } else if *i < prev_i {
                //go backwards from seq from i until i-1 found, thats the parent
                let mut parent_found = false;
                let mut j_index = index - 1;
                while !parent_found {
                    if seq[j_index] == i - 1 {
                        parent = j_index;
                        parent_found = true;
                    } else {
                        j_index -= 1;
                    }
                }
            } //else parent stays the same
            if index != 0 {
                //prevents edge (0,0)
                new_tree.add_edge(parent, index)
            }
            prev_i = *i;
        }
        //update tree node type
        for i in 0..new_tree.nodes.len() {
            if new_tree.has_no_children(i) {
                let mut node = new_tree.get_node_mut(i);
            }
        }
        //return the results
        DTree {
            tree: new_tree,
            id_i_map: new_id_i_map,
        }
    }
    // --------------- trim dependency trees
    pub fn remove_i_generic(&self, i: usize) -> Vec<DTree> {
        //checks all the way to remove i up to topological isomorphism
        //labels dont matter and are not taken into account.
        if i >= self.tree.nodes.len() {
            panic!("cannot reduce past singleton");
        }
        //(removed leaves so far, the resulting tree)
        let mut wips: Vec<DTree> = vec![];
        wips.push(self.clone());
        for _ in 0..i {
            let mut next_wips = vec![];
            while wips.len() > 0 {
                let wip = wips.pop().unwrap();
                for leaf in wip.tree.get_leaves() {
                    let mut cur = wip.clone();
                    cur.remove_leaf(leaf);
                    if next_wips
                        .iter()
                        .all(|t: &DTree| !t.is_shape_iso(0, &cur, 0))
                    {
                        //check that cur is not shape isomorphic to any t already in next wips
                        next_wips.push(cur);
                    }
                }
            }
            wips = next_wips;
        }
        wips
    }
    pub fn remove_i_labelled(&self, i: usize) -> Vec<DTree> {
        //TODO pub for testing
        //all the ways to remove i from a labelled network, considers labels
        if i >= self.tree.nodes.len() {
            panic!("cannot reduce past singleton");
        }
        //(removed leaves so far, the resulting tree)
        let mut wips: Vec<(BTreeSet<usize>, DTree)> = vec![];
        wips.push((BTreeSet::new(), self.clone()));
        for _ in 0..i {
            let mut next_wips = vec![];
            while wips.len() > 0 {
                let wip = wips.pop().unwrap();
                for leaf in wip.1.tree.get_leaves() {
                    let mut cur = wip.clone();
                    cur.0.insert(leaf);
                    if next_wips.iter().all(|(s, _): &(BTreeSet<usize>, DTree)| {
                        !(s.is_superset(&cur.0) && s.is_subset(&cur.0))
                    }) {
                        //check that some other wip hasn't removed the same leaf set
                        cur.1.remove_leaf(leaf);
                        next_wips.push(cur);
                    }
                }
            }
            wips = next_wips;
        }
        let (_, results): (Vec<BTreeSet<usize>>, Vec<DTree>) = wips.into_iter().unzip();
        results
    }
    // -------------------- make rt subnets from sub dep tree
    pub fn rt_from_subd_maker<T: Clone + PartialEq + Eq + Hash + Ord + Display>(
        &self,
        dt: &DTree,
        n: &Network<Node<T>>,
        debug: bool,
    ) -> Vec<Network<Node<T>>> {
        //for some subdt self, make all rt_subnets and populate self.
        //start timer
        let mut start = Instant::now();
        let mut elapsed = start.elapsed();
        //returns the set of all rt subnetworks given a subdependency tree
        let mut req_rs = vec![]; //required rets (to be reduced) i.e. the first absent layer of rets below kept layer
        for i in 0..self.tree.nodes.len() {
            //go by labels and NOT indices
            if self.tree.nodes[i].is_some() {
                let i_label = self.tree.get_label(i);
                //how on earth would a label in subdt NOT be in dt?
                if let Some(dt_i) = dt.id_i_map.get(&i_label) {
                    for dt_child in dt.tree.get_children_i(*dt_i) {
                        let dt_child_label = dt.tree.get_label(dt_child);
                        //for the subdt node and its corresponding dt node
                        //compare childs. if absent in subdt...
                        if self.id_i_map.get(&dt_child_label).is_none() {
                            req_rs.push(dt_child_label);
                        }
                    }
                }
            }
        }
        let mut req_pairset = vec![];
        for r in &req_rs {
            let mut pair = vec![];
            for parent in n.get_parents_i(*r) {
                pair.push((parent, *r));
            }
            req_pairset.push(pair);
        }
        //make every iteration of the req f set
        let mut req_fsets = one_of_each(req_pairset);
        //make trivial f set
        let mut triv_fset = vec![];
        for r in req_rs {
            //println!("req r: {r}");
            let triv_rs = n.get_rets_below(r);
            for triv_r in &triv_rs {
                let ps = n.get_parents_i(*triv_r);
                let mut p = None;
                for cand_p in &ps {
                    if ps
                        .iter()
                        .any(|other_p| cand_p != other_p && !n.is_below(*other_p, *cand_p))
                    {
                        p = Some(*cand_p);
                    }
                }
                if p.is_some() {
                    triv_fset.push((p.unwrap(), *triv_r));
                } else {
                    panic!("there were no available edges to pick for {triv_r}");
                }
            }
        }
        //for each iteration, add the trivial f set
        req_fsets
            .iter_mut()
            .for_each(|mut x| x.extend(triv_fset.clone()));
        //make rt subnets and populate subdt
        if debug {
            start = Instant::now();
        }
        let mut n_nets = 0;
        let mut results = vec![];
        for set in req_fsets {
            //println!("{:?}",set);
            //make the subnet, set is good so must return Some
            //TODO set is supposed to be good!!! why would it return None?
            if let Some(cand_net) = make_rt_subnet(n, &mut set.clone(), debug) {
                results.push(cand_net);
                n_nets += 1;
            }
        }
        if debug {
            elapsed = start.elapsed();
            //println!("{:?} to make {n_nets} nets",elapsed);
        }
        results
    }
    // ------------------ generic helpers
    fn make_key(&self) -> String {
        //keys are returned in order
        //keys are the ret id in net.
        self.id_i_map
            .keys()
            .cloned()
            .map(|i| i.to_string())
            .collect::<String>()
    }
    fn remove_leaf(&mut self, leaf_i: usize) -> usize {
        let binding = self.tree.get_labels(leaf_i);
        let removed_id = binding.first().unwrap();
        self.id_i_map.remove(removed_id);
        self.tree.remove_leaf(leaf_i);
        *removed_id
    }
    fn is_shape_iso(&self, self_cur: usize, t: &DTree, t_cur: usize) -> bool {
        let self_child_count = self.tree.get_children_i(self_cur).len();
        let t_child_count = t.tree.get_children_i(t_cur).len();
        if self_child_count == 0 && t_child_count == 0 {
            //base case both are leaves
            true
        } else if self_child_count == 0 || t_child_count == 0 {
            //base fail one is leaf and other is not
            false
        } else if self_child_count == t_child_count {
            //recursive case true
            let s_childs = self.tree.get_children_i(self_cur);
            let t_childs = t.tree.get_children_i(t_cur);
            let mut both = s_childs.iter().zip(t_childs.iter());
            both.all(|(s_c, t_c)| self.is_shape_iso(*s_c, t, *t_c))
        } else {
            //recursive case fail
            false
        }
    }
    fn get_id_from_i(&self, i: usize) -> Option<usize> {
        self.id_i_map
            .iter()
            .find_map(|(key, &val)| if val == i { Some(key) } else { None })
            .copied()
    }
}
impl std::fmt::Display for DTree {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut result = String::new();
        result = self.tree.to_string();
        write!(f, "{}", result)
    }
}
