//! Network Module
//!
//! Networks are binary orchards whose leaves are labelled.
//!
//! `network` module contains all network functionality including modules for
//! building and modifying networks
//!
//! This module contains 2 mains types: `Network<Node<T>>`, `Node<T>`
//!
use rand::distributions::Uniform;
use rand::Rng;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;

use crate::network::cherry::CherryType;
use crate::network::distance::DTree;

pub mod cherry;
pub mod distance;
mod heuristic;
mod newick;
mod random;

///An orchard network
pub struct Network<T> {
    /// The set of nodes, the root is always at index 0.
    /// the index is the id of node stored there
    nodes: Vec<Option<T>>,

    ///a subset of the reticulation edges (one corresponding to each reticulation), the removal of which leaves a base tree, these edges DO NOT neccessarily represent HGT.
    addn_edges: Vec<(usize, usize)>,

    ///The non-trivial biconnected components where (key,value) = (node ID, biconnected component ID)
    nontriv_bicons: HashMap<usize, usize>,
}

///A node which stores adjacencies and node type
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Node<T: Clone + PartialEq + Eq + Hash + Ord + Display> {
    ///a node id, also represents index in network node array
    i: usize,

    ///the set of labels of all leaves descending this node
    cluster: BTreeSet<T>,

    ///the labels should this node be a leaf. supports multi-labelling
    labels: BTreeSet<T>,

    ///the adjacencies including arc direction of this node
    adjs: Vec<Adj>,
}
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Adj {
    is_to: bool, //true when adj is child to v
    i: usize,    //ID of adjacent node;
}
//   ------------------------------- NETWORK Implementation
impl<T: Clone + PartialEq + Eq + Hash + Ord + Display> Network<Node<T>> {
    //TODO: keep cluster updated with edge add/delete operation??
    // ----- Constructors and Building
    ///Create a new empty network with enough initial space to construct a tree
    /// # Arguments
    /// * `size` - the number of leaves in the network
    pub fn new(size: usize) -> Self {
        //size is # leaves
        let leaves;
        if size < 2 {
            leaves = 2;
        } else {
            leaves = size;
        }
        let v = 2 * leaves - 1;
        let mut nodes = vec![];
        nodes.resize_with(v, Default::default);
        Network {
            nodes: nodes,
            addn_edges: vec![],
            nontriv_bicons: HashMap::new(),
        }
    }
    pub fn find_heuristic_sequence(
        n1: Network<Node<String>>,
        n2: Network<Node<String>>,
        num_iterations: usize,
    ) -> Vec<Vec<(String, String, CherryType)>> {
        heuristic::find_heuristic_sequence(n1, n2, num_iterations)
    }
    pub fn find_heuristic_sequence_multi_label(
        n1: Network<Node<String>>,
        n2: Network<Node<String>>,
        num_iterations: usize,
    ) -> Vec<Vec<(BTreeSet<String>, BTreeSet<String>, CherryType)>> {
        /* let leaves1 = n1.get_leaves_count();
        let leaves2 = n2.get_leaves_count();
        let rets1 = n1.get_ret_nodes().len();
        let rets2 = n2.get_ret_nodes().len();
        println!("Network 1 has {leaves1} leaves and {rets1} reticulations.");
        println!("Network 2 has {leaves2} leaves and {rets2} reticulations."); */
        heuristic::find_heuristic_sequence_multi_label(n1, n2, num_iterations)
    }
    fn add_leaf(&mut self, i: usize, label: T) {
        let mut labels = BTreeSet::new();
        labels.insert(label);
        let new_leaf = Node::new_leaf(i, labels);
        self.add_node(new_leaf);
    }
    fn add_empty_leaf(&mut self, i: usize) {
        let labels = BTreeSet::new();
        let new_leaf = Node::new_leaf(i, labels);
        self.add_node(new_leaf);
    }
    fn add_multi_leaf(&mut self, i: usize, labels: BTreeSet<T>) {
        let new_leaf = Node::new_leaf(i, labels);
        self.add_node(new_leaf);
    }
    //add a tree or reticulation node at index i
    fn add_node_i(&mut self, i: usize) {
        let new_node = Node::new_node(i);
        self.add_node(new_node);
    }
    fn add_node(&mut self, new_node: Node<T>) {
        //add new node to the network nodes array
        //node as arg
        let index = new_node.get_id();
        if index >= self.nodes.len() {
            self.nodes.resize(index + 1, None)
        }
        //nodes being added when building base tree
        self.nodes.replace_w_some(new_node, index);
    }
    fn add_edge(&mut self, outnode_i: usize, innode_i: usize) {
        //add edge to network, updates adjacencies of given nodes
        //TODO doesnt update node types
        //add edge to network, updates adjacencies of given nodes
        //TODO doesnt update node types
        let parent = self.nodes[outnode_i].as_mut().unwrap();
        if parent.is_ret() && !parent.get_children_i().is_empty() {
            panic!("reticulation cannot have more than 1 child.");
        }
        let outnode = self.get_node_mut(outnode_i);
        outnode.add_child(innode_i);
        let innode = self.get_node_mut(innode_i);
        innode.add_parent(outnode_i);
    }
    fn remove_edge(&mut self, parent_i: usize, child_i: usize) {
        //note: error checking done in get_node_mut
        let parent = self.get_node_mut(parent_i);
        parent.remove_adj(child_i);
        let child = self.get_node_mut(child_i);
        child.remove_adj(parent_i);
    }
    fn subdivide_edge(&mut self, (a, b): (usize, usize), node_i: usize) {
        //subdivide given edge with given node
        //subdivide given edge with given node
        self.remove_edge(a, b);
        self.add_edge(a, node_i);
        self.add_edge(node_i, b);
    }
    fn compress_path(&mut self, a: usize, b: usize, c: usize) {
        //compress path abc to ac (assuming b is unary)
        //compress path abc to ac (assuming b is unary)
        self.remove_edge(a, b);
        self.remove_edge(b, c);
        self.add_edge(a, c);
        self.nodes.replace_w_none(b);
    }
    fn add_reticulation(
        &mut self,
        (a, a_child): (usize, usize),
        (b, b_child): (usize, usize),
        bicon_id: usize,
    ) -> usize {
        //adds a reticulation on given outedge and inedge
        //adds a reticulation on given outedge and inedge
        //a is above(out edge), b is below(in edge)
        let yi = self.nodes.len();
        let xi = yi + 1;
        let a_cluster = self.get_cluster(a_child);
        let b_cluster = self.get_cluster(b_child);
        let mut ab_cluster = a_cluster.union(&b_cluster).cloned().collect();
        let b_cluster = self.get_cluster(b_child);
        //add vertex y below a and above achild
        self.add_node_i(yi);
        self.update_cluster(yi, &mut ab_cluster);
        self.subdivide_edge((a, a_child), yi);
        //add reticulation x below b and above bchild
        self.add_node_i(xi);
        self.update_cluster(xi, &mut b_cluster.clone());
        self.subdivide_edge((b, b_child), xi);
        //reticulation edge
        self.add_edge(yi, xi);
        self.add_addn_edge((yi, xi));
        //update clusters
        self.bubble_cluster_up(a, &b_cluster);
        //update the nontriv_bicon list
        let bicon = self.get_bicon(xi);
        for node_i in bicon {
            self.nontriv_bicons.insert(node_i, bicon_id);
        } //establishes they are all in bicon r.
        xi
    }
    fn remove_leaf(&mut self, leaf_i: usize) {
        //used on DTrees, return value is id
        let leaf_p = self.get_parents_i(leaf_i)[0];
        self.remove_edge(leaf_p, leaf_i);
        self.nodes.replace_w_none(leaf_i);
    }
    fn remove_subnet(&mut self, v: usize) {
        //if self.nontriv_bicons.contains_key(&v) {
        //panic!("Cannot remove subnet of vertex in a non-trivial biconnected component.")
        //}//possible this test doesnt work since called on retic? is it?
        //if self.nontriv_bicons.contains_key(&v) {
        //panic!("Cannot remove subnet of vertex in a non-trivial biconnected component.")
        //}//possible this test doesnt work since called on retic? is it?
        let mut lower_nodes = vec![];
        for op_node in &self.nodes {
            if let Some(node) = op_node {
                let n = node.get_id();
                if self.is_below(n, v) {
                    if n != v {
                        lower_nodes.push(n);
                    }
                }
            }
        }
        let u = self.get_parents_i(v)[0];
        for n in lower_nodes {
            self.nodes.replace_w_none(n);
        }
        self.remove_edge(u, v);
        self.nodes.replace_w_none(v);
    }
    // ------ Getters
    fn get_n(&self) -> usize {
        //returns number of nodes in network
        let mut result = 0;
        for n in &self.nodes {
            if n.is_some() {
                result += 1;
            }
        }
        result
    }
    fn get_labels(&self, i: usize) -> BTreeSet<T> {
        self.get_node_ref(i).get_labels()
    }
    fn get_label(&self, i: usize) -> T {
        //assuming there is only 1 label
        self.get_labels(i).first().unwrap().clone()
    }
    fn update_labels(&mut self, target_node: usize, removed_node: usize) {
        // get copy of labels of node being deleted
        let mut removed_node_labels = self.get_node_mut(removed_node).labels.clone();
        self.get_node_mut(target_node)
            .labels
            .append(&mut removed_node_labels);
    }
    fn get_cluster(&self, i: usize) -> BTreeSet<T> {
        self.get_node_ref(i).get_cluster()
    }
    fn get_edges(&self) -> Vec<(usize, usize)> {
        if !self.is_tree() {
            panic!("will contain duplicates if called on network");
        }
        let mut result: Vec<(usize, usize)> = vec![];
        let mut work_list = vec![0];
        while work_list.len() > 0 {
            let cur_i = work_list.pop().unwrap();
            let children = self.get_children_i(cur_i);
            let mut out_edges_i = vec![];
            for child in children {
                work_list.push(child);
                out_edges_i.push((cur_i, child));
            }
            result.append(&mut out_edges_i);
        }
        result
    }
    fn get_ret_nodes(&self) -> Vec<usize> {
        let mut result = vec![];
        for (_, r) in &self.addn_edges {
            result.push(*r);
        }
        result
    }
    fn get_retic_edge_pairs(&self) -> Vec<Vec<(usize, usize)>> {
        let mut result: Vec<Vec<(usize, usize)>> = vec![];
        for (v, r) in &self.addn_edges {
            let parents = self.get_parents_i(*r);
            let a = parents[0];
            let b = parents[1];
            if a != *v {
                result.push(vec![(a, *r), (*v, *r)]);
            } else {
                result.push(vec![(b, *r), (*v, *r)]);
            }
        }
        result
    }
    fn get_rets_below(&self, v: usize) -> BTreeSet<usize> {
        //a ret is not below itself??
        //a ret is not below itself??
        let mut result = BTreeSet::new();
        for r in self.get_ret_nodes() {
            if self.is_below(r, v) && r != v {
                result.insert(r);
            }
        }
        result
    }
    fn get_children_i(&self, i: usize) -> Vec<usize> {
        self.get_node_ref(i).get_children_i()
    }
    fn get_parents_i(&self, i: usize) -> Vec<usize> {
        self.get_node_ref(i).get_parents_i()
    }
    fn get_leaves(&self) -> Vec<usize> {
        let mut results = vec![];
        for i in 0..self.nodes.len() {
            if self.nodes[i].is_some() {
                let n = self.get_node_ref(i);
                let mut is_leaf = true;
                for adj in &n.adjs {
                    if adj.is_to {
                        is_leaf = false;
                    }
                }
                if is_leaf {
                    results.push(i);
                }
            }
        }
        results
    }
    fn get_leaves_count(&self) -> usize {
        self.get_node_ref(0).cluster.len()
    }
    fn get_node_ref(&self, i: usize) -> &Node<T> {
        if self.nodes.len() <= i {
            panic!("cannot get node: i too large");
        }
        match self.nodes[i].as_ref() {
            Some(node) => node,
            None => panic!("no node at i={i}"),
        }
    }
    fn get_node_mut(&mut self, i: usize) -> &mut Node<T> {
        if self.nodes.len() <= i {
            panic!("cannot get node: i too large");
        }
        match self.nodes[i].as_mut() {
            Some(node) => node,
            None => panic!("no node at i={i}"),
        }
    }
    fn get_bicon(&self, ret_node_i: usize) -> Vec<usize> {
        //set up
        let parents = self.get_node_ref(ret_node_i).get_parents_i();
        let mut p1i = parents[0];
        let mut p2i = parents[1];
        if self.is_below(p2i, p1i) {
            //then p1 is bicon root
            //ensures cur is below bicon root
            let temp = p2i;
            p2i = p1i;
            p1i = temp;
        }
        let mut cur = p1i;
        let mut result = vec![ret_node_i];
        let mut up = true;
        'outer: loop {
            if self.get_node_ref(cur).is_leaf() {
                panic!("there can't be a leaf in a nontrivial bicon.");
            }
            result.push(cur);
            if cur == p2i {
                //base case
                break 'outer;
            }
            if self.is_below(p2i, cur) {
                //reached bicon root
                up = false;
            }
            if up {
                cur = self.get_parents_i(cur)[0];
            } else {
                let childs = self.get_children_i(cur);
                let mut bicon_child = 0;
                for (i, child) in childs.iter().enumerate() {
                    if self.is_below(ret_node_i, *child) && !result.contains(child) {
                        bicon_child = i;
                    }
                }
                cur = childs[bicon_child];
            }
        }
        result
    }
    fn get_ret_of_bicon(&self, v: usize) -> usize {
        let mut result = 0;
        if let Some(v_bicon_id) = self.nontriv_bicons.get(&v) {
            for (node_i, bicon_i) in &self.nontriv_bicons {
                if *bicon_i == *v_bicon_id && self.is_ret_node(*node_i) {
                    result = *node_i;
                }
            }
        } else {
            panic!("{v} is not in a nontriv bicon");
        }
        if result == 0 {
            panic!("didn't find the ret in the bicon");
        } else {
            result
        }
    }
    fn get_component_paths(&self, v: usize, r: usize) -> (Vec<usize>, Vec<usize>) {
        //v is bicon root, r is bicon ret
        if !self.is_bicon_root_of(v, r) {
            panic! {"{v} is not bicon root of {r}"}
        }
        //returns left and right component path in vertices (doesnt include v or r) in order
        let mut l_result = vec![];
        let mut r_result = vec![];
        //a node is below itself so the following should be included
        let l_v = self.get_children_i(v)[0];
        let r_v = self.get_children_i(v)[1];
        let bicon = self.get_bicon(r);
        for u in bicon {
            if u != v && u != r {
                if self.is_below(u, l_v) {
                    l_result.push(u);
                } else if self.is_below(u, r_v) {
                    r_result.push(u);
                } else {
                    panic!("member of bicon ({u}) wasn't on a left or right path?");
                }
            }
        }
        //topological sort of both paths!
        let mut sorted_l_result = vec![];
        let mut sorted_r_result = vec![];
        //sorted_l_result.push(l_result.pop());
        while l_result.len() > 0 {
            let mut index = 0;
            let candidate = l_result.pop().unwrap();
            for v in &sorted_l_result {
                if self.is_below(candidate, *v) {
                    index += 1;
                } else {
                    break;
                }
            }
            sorted_l_result.insert(index, candidate);
        }
        while r_result.len() > 0 {
            let mut index = 0;
            let candidate = r_result.pop().unwrap();
            for v in &sorted_r_result {
                if self.is_below(candidate, *v) {
                    index += 1;
                } else {
                    break;
                }
            }
            sorted_r_result.insert(index, candidate);
        }
        (sorted_l_result, sorted_r_result)
    }
    fn get_potential_bicons_leafs(&self) -> Vec<Vec<usize>> {
        //get all bicons for all reticulations
        let (_, rets): (Vec<_>, Vec<_>) = self.addn_edges.iter().cloned().unzip();
        let bicons: Vec<Vec<usize>> = rets.into_iter().map(|ret| self.get_bicon(ret)).collect();
        //(k:bicon_node_i,v:potential_bicon_index)
        let mut group_ids = HashMap::new();
        let mut i = 1; //0th index for nodes above any bicon
        for bicon in bicons {
            for v in bicon {
                group_ids.insert(v, i);
                i += 1;
            }
        }
        let mut potential_bicons: Vec<Vec<usize>> = vec![];
        potential_bicons.resize_with(i, || vec![]);
        let mut processing: Vec<(usize, usize)> = vec![(0, 0)];
        while processing.len() > 0 {
            //works because there are no leaves in non-triv bicons
            let cur = processing.pop().unwrap();
            let cur_node_i = cur.0;
            let cur_bicon_id = cur.1;
            if self.is_leaf_node(cur_node_i) {
                //base case, add leaf to potential bicons
                potential_bicons[cur_bicon_id].push(cur_node_i);
            } else if let Some(new_bicon_id) = group_ids.get(&cur_node_i) {
                //its a bicon node, continue down all children w/ new id
                //cur node is not the outnode of an addn_edge.
                for child in self.get_children_i(cur_node_i) {
                    if !self.is_addn_edge(cur_node_i, child) {
                        processing.push((child, *new_bicon_id));
                    }
                }
            } else {
                //regular tree node, continue down all children
                for child in self.get_children_i(cur_node_i) {
                    processing.push((child, cur_bicon_id));
                }
            }
        }
        potential_bicons.retain(|group| group.len() > 1);
        potential_bicons
    }
    fn get_open_id(&self) -> usize {
        for (i, n) in self.nodes.iter().enumerate() {
            if n.is_none() {
                return i;
            }
        }
        self.nodes.len()
    }
    fn get_open_bicon_id(&self) -> usize {
        let mut values = BTreeSet::new();
        let mut new_id = 0;
        for (_, v) in &self.nontriv_bicons {
            values.insert(v);
        }
        for i in 0.. {
            if !values.contains(&i) {
                new_id = i;
                break;
            }
        }
        new_id
    }
    // ------ Setters
    fn add_nontriv_bicon(&mut self, bicon: Vec<usize>, ret_id: usize) {
        for v in bicon {
            self.nontriv_bicons.insert(v, ret_id);
        }
    }
    fn remove_nontriv_bicon(&mut self, ret: usize) {
        let Some(&ret_id) = self.nontriv_bicons.get(&ret) else {
            panic! {"reticulation is not in a bicon?"}
        };
        self.nontriv_bicons.retain(|_, v| *v != ret_id);
    }
    fn remove_addn_edge(&mut self, ret: usize) {
        self.addn_edges.retain(|(_, b)| *b != ret);
    }
    fn add_addn_edge(&mut self, edge: (usize, usize)) {
        self.addn_edges.push(edge);
    }
    fn update_tree_clusters(&mut self) {
        //propogates labels from leaves up to each node to its cluster var
        //only tested on tree structure, tests removed to be able to update "in-progress" newicks
        let mut new_cluster = self.update_tree_clusters_helper(0);
        self.update_cluster(0, &mut new_cluster)
    }
    fn update_cluster(&mut self, node_i: usize, cluster: &mut BTreeSet<T>) {
        self.get_node_mut(node_i).update_cluster(cluster);
    }
    fn remove_cluster(&mut self, node_i: usize, cluster: &BTreeSet<T>) -> bool {
        //for each item in the set, remove
        let mut result = false; //default false if cluster is empty
        for label in cluster {
            if !self.get_node_mut(node_i).remove_label(label) {
                result = false //there was no such label to remove
            } else {
                result = true //true means keep going
            }
        }
        result
    }
    fn update_tree_clusters_helper(&mut self, cur_node_i: usize) -> BTreeSet<T> {
        let mut new_cluster: BTreeSet<T> = BTreeSet::new();
        if self.is_leaf_node(cur_node_i) {
            //base case: leaf
            for label in self.get_labels(cur_node_i) {
                new_cluster.insert(label.clone());
            }
        } else {
            //recursive case: internal node
            for child_i in &self.get_children_i(cur_node_i) {
                // if !self.is_addn_edge(cur_node_i, *child_i) {
                new_cluster.extend(self.update_tree_clusters_helper(*child_i));
                // }
            }
            self.update_cluster(cur_node_i, &mut new_cluster.clone());
        }
        new_cluster
    }
    fn bubble_cluster_up(&mut self, cur_node_i: usize, new_cluster: &BTreeSet<T>) {
        //starting at index given, keep attaching the given cluster until root.
        let mut new_cluster_clone = new_cluster.clone();
        self.update_cluster(cur_node_i, &mut new_cluster_clone);
        for parent_i in self.get_parents_i(cur_node_i) {
            self.bubble_cluster_up(parent_i, &new_cluster); //don't need to worry about duplicates, its a set
        }
    }
    fn bubble_cluster_up_remove(&mut self, cur_node_i: usize, old_cluster: &BTreeSet<T>) {
        if self.remove_cluster(cur_node_i, old_cluster) {
            for parent in self.get_parents_i(cur_node_i) {
                self.bubble_cluster_up_remove(parent, old_cluster);
            }
        }
    }
    pub fn find_cluster(&self) -> Vec<T> {
        //for DTrees only!!!
        let mut result = vec![];
        for i in 0..self.nodes.len() {
            result.push(self.get_label(i))
        }
        result
    }
    // ------ Functionality
    fn has_no_children(&self, node_i: usize) -> bool {
        let adjs = &self.nodes[node_i].as_ref().unwrap().adjs;
        for adj in adjs {
            if adj.is_to {
                return false;
            }
        }
        true
    }
    fn is_ret_node(&self, node_i: usize) -> bool {
        self.get_node_ref(node_i).is_ret()
    }
    fn is_leaf_node(&self, node_i: usize) -> bool {
        self.get_node_ref(node_i).is_leaf()
    }
    fn is_tree_node(&self, node_i: usize) -> bool {
        self.get_node_ref(node_i).is_tree()
    }
    fn is_tree(&self) -> bool {
        self.addn_edges.is_empty()
    }
    fn is_addn_edge(&self, ai: usize, bi: usize) -> bool {
        //i.e. true if added reticulation edge to base tree
        self.addn_edges.contains(&(ai, bi))
    }
    fn is_trivial(&self, node_i: usize) -> bool {
        !self.nontriv_bicons.contains_key(&node_i)
    }
    fn is_child(&self, potential_child: usize, potential_parent: usize) -> bool {
        self.get_node_ref(potential_child)
            .is_child_of(potential_parent)
    }
    fn is_above(&self, cur_node_i: usize, target_node_i: usize) -> bool {
        //returns true is cur is above target
        //a node is above itself
        (cur_node_i == target_node_i) || self.is_below(target_node_i, cur_node_i)
    }
    fn is_below(&self, cur_node_i: usize, target_node_i: usize) -> bool {
        //returns true of cur is below target
        //a node is below itself
        let cur_cluster = self.get_cluster(cur_node_i);
        let target_cluster = self.get_cluster(target_node_i);
        let is_child = self.is_child(cur_node_i, target_node_i);
        if cur_node_i == target_node_i {
            true
        } else if is_child {
            true
        } else if cur_cluster == target_cluster && !is_child {
            false
        } else if cur_cluster.is_subset(&target_cluster) {
            true
        } else {
            false
        }
    }
    fn is_bicon_root_of(&self, v: usize, r: usize) -> bool {
        //returns if v is the root of the bicon containing ret r
        let mut bicon = self.get_bicon(r);
        bicon.retain(|x| *x != v);
        bicon.iter().all(|u| self.is_above(v, *u))
    }
    fn is_iso(&self, n: &Self, s_i: usize, n_i: usize) -> bool {
        //checks equality through labels
        if self.get_label(s_i) == n.get_label(n_i) {
            let s_chs = self.get_children_i(s_i);
            let n_chs = n.get_children_i(n_i);
            if s_chs.len() == n_chs.len() {
                let mut childs_match = vec![];
                for s_ch in s_chs {
                    let s_ch_label = self.get_label(s_ch);
                    for n_ch in &n_chs {
                        let n_ch_label = n.get_label(*n_ch);
                        let mut result = false;
                        if s_ch_label == n_ch_label {
                            result = self.is_iso(n, s_ch, *n_ch);
                            break;
                        } //assumes all labels are distinct
                        childs_match.push(result);
                    }
                }
                //true when empty
                childs_match.iter().all(|x| *x)
            } else {
                false
            }
        } else {
            false
        }
    }
    // ------------------------ NETWORK TESTS
    /// Primarily for testing purposes, will return true if there are no
    /// (directed) cycles in the network and false otherwise.
    pub fn is_acyclic(&mut self) -> bool {
        //Attempt a topological sort of the network,
        // and see if remaining graph has edges
        //NOTE possibly doesnt test for cycle containing the root??
        //Q ← Set of all nodes with no incoming edges(the root)
        let mut q_list = vec![0 as usize];
        while q_list.len() > 0 {
            let node_i = q_list.pop().unwrap();
            let n_childs = self.get_children_i(node_i);
            for n_child_i in n_childs {
                self.remove_edge(node_i, n_child_i);
                // if m has no other incoming edges then push children to Q
                if self.get_parents_i(n_child_i).len() == 0 {
                    q_list.push(n_child_i);
                }
            }
        }
        if self.has_edges() {
            //graph has a cycle
            false
        } else {
            true
        }
    }
    fn has_edges(&self) -> bool {
        //for testing purposes (acyclicity test)
        let mut result = false;
        for (i, candidate) in self.nodes.iter().enumerate() {
            if let Some(_) = candidate {
                if self.get_node_ref(i).has_adj() {
                    result = true;
                }
            }
        }
        result
    }
    /// Primarily for testing purposes, will return true if the network is
    /// level-1 and false otherwise.
    pub fn is_level_1(&mut self) -> bool {
        //conjecture: the biconnected component(not well defined so "any", i.e. bicon found by the alg defined here) of some 2 different reticulations would be non-disjoint if not level-1
        //cannot rely on info from nontriv-bicons, must go by topology
        //collect rets
        let mut rets: Vec<usize> = vec![];
        for (_, ret) in &self.addn_edges {
            rets.push(*ret);
        }
        let mut bicons = vec![];
        for ret in rets {
            bicons.push(BTreeSet::from_iter(self.get_bicon(ret)));
        }
        bicons.iter().all(|bicon: &BTreeSet<usize>| {
            for other_bicon in &bicons {
                return bicon == other_bicon || bicon.is_disjoint(&other_bicon);
            }
            true
        })
    }
    /// Primarily for testing purposes, will return true of the network is
    /// binary, and false otherwise.
    pub fn is_binary(&mut self) -> bool {
        self.nodes.iter().all(|entry: &Option<Node<T>>| {
            if let Some(node) = entry {
                if (node.is_tree() || node.i == 0) && !node.is_leaf() {
                    if node.i == 0 {
                        if self.get_parents_i(node.i).len() == 0
                            && self.get_children_i(node.i).len() == 2
                        {
                            true
                        } else {
                            false
                        }
                    } else {
                        if self.get_parents_i(node.i).len() == 1
                            && self.get_children_i(node.i).len() == 2
                        {
                            true
                        } else {
                            false
                        }
                    }
                } else if node.is_ret() {
                    if self.get_parents_i(node.i).len() == 2
                        && self.get_children_i(node.i).len() == 1
                    {
                        true
                    } else {
                        false
                    }
                } else if node.is_leaf() {
                    if self.get_parents_i(node.i).len() == 1
                        && self.get_children_i(node.i).len() == 0
                    {
                        true
                    } else {
                        false
                    }
                } else {
                    panic!("this isn't any node type?? {0}", node.i);
                }
            } else {
                true
            }
        })
    }
    // ---- printing the String Network
    fn display_helper(&self, cur_node_i: usize) -> String {
        let mut result = String::new();
        let cur_node = &self.get_node_ref(cur_node_i);
        if !self.is_leaf_node(cur_node_i) {
            result.push('(');
            for child_i in &self.get_children_i(cur_node_i) {
                if !self.is_addn_edge(cur_node_i, *child_i) {
                    result.push_str(&Self::display_helper(&self, *child_i));
                } else {
                    result.push('#');
                    result.push_str(child_i.to_string().as_str());
                }
                result.push(',');
            }
            result.pop();
            result.push(')');
        }
        result.push_str(&cur_node.to_string());
        result
    }
}

impl Network<Node<String>> {
    //TODO pub for testing only!!
    pub fn dep_tree_maker(n: &Network<Node<String>>) -> DTree {
        DTree::dep_tree_maker(n)
    }
    pub fn dep_tree_from_level_seq(seq: Vec<usize>) -> DTree {
        DTree::build_from_level_seq(seq)
    }
    // ----------------------- public function
    /// Generates new random network
    /// # Arguments
    /// * `size` - the number of desired leaves in network
    /// * `r` - the maximum numer of reticulations in the network
    /// (occurance of this amount is not guaranteed)
    pub fn new_random(size: usize, r: usize, exact: bool) -> Self {
        random::new_random(size, r, exact)
    }
    /// Given a network n, randomly modifies n and returns the original
    /// and the modified network. Random modifications are done by cherry
    /// operations the number of which (likely) correspond to the construction,
    /// deconstruction, or tail distance.
    /// There will be 2 leaf labels in common in the outputs
    /// # Arguments
    /// * `net` - a network
    /// * `dist` - the desired distance (in cherry operations) between
    /// the two output networks
    pub fn random_modify(net: Self, dist: usize) -> (Self, Self) {
        cherry::random_modify(net, dist)
    }
    ///
    ///  # Arguments
    ///  * `newick_string` - a string slice of a potential newick string
    pub fn parse_newick(newick_string: &str) -> Self {
        //newick_string.retain(|c| !c.is_whitespace() );
        if newick_string.len() == 0 {
            panic!("cannot parse empty string");
        } else {
            let mut cur_net = Network::new(0);
            newick::network(newick_string, &mut cur_net);
            cur_net
        }
    }
    /// Primarily for testing purposes, will return true of the network is
    /// orchards, and false otherwise.
    pub fn is_orchard(&mut self) -> bool {
        cherry::is_orchard(self)
    }
}
//  ----------------------  Network Display
impl<T: Clone + PartialEq + Eq + Hash + Ord + Display> std::fmt::Display for Network<Node<T>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let result;
        if self.nodes[0].is_none() {
            result = String::from("Network without root");
        } else {
            result = Self::display_helper(&self, 0);
        }
        write!(f, "{}", result)
    }
}
//   ------------------------------- NODE Implementation

//   ------------------------------- NODE
impl<T: Clone + PartialEq + Eq + Hash + Ord + Display> Node<T> {
    //Constructors
    fn new(i: usize, cluster: BTreeSet<T>, labels: BTreeSet<T>, adjs: Vec<Adj>) -> Self {
        Node {
            i: i,
            cluster: cluster,
            labels: labels,
            adjs: adjs,
        }
    }
    fn new_leaf(i: usize, labels: BTreeSet<T>) -> Self {
        Node {
            i: i,
            cluster: labels.clone().into_iter().collect(),
            labels: labels,
            adjs: vec![],
        }
    }
    fn new_node(i: usize) -> Self {
        Node {
            i: i,
            cluster: BTreeSet::new(),
            labels: BTreeSet::new(),
            adjs: vec![],
        }
    }
    //getters
    fn get_labels(&self) -> BTreeSet<T> {
        self.labels.clone()
    }
    fn get_cluster(&self) -> BTreeSet<T> {
        self.cluster.clone()
    }
    fn get_id(&self) -> usize {
        self.i
    }
    fn get_children_i(&self) -> Vec<usize> {
        let children_i: Vec<usize> = self
            .adjs
            .iter()
            .filter(|adj| adj.is_to)
            .map(|adj| adj.i)
            .collect();
        children_i
    }
    fn get_parents_i(&self) -> Vec<usize> {
        let parents_i: Vec<usize> = self
            .adjs
            .iter()
            .filter(|adj| !adj.is_to)
            .map(|adj| adj.i)
            .collect();
        parents_i
    }
    //setters
    fn remove_adj(&mut self, adj_i: usize) {
        let mut index = 0;
        let mut found = false;
        for (i, adj) in self.adjs.iter().enumerate() {
            if adj.i == adj_i {
                found = true;
                index = i;
                break;
            }
        }
        if found {
            self.adjs.remove(index);
        } else {
            panic!("{}'s adj to {adj_i} : no such adjacency to remove", self.i);
        }
    }
    fn update_cluster(&mut self, new_cluster: &mut BTreeSet<T>) {
        self.cluster.append(new_cluster);
    }
    fn remove_label(&mut self, old_label: &T) -> bool {
        self.cluster.remove(old_label)
    }
    fn add_child(&mut self, child_i: usize) {
        self.adjs.push(Adj {
            is_to: true,
            i: child_i,
        });
    }
    fn add_parent(&mut self, parent_i: usize) {
        self.adjs.push(Adj {
            is_to: false,
            i: parent_i,
        });
    }
    //functionality
    fn is_leaf(&self) -> bool {
        self.get_parents_i().len() == 1 && !(self.get_children_i().len() > 0)
    }
    fn is_tree(&self) -> bool {
        self.get_children_i().len() == 2
    }
    fn is_ret(&self) -> bool {
        self.get_parents_i().len() > 1
    }
    fn is_child_of(&self, potential_parent: usize) -> bool {
        let mut result = false;
        for adj in &self.adjs {
            if adj.i == potential_parent && adj.is_to == false {
                result = true;
            }
        }
        result
    }
    fn has_adj(&self) -> bool {
        //for testing purposes (acyclicity test)
        if self.adjs.len() > 0 {
            true
        } else {
            false
        }
    }
}
impl<T: Clone + PartialEq + Eq + Hash + Ord + Display> std::fmt::Display for Node<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut result = String::new();
        //  result.push_str("i:");
        if self.is_ret() {
            result.push('#');
            result.push_str(&self.i.to_string().as_str()); //i, not label
        } else {
            result.push_str(&self.i.to_string().as_str()); //i, not label
        }
        // result.push_str(&self.i.to_string().as_str()); //i, not label
        let binding = self.get_labels();
        if let Some(label) = binding.first() {
            //     result.push_str(" l:");
            result.push_str(label.to_string().as_str());
        }
        write!(f, "{}", result)
    }
}

//  -----------------------------------  REPLACEABLE
trait Replaceable<T: Clone + PartialEq + Eq + Hash + Ord + Display> {
    fn replace_w_some(&mut self, node: Node<T>, index: usize);
    fn replace_w_none(&mut self, index: usize);
}
impl<T: Clone + PartialEq + Eq + Hash + Ord + Display> Replaceable<T> for Vec<Option<Node<T>>> {
    fn replace_w_some(&mut self, node: Node<T>, index: usize) {
        if index >= self.len() {
            panic!("trying to replace with some on index >= length")
        };
        if self[index].is_some() {
            panic!("trying to replace some with some at index {index}");
        }
        self.remove(index);
        self.insert(index, Some(node));
    }
    fn replace_w_none(&mut self, index: usize) {
        if index >= self.len() {
            panic!("trying to replace with none on index >= length")
        };
        if self[index].is_none() {
            panic!("trying to replace none with none at index {index}");
        }
        self.remove(index);
        self.insert(index, None);
    }
}
// ------------------------------------- Clone
impl<T> Clone for Network<T>
where
    T: Clone + PartialEq + Eq + Hash + Display,
{
    fn clone(&self) -> Self {
        let cloned_nodes = self.nodes.clone();
        let cloned_edges = self.addn_edges.clone();

        // Clone the nontriv_bicons HashMap
        let cloned_nontriv_bicons = self
            .nontriv_bicons
            .iter()
            .map(|(&k, &v)| (k, v))
            .collect::<HashMap<usize, usize>>();

        Network {
            nodes: cloned_nodes,
            addn_edges: cloned_edges,
            nontriv_bicons: cloned_nontriv_bicons,
        }
    }
}
// ---------------------------- OTHER
fn gen_label() -> String {
    //generates lowercase alpha strings
    let s: String = rand::thread_rng()
        .sample_iter(Uniform::new(char::from(97), char::from(122)))
        .take(7)
        .map(char::from)
        .collect();
    s
}
