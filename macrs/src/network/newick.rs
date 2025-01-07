//
use crate::network::Network;
use crate::network::Node;
use crate::network::Replaceable;

/// state machine parses given newick string, 
/// and the mutable reference to a network that it updates as 
/// the string is parsed
/// production rules for Extended Newick Format
//Network → Branch ";"
//Subnet → Leaf | Internal
//Leaf → Name Hybrid
//Internal → "(" BranchSet ")" Name
//BranchSet → Branch | Branch "," BranchSet
//Branch → Subnet Length
//Name → empty | string
//Length → empty | ":" number
//Hybrid → empty | "#" Type integer
//Type → empty | string 
//
//example: (A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;   
// ------------------------------ main production methods
//
pub fn network(newick_string: &str, cur_net: &mut Network<Node<String>>) {
    //Tree -> "(" BranchSet ")" Name Length ";"
    //constructs to root node, peels parens and semi-colon, gives to branchset
    let newick_string = newick_string.trim();
    let (bare_newick_string, potential_last_char) = newick_string.split_at(newick_string.len()-1);
    if let Some(last_char) = potential_last_char.chars().next() {
        if last_char != ';' {
            panic!("Not in newick format: no trailing semicolon.")
        }
        trim_length(newick_string);
        let (trimmed, _) = split_name(bare_newick_string);
        let bare_branch_set = peel_parens(trimmed);
        cur_net.add_node_i(0);
        branch_set(bare_branch_set,cur_net, 0);
    } else {
        panic!("cannot parse empty string");//this error is already handled
    }
    //update tree cluster
    cur_net.update_tree_clusters();
    //join reticulations
    join_rets(cur_net);
}
fn subnet(newick_string: &str, cur_net: &mut Network<Node<String>>, parent: usize) {
    //Subnet -> Leaf | Internal 
    if newick_string.contains("(") {
        internal(newick_string, cur_net, parent)
    } else {
        leaf(newick_string, cur_net, parent)
    }
}
fn internal(newick_string: &str, cur_net: &mut Network<Node<String>>, parent: usize) {
    //Internal → "(" BranchSet ")" Name Hybrid
    let (bare_subnet, node_name) = split_name(newick_string);
    let i = cur_net.get_open_id();
    if node_name.contains("#") {
        //add reticulation node 
        let ret_id: usize = get_ret_id(node_name).parse().unwrap();
        cur_net.add_node_i(i);
        //use add_edges to temp store ret info for rejoining later
        cur_net.add_addn_edge((ret_id, i));
    } else {
        //add internal node
        cur_net.add_node_i(i);
    }
    cur_net.add_edge(parent, i);
    //process children
    let branch_set_str = peel_parens(bare_subnet);
    branch_set(branch_set_str, cur_net, i);
}
fn branch_set(newick_string: &str, cur_net: &mut Network<Node<String>>, parent: usize) {
    //BranchSet → Branch | Branch "," BranchSet
    if let Some((branch_str,branch_set_str)) = split_first_branch(newick_string) {
        branch(branch_str, cur_net, parent);
        branch_set(branch_set_str, cur_net, parent);
    } else {
        branch(newick_string, cur_net, parent);
    }
}
fn branch(newick_string: &str, cur_net: &mut Network<Node<String>>, parent: usize) {
    //Branch → Subnet Length
    let bare_subnet = trim_length(newick_string);
    subnet(bare_subnet, cur_net, parent)
}
fn leaf(newick_string: &str, cur_net: &mut Network<Node<String>>, parent: usize) {
    //Leaf → Name Hybrid
    //Name → empty | string
    //Hybrid → empty | "#" Type integer
    let i = cur_net.get_open_id();
    if newick_string.contains("#") {
        cur_net.add_empty_leaf(i);
        let ret_id: usize = get_ret_id(newick_string).parse().unwrap();
        //use add_edges to temp store ret info for rejoining later
        cur_net.add_addn_edge((ret_id, i));
    } else {
        if newick_string.len() == 0 {
            panic!("all leaves must be labelled.");
        }
        let newick_string = newick_string.trim();
        cur_net.add_leaf(i, newick_string.to_owned());
    }
    cur_net.add_edge(parent, i);
}
//
//  -------------------------------    helper functions
//
fn trim_length(s: &str) -> &str {
    //only trim if ':' is found
    let mut trimmable = false;
    let mut split_at = 0;
    for (i,c) in s.chars().rev().enumerate() {
        match c {
            ')' => break,
            ':' => {trimmable = true; split_at = i; break},
            '0'..='9' | '.' => continue,
            _ => break,
        }
    }  
    if trimmable {
        let (return_s,_) = s.split_at(s.len()-split_at-1);
        return_s
    } else {
        s
    }
}
fn split_name(s: &str) -> (&str, &str) {
    if let Some(i) = s.rfind(')') {
        s.split_at(i+1)
    } else {
        panic!("no subnet to split from");
    }
}
fn peel_parens(s: &str) -> &str {
    if let Some(s) = s.strip_prefix("(") {
        if let Some(s) = s.strip_suffix(")") {
            s
        } else {
            panic!("no parens to peel {s}");
        }
    } else {
        panic!("no parens to peel {s}");
    }
}
fn split_first_branch(s: &str) -> Option<(&str, &str)> {
    //go through char by char, 
    //if not in an inner paren, split on comma
    let mut split_index = 0;
    let mut is_splittable = false;
    let mut is_inner = 0;
    for (i,c) in s.chars().enumerate() {
        match c {
            '(' => is_inner +=1,
            ')' => is_inner -=1,
            ',' => if is_inner == 0 {split_index = i;is_splittable = true;break;},
            _ => continue,
        }
    }
    if is_splittable {
        let (f,pre_l) = s.split_at(split_index);
        let l = pre_l.strip_prefix(",").unwrap();
        Some((f,l))
    } else {
        None
    }
}
fn get_ret_id(node_name: &str) -> &str {
    //Hybrid → empty | "#" Type integer
    node_name.trim_start_matches(|c:char| !c.is_numeric() ) 
}
fn join_rets(net: &mut Network<Node<String>>) {
    //first, find an index in ret pairs to associate each ret with
    let mut ret_pairs = vec![]; 
    ret_pairs.resize(net.addn_edges.len()/2, (0,0));
    let mut reindex_rets = vec![];
    for (ret_id,_) in &net.addn_edges {
        if !reindex_rets.contains(ret_id) {
            reindex_rets.push(*ret_id);
        }
    }
    reindex_rets.sort();
    //associate each corresonding leaf/ret
    for (ret_id,v) in &net.addn_edges { 
        //binary search returns index of matched element
        let ret_index = match reindex_rets.binary_search(ret_id) {
            Ok(index) => index,
            Err(error) => panic!("cannot find ret_id when indexing.error index:{:?}", error),
        };
        if net.is_leaf_node(*v) {
            ret_pairs[ret_index].0 = *v;
        } else {
            ret_pairs[ret_index].1 = *v;
        }
    }
    //clear temp data from addn_edges
    net.addn_edges.clear();
    //remove the leaf and join that edge to the ret, updating network struct
    for (leaf,ret) in ret_pairs {
        //delete leaf
        let p = net.get_parents_i(leaf)[0];
        net.remove_edge(p,leaf);
        net.nodes.replace_w_none(leaf);
        //add edge to reticulation proper
        net.add_edge(p,ret);
        //make it an addn_edge
        net.add_addn_edge((p,ret));
        //update cluster
        net.bubble_cluster_up(p, &net.get_cluster(ret));
        //update non-triv bicons
        let bicon_id = net.get_open_bicon_id();
        for node_i in net.get_bicon(ret) {
            net.nontriv_bicons.insert(node_i, bicon_id);
        }
    }
}