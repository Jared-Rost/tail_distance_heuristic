//! Cherry Submodule (of Network)
//!
//! Cherries are sibling leaves (simple) and leaf parents seperated by
//! only a reticulation node (retic)
//!
//! `cherry` submodule contains all functionality regarding cherry operations
//! such as cherry expansions and reductions on a network, as well as the
//! random modification of a network as per the construction, reconstruction,
//! or tail distance but not the mixed distance.
//!
//!ca
//!ca
//! This module contains the type: `Cherry`
//!
use rand::{distributions::Uniform, Rng};
use std::cmp;

use crate::network::gen_label;
use crate::network::Network;
use crate::network::Node;
use crate::network::Replaceable;
use std::time::{Duration, Instant};

// NOTE: at least a cherry must remain in MACRS as there are no special cases of reduce/expand implemented

#[derive(Debug, Clone, PartialEq, Hash, PartialOrd, Eq, Ord)]
pub struct Cherry {
    pub cherry: (usize, usize),
    pub cherry_type: CherryType,
}
#[derive(Debug, PartialEq, Clone, PartialOrd, Eq, Hash, Ord, Copy)]
pub enum CherryType {
    Retic,
    Simple,
}
// ------ Cherry Operations
pub fn get_reducible_cherries(net: &Network<Node<String>>) -> Vec<Cherry> {
    get_reducible_cherries_helper(net, 0, vec![])
}
fn get_reducible_cherries_helper(
    net: &Network<Node<String>>,
    cur_node: usize,
    mut cur_result: Vec<Cherry>,
) -> Vec<Cherry> {
    if net.is_ret_node(cur_node) {
        let child = net.get_children_i(cur_node)[0];
        if !net.is_leaf_node(child) {
            cur_result.append(&mut get_reducible_cherries_helper(net, child, vec![]));
        }
    } else if net.is_leaf_node(cur_node) {
        //do nothing
    } else {
        let children = net.get_children_i(cur_node);
        let c1 = children[0];
        let c2 = children[1];
        let c1_isleaf = net.is_leaf_node(c1);
        let c2_isleaf = net.is_leaf_node(c2);
        let c1_isret = net.is_ret_node(c1);
        let c2_isret = net.is_ret_node(c2);
        if c1_isleaf {
            if c2_isleaf {
                //simple cherry found
                cur_result.push(Cherry {
                    cherry: (c1, c2),
                    cherry_type: CherryType::Simple,
                });
                cur_result.push(Cherry {
                    cherry: (c2, c1),
                    cherry_type: CherryType::Simple,
                });
            }
            if c2_isret {
                let c2_child = net.get_children_i(c2)[0];
                if net.is_leaf_node(c2_child) {
                    //reticulated cherry found
                    cur_result.push(Cherry {
                        cherry: (c2_child, c1),
                        cherry_type: CherryType::Retic,
                    });
                }
            }
        } else if c1_isret {
            if c2_isleaf {
                let c1_child = net.get_children_i(c1)[0];
                if net.is_leaf_node(c1_child) {
                    //reticulated cherry found
                    cur_result.push(Cherry {
                        cherry: (c1_child, c2),
                        cherry_type: CherryType::Retic,
                    });
                }
            }
        }
        if !net.is_addn_edge(cur_node, c1) {
            cur_result.append(&mut get_reducible_cherries_helper(net, c1, vec![]));
        }
        if !net.is_addn_edge(cur_node, c2) {
            cur_result.append(&mut get_reducible_cherries_helper(net, c2, vec![]));
        }
    }
    cur_result
}
pub fn reduce_cherry(net: &mut Network<Node<String>>, cherry: &Cherry) {
    if net.get_n() <= 3 {
        panic!("cannot reduce to less than a cherry network.");
    }
    match cherry.cherry_type {
        CherryType::Simple => reduce_simple_cherry(net, cherry),
        CherryType::Retic => reduce_retic_cherry(net, cherry),
    }
}
/**
 * reduce_cherry_multi_label
 *
 * Given a network and a cherry to reduce that network by, reduce the network and ensure that if the cherry is simple then all labels are preserved in the leaf that was not reduced
 * @param &Network<Node<String>> net - input network
 * @param &Cherry cherry - cherry to reduce
 */
pub fn reduce_cherry_multi_label(net: &mut Network<Node<String>>, cherry: &Cherry) {
    if net.get_n() <= 3 {
        panic!("cannot reduce to less than a cherry network.");
    }
    match cherry.cherry_type {
        CherryType::Simple => reduce_simple_cherry_multi_label(net, cherry),
        CherryType::Retic => reduce_retic_cherry(net, cherry),
    }
}
fn reduce_simple_cherry(net: &mut Network<Node<String>>, cherry: &Cherry) {
    if cherry.cherry_type == CherryType::Retic {
        panic!("can't simple reduce a reticulated cherry.");
    }
    let (x, y) = cherry.cherry;
    let parent = net.get_parents_i(x)[0];
    let grand_parent = net.get_parents_i(parent)[0];
    let x_cluster = net.get_cluster(x);
    net.nodes.replace_w_none(x);
    net.compress_path(grand_parent, parent, y);
    net.bubble_cluster_up_remove(grand_parent, &x_cluster);
}
/**
 * reduce_simple_cherry_multi_label
 *
 * Reduce the network by the provided simple cherry and ensure that labels of the removed leaf is added to the list of labels associated with the preserved leaf
 * @param &Network<Node<String>> net - input network
 * @param &Cherry cherry - cherry to reduce
 */
fn reduce_simple_cherry_multi_label(net: &mut Network<Node<String>>, cherry: &Cherry) {
    if cherry.cherry_type == CherryType::Retic {
        panic!("can't simple reduce a reticulated cherry.");
    }
    let (x, y) = cherry.cherry;
    net.update_labels(y, x); // move all of x labels over to y
    let parent = net.get_parents_i(x)[0];
    let grand_parent = net.get_parents_i(parent)[0];
    let x_cluster = net.get_cluster(x);
    net.nodes.replace_w_none(x);
    net.compress_path(grand_parent, parent, y);
    net.bubble_cluster_up_remove(grand_parent, &x_cluster);
}
fn reduce_retic_cherry(net: &mut Network<Node<String>>, cherry: &Cherry) {
    if cherry.cherry_type == CherryType::Simple {
        panic!("can't reticulated reduce a simple cherry.");
    }
    let (x, y) = cherry.cherry;
    let y_parent = net.get_parents_i(y)[0];
    let y_grand_parent = net.get_parents_i(y_parent)[0];
    let x_parent = net.get_parents_i(x)[0];
    let x_grand_parents = net.get_parents_i(x_parent);
    let x_grand_parent;
    if x_grand_parents[0] == y_parent {
        x_grand_parent = x_grand_parents[1];
    } else {
        x_grand_parent = x_grand_parents[0];
    }
    //compress path will impliclty remove retic edge
    let mut x_cluster = net.get_cluster(x);
    net.compress_path(y_grand_parent, y_parent, y);
    net.compress_path(x_grand_parent, x_parent, x);
    net.remove_addn_edge(x_parent);
    net.remove_nontriv_bicon(x_parent);
    net.bubble_cluster_up_remove(y_grand_parent, &x_cluster);
    net.bubble_cluster_up(x_grand_parent, &mut x_cluster);
}
fn expand_simple_cherry(net: &mut Network<Node<String>>, (xi, yi): (usize, usize)) {
    //error check: (x,y) network must contain x and not y
    let x_none = xi >= net.nodes.len() || net.nodes[xi].is_none();
    let y_present = net.is_leaf_node(yi);
    if !(x_none && y_present) {
        panic!("to expand (x,y), we have assumed y is leaf and present and x is not");
    }
    //if x_present { (x,y) = (cherry.cherry.1, cherry.cherry.0);}
    let parent_y = net.get_parents_i(yi)[0];
    net.remove_edge(parent_y, yi);
    //construct new nodes
    net.add_leaf(xi, gen_label());
    let new_parent = net.get_open_id();
    net.add_node_i(new_parent);
    //connect edges
    net.add_edge(new_parent, xi);
    net.add_edge(new_parent, yi);
    net.add_edge(parent_y, new_parent);
    net.bubble_cluster_up(new_parent, &net.get_cluster(xi));
    let mut y_cluster_copy = net.get_cluster(yi);
    net.get_node_mut(new_parent)
        .update_cluster(&mut y_cluster_copy);
}
fn expand_retic_cherry(net: &mut Network<Node<String>>, (xi, yi): (usize, usize)) {
    let x_present = net.nodes[xi].is_some() && net.is_leaf_node(xi);
    let y_present = net.nodes[yi].is_some() && net.is_leaf_node(yi);
    if !(x_present && y_present) {
        panic!("to expand reticulated cherry we have assumed both x and y present and leaves");
    }
    let xp = net.get_parents_i(xi)[0];
    let yp = net.get_parents_i(yi)[0];
    net.add_reticulation((yp, yi), (xp, xi), net.get_open_bicon_id());
    //delete edges
    //net.remove_edge(yp,yi);
    //net.remove_edge(xp,xi);
    //add new nodes
    //let r = net.get_open_id();
    //net.add_ret(r);
    //let new_yp = net.get_open_id();
    //net.add_tree(new_yp);
    //add edges
    //net.add_edge(xp,r);
    //net.add_edge(r,xi);
    //net.add_edge(new_yp,r);
    //net.add_edge(yp,new_yp);
    //net.add_edge(new_yp,yi);
    //update network
    //net.add_addn_edge((new_yp,r));
    //update clusters
    //let mut cluster = net.get_cluster(xi);
    //net.get_node_mut(r).update_cluster(&mut cluster);
    //cluster = net.get_cluster(xi);
    //net.bubble_cluster_up(new_yp, &mut cluster);
    //update nontriv bicon
    //let ret_id = net.get_open_ret_id();
    //let new_bicon = net.get_bicon(r);
    //net.add_nontriv_bicon(new_bicon, ret_id);
}
pub fn update_reducible_cherries(
    net: &Network<Node<String>>,
    cherry: &Cherry,
    is_reduction: bool,
    cur_reducible_cherries: Vec<Cherry>,
) -> Vec<Cherry> {
    //match on this bool to allow for calling on expansion in the future
    match is_reduction {
        true => update_rc_reduction(net, cherry, cur_reducible_cherries),
        false => cur_reducible_cherries, //update_rc_expansion(cherry, cur_reducible_cherries),
    }
}
fn update_rc_reduction(
    net: &Network<Node<String>>,
    cherry: &Cherry,
    mut cur_reducible_cherries: Vec<Cherry>,
) -> Vec<Cherry> {
    let (x, y) = cherry.cherry;
    if cherry.cherry_type == CherryType::Simple {
        cur_reducible_cherries.retain(|c| c.cherry != (y, x));
        cur_reducible_cherries.retain(|c| c.cherry != (x, y));
        let parent = net.get_parents_i(y)[0];
        let parent_childs = net.get_children_i(parent);
        if net.is_ret_node(parent) {
            for g_parent in net.get_parents_i(parent) {
                for child in net.get_children_i(g_parent) {
                    if child != parent {
                        if net.is_leaf_node(child) {
                            cur_reducible_cherries.push(Cherry {
                                cherry: (y, child),
                                cherry_type: CherryType::Retic,
                            });
                        }
                    }
                }
            }
        } else if net.is_ret_node(parent_childs[0]) {
            //finding if a sibling is ret, y itself cannot be
            let new_x = net.get_children_i(parent_childs[0])[0];
            if net.is_leaf_node(new_x) {
                cur_reducible_cherries.push(Cherry {
                    cherry: (new_x, y),
                    cherry_type: CherryType::Retic,
                });
            }
        } else if net.is_ret_node(parent_childs[1]) {
            //finding if a sibling is ret, y itself cannot be
            let new_x = net.get_children_i(parent_childs[1])[0];
            if net.is_leaf_node(new_x) {
                cur_reducible_cherries.push(Cherry {
                    cherry: (new_x, y),
                    cherry_type: CherryType::Retic,
                });
            }
        } else {
            for child in net.get_children_i(parent) {
                if child != y && net.is_leaf_node(child) {
                    cur_reducible_cherries.push(Cherry {
                        cherry: (child, y),
                        cherry_type: CherryType::Simple,
                    });
                    cur_reducible_cherries.push(Cherry {
                        cherry: (y, child),
                        cherry_type: CherryType::Simple,
                    });
                }
            }
        }
    } else {
        //was reticulate reduction
        cur_reducible_cherries.retain(|c| c.cherry.0 != x);
        let mut parents = net.get_parents_i(x);
        parents.append(&mut net.get_parents_i(y));
        if parents[0] == parents[1] {
            cur_reducible_cherries.push(Cherry {
                cherry: (x, y),
                cherry_type: CherryType::Simple,
            });
            cur_reducible_cherries.push(Cherry {
                cherry: (y, x),
                cherry_type: CherryType::Simple,
            });
        } else {
            for parent in parents {
                let childs = net.get_children_i(parent);
                if net.is_leaf_node(childs[0]) && net.is_leaf_node(childs[1]) {
                    cur_reducible_cherries.push(Cherry {
                        cherry: (childs[0], childs[1]),
                        cherry_type: CherryType::Simple,
                    });
                    cur_reducible_cherries.push(Cherry {
                        cherry: (childs[1], childs[0]),
                        cherry_type: CherryType::Simple,
                    });
                }
            }
        }
    }
    cur_reducible_cherries
}
//fn update_rc_expansion(
//          cherry: &Cherry,
//        mut cur_reducible_cherries: Vec<Cherry>
//  ) -> Vec<Cherry> {
//"generation of orchard and tree-child networks, Cardona et al"
//cur_reducible_cherries
//}
// ---------------------------   testing
pub fn is_orchard(net: &mut Network<Node<String>>) -> bool {
    //if there are always cherries to reduce until network is a cherry
    let mut rc = get_reducible_cherries(&net);
    while net.get_n() > 3 {
        if rc.len() == 0 {
            return false;
        }
        let cherry = rc.swap_remove(0);
        reduce_cherry(net, &cherry);
        rc = update_reducible_cherries(net, &cherry, true, rc);
    }
    if rc.len() > 2 {
        //there will be 2 simple cherries left: (a,b), (b,a)
        false
    } else {
        true
    }
}

// ---------------------------   RANDOM MODIFY
pub fn random_modify(
    net: Network<Node<String>>,
    dist: usize,
) -> (Network<Node<String>>, Network<Node<String>>) {
    //reduces input network by half distance or to a cherry (min of)
    //then expands the rest of the distance
    //the number of reticulations reduced is random according to
    //proportion of reducible cherries at any step are reticulated
    //the number of reticulations expanded is (up to) as many are reduced + 1
    let mut rng = rand::thread_rng();
    let rand_dist = Uniform::from(0usize..99usize); //100%len() not valid index
                                                    //NOTE:
    let mut new_net = net.clone();
    let decon = cmp::min((net.get_n() - 1) / 2 - 1, dist / 2);
    let recon = dist - decon;
    //calculate all reducible pairs
    let mut r_c = get_reducible_cherries(&new_net);
    let mut simple_count = 0;
    for _ in 0..decon {
        let cherry_i = r_c.len() * rng.sample(rand_dist) / 100;
        let cherry = r_c.swap_remove(cherry_i);
        reduce_cherry(&mut new_net, &cherry);
        r_c = update_reducible_cherries(&new_net, &cherry, true, r_c);
        if cherry.cherry_type == CherryType::Simple {
            simple_count += 1;
        }
    }
    //randomly choose positions for reticulation expansions
    let mut retics = recon - simple_count;
    let mut ret_pos = vec![];
    while retics > 0 {
        ret_pos.push(rng.sample(rand_dist) * recon / 100);
        retics -= 1;
    }
    ret_pos.sort();
    ret_pos.dedup();
    //reconstruct
    let mut i = 0; //i only incremented when a expand by a cherry
    while i < recon {
        if ret_pos.contains(&i) {
            //expand retic
            if let Some((v1, v2)) = choose_expandable_cherry(&new_net) {
                expand_retic_cherry(&mut new_net, (v1, v2));
                i += 1;
            } else {
                ret_pos.retain(|x| *x != i);
                ret_pos.push(i + 1); //add a simple cherry and try again
                                     //do not increment i
            }
        } else {
            //expand simple
            let mut yi = 0;
            while new_net.nodes[yi].is_none() || !new_net.is_leaf_node(yi) {
                yi = rng.sample(rand_dist) * new_net.nodes.len() / 100;
            }
            let xi;
            xi = new_net.get_open_id();
            expand_simple_cherry(&mut new_net, (xi, yi));
            i += 1;
        }
    }
    (net, new_net)
}
fn choose_expandable_cherry(net: &Network<Node<String>>) -> Option<(usize, usize)> {
    let mut rng = rand::thread_rng();
    let rand_dist = Uniform::from(0usize..99usize);
    let mut potential_bicons = net.get_potential_bicons_leafs();
    let n: usize = potential_bicons.iter().map(|l| l.len()).sum();
    if potential_bicons.len() > 0 {
        let v1: usize;
        let v2: usize;
        let mut cur_group;
        let mut g = 0; //g is index of cur group
        let mut v_meta_i = rng.sample(rand_dist) * n / 100;
        'v_1: loop {
            for group in &mut potential_bicons {
                for (i, _) in group.iter().enumerate() {
                    if v_meta_i == 0 {
                        v1 = group.swap_remove(i);
                        break 'v_1;
                    }
                    v_meta_i -= 1;
                }
                g += 1;
            }
        }
        cur_group = potential_bicons.swap_remove(g);
        let v_i = rng.sample(rand_dist) * cur_group.len() / 100;
        v2 = cur_group.swap_remove(v_i);
        Some((v1, v2))
    } else {
        None
    }
}
