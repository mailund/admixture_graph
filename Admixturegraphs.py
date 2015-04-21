# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Admixture graph fitting

# <headingcell level=2>

# Graph representation

# <markdowncell>

# Representation of admixture graph nodes. If there are more than one parent it is an admixture node, otherwise it is just a merge of populations.

# <codecell>

class Node(object):
    def __init__(self, name = ''):
        self.name = name
        self.parents = []
        self.children = []
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return 'Node({})'.format(self.name)
        
    def connect_child(self, child):
        self.children.append(child)
        child.parents.append(self)
        
    def connect_parent(self, parent):
        self.parents.append(parent)
        parent.children.append(self)

# <markdowncell>

# Just a few functions for debugging, making it easier to read paths and sets.

# <codecell>

def pretty_edge(edge):
    src, dst = edge
    return '({},{})'.format(src.name, dst.name)

def pretty_path(path):
    return '[{}]'.format(' '.join(pretty_edge(edge) for edge in path))

def pretty_paths(paths):
    return '; '.join(pretty_path(path) for path in paths)

def pretty_cost(s):
    return '{{{}}}'.format(','.join('({},{})'.format(src,dst) for src,dst in s))

# <markdowncell>

# We need a mapping of all edges in the graph in order to build the linear system for estimating branch lengths.

# <codecell>

def enumerate_edges(root_leaf):
    seen = set()
    edges = set()
    
    def dfs(node):
        if node in seen:
            return
        seen.add(node)
        
        for parent in node.parents:
            if (node,parent) not in edges and (parent,node) not in edges:
                edges.add((node,parent))
            dfs(parent)
            
        for child in node.children:
            if (child,node) not in edges and (node,child) not in edges:
                edges.add((child,node))
            dfs(child)
            
    dfs(root_leaf)
    
    numbered_edges = {}
    for i, (s, d) in enumerate(edges):
        numbered_edges[(s,d)] = i
        numbered_edges[(d,s)] = i
        
    return numbered_edges

# <headingcell level=2>

# Computing expected f statistics from an admixture graph

# <markdowncell>

# Given a source and a destination we compute all the paths between the two nodes.
# 
# The f statistics are based on the overlap between paths between various configurations of nodes, so we also need a way to get the edges found in both of two paths; here the direction of the edge matters since we need to add a positive contribution when the direction is the same in the two paths and a negative contribution when the paths move in opposite direction.
# 
# When we compute the f statistics and there are several paths between two nodes, these paths are weighted with the admixture proportions along each path. These we need to extract for the individual paths. The weight of a pair of paths is the product of the weight of each.

# <codecell>

def generate_paths(src, dst):
    paths = []
    
    def U(s,e):
        x = s.copy()
        x.add(e)
        return x
    
    def dfs(node, seen, current_path, going_up = True):
        if node == dst:
            paths.append(current_path)
        else:
            if going_up:
                for parent in node.parents:
                    if parent not in seen:
                        dfs(parent, U(seen, parent), current_path + [(node, parent)], True)
            for child in node.children:
                if child not in seen:
                    dfs(child, U(seen, child), current_path + [(node, child)], False)

    dfs(src, {src}, [])        
    return paths

def get_overlapping_edges(path1, path2):
    edges1 = set(path1)
    edges2 = set(path2)
    overlap = edges1.intersection(edges2)
    
    reverse_edges2 = set((dst, src) for src, dst in path2)
    reverse_overlap = edges1.intersection(reverse_edges2)
    
    return overlap, reverse_overlap

def is_admixture_edge(edge):
    src, dst = edge
    if dst in src.parents and len(src.parents) > 1:
        return True
    if dst in src.children and len(dst.parents) > 1:
        return True
    return False

def path_admixture_probability(path):
    admix_edges = set()
    for edge in path:
        if is_admixture_edge(edge):
            admix_edges.add(edge)
    return admix_edges

# <markdowncell>

# The f statistics can now be computed. The f2 and f3 statistics can be phrased in terms of the f4 statistic so that is the only one that requires some work. It is computed by getting all the paths between the relevant pairs of leaves, together with their admixture weight. Then take the overlap of edges between each pair of paths, distinguishing between edges where the paths go in the same direction or the opposite direction. For each edge in overlapping paths we add or subtract -- depending on directions in the paths -- the admixture weight for the pair of paths for the coefficient of that edge. The result is a polynomial/vector describing the statistcs as a polynomial over edges. The inner product of this vector with a vector of branch lengths would give the relevant statistics.

# <codecell>

def f4_edges(a, b, c, d):
    '''f4(a,b; c,d) -- overlap of paths from a to b with paths from c to d.'''
    ab_paths = generate_paths(a, b)
    cd_paths = generate_paths(c, d)
    
    def U(s1, s2):
        return s1.union(s2)
    
    weighted_edges = []
    for ab in ab_paths:
        ab_probability = path_admixture_probability(ab)
        for cd in cd_paths:
            cd_probability = path_admixture_probability(cd)
            positive, negative = get_overlapping_edges(ab, cd)
            weighted_edges.append( (U(ab_probability, cd_probability), positive, negative) )
    
    return weighted_edges

def f3_edges(a, b, c):
    return f4_edges(a, b, a, c)

def f2_edges(a, b):
    return f4_edges(a, b, a, b)

def edges_to_polynomial(weighted_edges, edge_mapping, admixture_proportions):
    from numpy import prod, zeros
    
    def calc_cost(cost):
        if len(cost) == 0:
            return 1.0
        else:
            return prod([admixture_proportions[edge] for edge in cost])
    
    polynomia = zeros(len(edge_mapping) / 2)
    for cost, pos, neg in weighted_edges:
        cost = calc_cost(cost)
        for edge in pos:
            polynomia[edge_mapping[edge]] += cost
        for edge in neg:
            polynomia[edge_mapping[edge]] -= cost
    return polynomia

def f4(a, b, c, d, edge_mapping, admixture_proportions):
    return edges_to_polynomial(f4_edges(a, b, c, d), edge_mapping, admixture_proportions)

def f3(a, b, c, edge_mapping, admixture_proportions):
    return f4(a, b, a, c, edge_mapping, admixture_proportions)

def f2(a, b, edge_mapping, admixture_proportions):
    return f4(a, b, a, b, edge_mapping, admixture_proportions)

# <headingcell level=2>

# Constructing the equation system

# <markdowncell>

# The set of f2 and f3 statistics for all pairs and triplets starting with a given leaf node defines the space of equations we need to consider.

# <codecell>

def get_all_pairs(leaf_root, leaves):
    pairs = [(leaf_root, leaf) for leaf in leaves if leaf is not leaf_root]
    return pairs

def get_all_triplets(leaf_root, leaves):
    triplets = []
    for j in range(len(leaves)):
        for i in range(j):
            li = leaves[i]
            lj = leaves[j]
            if leaf_root is not li and leaf_root is not lj:
                triplets.append((leaf_root, li, lj))
    return triplets

def enumerate_equation_system(leaf_root, leaves):
    pairs = get_all_pairs(leaf_root, leaves)
    triplets = get_all_triplets(leaf_root, leaves)
    
    index_to_nodes = pairs + triplets
    nodes_to_index = dict((x,idx) for idx, x in enumerate(index_to_nodes))
    
    return index_to_nodes, nodes_to_index

def build_model_matrix(index_to_nodes, edge_mapping, admixture_costs):
    nrows, ncols = len(index_to_nodes), len(edge_mapping)/2
    model_matrix = zeros((nrows,ncols))
    for row, nodes in enumerate(index_to_nodes):
        if len(nodes) == 2:
            model_matrix[row,:] = f2(nodes[0], nodes[1], edge_mapping, admixture_costs)
        elif len(nodes) == 3:
            model_matrix[row,:] = f3(nodes[0], nodes[1], nodes[2], edge_mapping, admixture_costs)
        else:
            assert False, "We should not get here."
    return model_matrix

# <codecell>

def pretty_print_equation_system(index_to_nodes):
    for i, nodes in enumerate(index_to_nodes):
        if len(nodes) == 2:
            print i, ':', 'f2({},{})'.format(nodes[0].name, nodes[1].name)
        elif len(nodes) == 3:
            print i, ':', 'f3({};{},{})'.format(nodes[0].name, nodes[1].name, nodes[2].name)
        else:
            assert False, "We should not get here."

# <headingcell level=2>

# Testing...

# <codecell>

#     R
#    / \
#   S   T
#  / \ /\
# A   B  C
R = Node('R')
S = Node('S'); R.connect_child(S)
T = Node('T'); R.connect_child(T)
A = Node('A'); S.connect_child(A)
B = Node('B'); S.connect_child(B); T.connect_child(B)
C = Node('C'); T.connect_child(C)


edges = enumerate_edges(A)

AB = generate_paths(A,B); BA = generate_paths(B,A)
AC = generate_paths(A,C); CA = generate_paths(C,A)
BC = generate_paths(B,C); CB = generate_paths(B,C)

for edge, number in edges.items():
    print pretty_edge(edge), '->', number

admixture_costs = {}
admixture_costs[(B,S)] = admixture_costs[(S,B)] = 0.1
admixture_costs[(B,T)] = admixture_costs[(T,B)] = 0.9

print pretty_paths(AB)
print pretty_paths(AC)
print f2(A, B, edges, admixture_costs)
print f3(A, B, C, edges, admixture_costs)
print

print pretty_paths(BA)
print pretty_paths(BC)
print f3(B, A, C, edges, admixture_costs)
print


# <codecell>

index_to_nodes, nodes_to_index = enumerate_equation_system(A, [A,B,C])
pretty_print_equation_system(index_to_nodes)
model_matrix = build_model_matrix(index_to_nodes, edges, admixture_costs)
print model_matrix

# <codecell>

from numpy.linalg import lstsq
test_stats = [2.1,2.0,-0.1]
print lstsq(model_matrix, test_stats)

# <codecell>

help(lstsq)

# <codecell>


