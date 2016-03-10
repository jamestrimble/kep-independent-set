import sys
import random

class OptimisationException(Exception):
    pass

class PoolOptimiser(object):
    EPSILON = 1e-7

    def __init__(self, pool, opt_criteria, max_cycle, max_chain):
        self.opt_criteria = opt_criteria
        self.pool = pool
        self.cycles = pool.find_cycles(max_cycle)
        self.chains = pool.find_chains(max_chain)

    def add_clique(self, adj_mat, node_ids):
        for i in range(len(node_ids) - 1):
            for j in range(i+1, len(node_ids)):
                x = node_ids[i]
                y = node_ids[j]
                adj_mat[x][y] = True
                adj_mat[y][x] = True

    def calc_hier_chain_score(self, chain):
        val = 0
        for oc in self.opt_criteria:
            val <<= 10
            if oc.sense == 'MAX':
                val += min(10, int(oc.chain_val(chain)))
            else:
                val += 1023 - int(oc.chain_val(chain))
        return val

    def calc_hier_cycle_score(self, cycle):
        val = 0
        for oc in self.opt_criteria:
            val <<= 10
            if oc.sense == 'MAX':
                val += min(10, int(oc.cycle_val(cycle)))
            else:
                val += 1023 - int(oc.cycle_val(cycle))
        return val

    def calc_hier_ndd_score(self, ndd):
        val = 0
        for oc in self.opt_criteria:
            val <<= 10
            if oc.sense == 'MAX':
                val += min(10, int(oc.altruist_val(ndd)))
            else:
                val += 1023 - int(oc.altruist_val(ndd))
        return val

    def are_adj_mat_rows_almost_equal(self, row1, row2, i, j):
        # Returns True iff row1 and row2 are equal except in positions i and j
        for k in range(len(row1)):
            if k != i and k != j and row1[k] != row2[k]:
                return False
        return True

    def reduce_instance(self, adj_mat, hier_scores, descriptions):
        assert len(adj_mat) == len(hier_scores)

        # For each node that is "almost equal" to another node, this gives the equivalence class.
        # i and j are "almost equal" if they are connected by an edge, and if
        # their neighbourhoods with i and j excluded are equal
        node_to_equiv_class = {}

        for i in range(len(adj_mat)-1):
            # TODO: I don't think we need the j loop if i is in node_to_equiv_class
            for j in range(i+1, len(adj_mat)):
                if adj_mat[i][j] and self.are_adj_mat_rows_almost_equal(adj_mat[i], adj_mat[j], i, j):
                    if i in node_to_equiv_class:
                        equiv_class = node_to_equiv_class[i]
                    else:
                        assert j not in node_to_equiv_class
                        equiv_class = NodeEquivClass()
                        equiv_class.add_node(i, hier_scores[i])
                        node_to_equiv_class[i] = equiv_class
                        node_to_equiv_class[j] = equiv_class
                    equiv_class.add_node(j, hier_scores[j])

        nodes_to_keep = []
        for i in range(len(adj_mat)):
            if i not in node_to_equiv_class or node_to_equiv_class[i].best_node_id == i:
                nodes_to_keep.append(i)

        reduced_hier_scores = [hier_scores[i] for i in nodes_to_keep]
        reduced_adj_mat = [[adj_mat[i][j] for j in nodes_to_keep] for i in nodes_to_keep]
        reduced_descriptions = [descriptions[i] for i in nodes_to_keep]

        return len(nodes_to_keep), reduced_hier_scores, reduced_adj_mat, reduced_descriptions

    def solve(self, invert_edges, reduce_nodes, node_order):
        # param node_order: 0=default, 1=random, 2=score ascending, 3=score descending
        #                   4=degree asc., 5=degree desc.
        
        patients = self.pool.patients
        paired_donors = self.pool.paired_donors
        altruists = self.pool.altruists

#        for c in self.cycles:
#            print c.index, "   ", " ".join(str(pdp.patient.index) for pdp in c.pd_pairs)

        chain_node_ids = range(len(self.chains))
        cycle_node_ids = [len(chain_node_ids) + i for i in range(len(self.cycles))]
        unused_ndd_node_ids = [len(chain_node_ids) + len(cycle_node_ids) + i for i in range(len(altruists))]
        num_nodes = len(chain_node_ids) + len(cycle_node_ids) + len(unused_ndd_node_ids)
        #print chain_node_ids
        #print cycle_node_ids
        #print unused_ndd_node_ids

        # We'll use "node" to denote a vertex in our MWIS instance
        # Nodes have zero-based indices, but we'll print them out using 1-based indexing
        patient_to_nodes = {patient: [] for patient in patients}
        paired_donor_to_node = {donor: [] for donor in paired_donors}
        ndd_to_node = {ndd: [] for ndd in altruists}

        # Each element of hier_scores will be the full hierarchy of scores for a node, compressed
        # into a single int
        hier_scores = [0] * num_nodes

        descriptions = [] # A description of each node, for printing in the comments

        for c, node_id in zip(self.chains, chain_node_ids):
            for pd_pair in c.pd_pairs:
                patient_to_nodes[pd_pair.patient].append(node_id)
                paired_donor_to_node[pd_pair.donor].append(node_id)
            ndd_to_node[c.altruist_edge.altruist].append(node_id)
            hier_scores[node_id] = self.calc_hier_chain_score(c)
            descriptions.append("{} {}".format(
                    c.altruist_edge.altruist.nhs_id,
                    " ".join("({},{})".format(pd_pair.patient.nhs_id, pd_pair.donor.nhs_id)
                                              for pd_pair in c.pd_pairs)))

        for c, node_id in zip(self.cycles, cycle_node_ids):
            for pd_pair in c.pd_pairs:
                patient_to_nodes[pd_pair.patient].append(node_id)
                paired_donor_to_node[pd_pair.donor].append(node_id)
            hier_scores[node_id] = self.calc_hier_cycle_score(c)
            descriptions.append(" ".join("({},{})".format(
                    pd_pair.patient.nhs_id, pd_pair.donor.nhs_id) for pd_pair in c.pd_pairs))

        for ndd, node_id in zip(altruists, unused_ndd_node_ids):
            ndd_to_node[ndd].append(node_id)
            hier_scores[node_id] = self.calc_hier_ndd_score(ndd)
            descriptions.append(str(ndd.nhs_id))

        adj_mat = [[False]*num_nodes for i in range(num_nodes)]
        
        for d in [patient_to_nodes, paired_donor_to_node, ndd_to_node]:
            for key, val in d.iteritems():
                self.add_clique(adj_mat, val)

        if reduce_nodes:
            old_num_nodes = None
            while old_num_nodes is None or old_num_nodes > num_nodes:
                print "c Reducing ...", num_nodes, "nodes"
                old_num_nodes = num_nodes
                num_nodes, hier_scores, adj_mat, descriptions = self.reduce_instance(
                        adj_mat, hier_scores, descriptions)
                
        if node_order in [1, 2, 3, 4, 5]:
            if node_order == 1:
                order = range(num_nodes)
                random.shuffle(order)
            elif node_order == 2:
                order = sorted(range(num_nodes), key=lambda i: hier_scores[i])
            elif node_order == 3:
                order = sorted(range(num_nodes), key=lambda i: hier_scores[i], reverse=True)
            elif node_order == 4:
                order = sorted(range(num_nodes), key=lambda i: sum(adj_mat[i]))
            elif node_order == 5:
                order = sorted(range(num_nodes), key=lambda i: sum(adj_mat[i]), reverse=True)
            hier_scores = [hier_scores[i] for i in order]
            adj_mat = [[adj_mat[i][j] for j in order] for i in order]
            descriptions = [descriptions[i] for i in order]

        if num_nodes < 30:
            print "c Adjacency matrix"
            for row in adj_mat:
                print "c " + " ".join("X" if x else "." for x in row)

        for i, desc in enumerate(descriptions):
            print "c {} {}".format(i+1, desc)

        if invert_edges:
            print "p edge {} {}".format(num_nodes,
                (num_nodes*num_nodes - num_nodes) / 2 - sum(sum(row) for row in adj_mat) / 2)
        else:
            print "p edge {} {}".format(num_nodes, sum(sum(row) for row in adj_mat) / 2)

        for i in range(num_nodes-1):
            for j in range(i+1, num_nodes):
                edge_exists = adj_mat[i][j]
                if invert_edges:
                    edge_exists = not edge_exists
                if edge_exists:
                    print "e", i+1, j+1

        for i, score in enumerate(hier_scores):
            print "n", i+1, score

class NodeEquivClass(object):
    def __init__(self):
        self.best_score = -1
        self.best_node_id = None
        self.node_ids = set()

    def contains_node(self, node_id):
        return node_id in self.node_ids

    def add_node(self, node_id, score):
        self.node_ids.add(node_id)
        if score > self.best_score:
            self.best_node_id = node_id
            self.best_score = score

    
