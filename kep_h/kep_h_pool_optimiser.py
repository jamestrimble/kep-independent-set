from gurobipy import *

class OptimisationException(Exception):
    pass

class PoolOptimiser(object):
    EPSILON = 1e-7

    def __init__(self, pool, opt_criteria, max_cycle, max_chain):
        self.opt_criteria = opt_criteria
        self.pool = pool
        self.cycles = pool.find_cycles(max_cycle)
        self.chains = pool.find_chains(max_chain)

    def _add_var_to_patients_and_donors(self, pd_pairs, mip_var):
        for pd_pair in pd_pairs:
            pd_pair.patient.mip_vars.append(mip_var)
            if len(pd_pair.donor.paired_patients) > 1:
                # If donor has more than one patient, it will need a constraint
                pd_pair.donor.mip_vars.append(mip_var)

    def _optimise(self, model, opt_criterion):
        z = quicksum(
            [chain.mip_var*opt_criterion.chain_val(chain) for chain in self.chains] +
            [cycle.mip_var*opt_criterion.cycle_val(cycle) for cycle in self.cycles] +
            [altruist.mip_var*opt_criterion.altruist_val(altruist)
                        for altruist in self.pool.altruists])
            
        model.setObjective(z)
        model.modelSense = GRB.MAXIMIZE if opt_criterion.sense=='MAX' else GRB.MINIMIZE
        model.optimize()

        return z, model.status

    def _create_vars_and_constraints(self, model, patients, paired_donors, altruists):
        for person in patients + paired_donors + altruists:
            person.mip_vars = []
        
        for chain in self.chains:
            chain.mip_var = model.addVar(vtype=GRB.BINARY, name='chain' + str(chain.index))
        for cycle in self.cycles:
            cycle.mip_var = model.addVar(vtype=GRB.BINARY, name='cycle' + str(cycle.index))
        for altruist in altruists:
            altruist.mip_var = model.addVar(
                    vtype=GRB.BINARY, name='altruist' + str(altruist.nhs_id))
            altruist.mip_vars.append(altruist.mip_var)

        model.update()

        for chain in self.chains:
            chain.altruist_edge.altruist.mip_vars.append(chain.mip_var)
            self._add_var_to_patients_and_donors(chain.pd_pairs, chain.mip_var)
        for cycle in self.cycles:
            self._add_var_to_patients_and_donors(cycle.pd_pairs, cycle.mip_var)

        for person in patients + paired_donors:
            model.addConstr(quicksum(person.mip_vars) <= 1)

        for altruist in altruists:
            model.addConstr(quicksum(altruist.mip_vars) == 1)

    def _enforce_objective(self, model, z, obj_val, sense):
        if sense=='MAX':
            model.addConstr(z >= obj_val)
        else:
            model.addConstr(z <= obj_val)

    def _items_in_optimal_solution(self, items):
        return [item for item in items if item.mip_var.X > 0.5]

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
                val += min(1023, int(oc.chain_val(chain)))
            else:
                val += 1023 - int(oc.chain_val(chain))
        return val

    def calc_hier_cycle_score(self, cycle):
        val = 0
        for oc in self.opt_criteria:
            val <<= 10
            if oc.sense == 'MAX':
                val += min(1023, int(oc.cycle_val(cycle)))
            else:
                val += 1023 - int(oc.cycle_val(cycle))
        return val

    def calc_hier_ndd_score(self, ndd):
        val = 0
        for oc in self.opt_criteria:
            val <<= 10
            if oc.sense == 'MAX':
                val += min(1023, int(oc.altruist_val(ndd)))
            else:
                val += 1023 - int(oc.altruist_val(ndd))
        return val

    def solve(self, max_solutions, invert_edges):
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

        for c, node_id in zip(self.chains, chain_node_ids):
            for pd_pair in c.pd_pairs:
                patient_to_nodes[pd_pair.patient].append(node_id)
                paired_donor_to_node[pd_pair.donor].append(node_id)
            ndd_to_node[c.altruist_edge.altruist].append(node_id)
            hier_scores[node_id] = self.calc_hier_chain_score(c)

        for c, node_id in zip(self.cycles, cycle_node_ids):
            for pd_pair in c.pd_pairs:
                patient_to_nodes[pd_pair.patient].append(node_id)
                paired_donor_to_node[pd_pair.donor].append(node_id)
            hier_scores[node_id] = self.calc_hier_cycle_score(c)

        for ndd, node_id in zip(altruists, unused_ndd_node_ids):
            ndd_to_node[ndd].append(node_id)
            hier_scores[node_id] = self.calc_hier_ndd_score(ndd)

        adj_mat = [[False]*num_nodes for i in range(num_nodes)]
        
        for d in [patient_to_nodes, paired_donor_to_node, ndd_to_node]:
            for key, val in d.iteritems():
                self.add_clique(adj_mat, val)

#        for row in adj_mat:
#            print " ".join("X" if x else "." for x in row)

        for i in range(num_nodes-1):
            for j in range(i+1, num_nodes):
                edge_exists = adj_mat[i][j]
                if invert_edges:
                    edge_exists = not edge_exists
                if edge_exists:
                    print "e", i+1, j+1

        for i, score in enumerate(hier_scores):
            print "n", i+1, score

