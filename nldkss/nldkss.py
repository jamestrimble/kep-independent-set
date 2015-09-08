import argparse
import sys
import json
from ortools.linear_solver import pywraplp
from nldkss_pool import *
from optimality_criteria import *
import pool_reader

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

    def _optimise(self, s, opt_criterion):
        z = s.Sum(
            [chain.mip_var*opt_criterion.chain_val(chain) for chain in self.chains] +
            [cycle.mip_var*opt_criterion.cycle_val(cycle) for cycle in self.cycles] +
            [altruist.mip_var*opt_criterion.altruist_val(altruist)
                        for altruist in self.pool.altruists])
            
        if opt_criterion.sense=='MAX':
            s.Maximize(z)
        else:
            s.Minimize(z)

        solve_status = s.Solve()

        return z, solve_status

    def _create_solver(self):
        return pywraplp.Solver('CoinsGridCLP',
                                 pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    def _create_vars_and_constraints(self, s, patients, paired_donors, altruists):
        for person in patients + paired_donors + altruists:
            person.mip_vars = []
        
        for chain in self.chains:
            chain.mip_var = s.IntVar(0, 1, 'chain' + str(chain.index))
        for cycle in self.cycles:
            cycle.mip_var = s.IntVar(0, 1, 'cycle' + str(cycle.index))
        for altruist in altruists:
            altruist.mip_var = s.IntVar(0, 1, 'altruist' + str(altruist.nhs_id))
            altruist.mip_vars.append(altruist.mip_var)

        for chain in self.chains:
            chain.altruist_edge.altruist.mip_vars.append(chain.mip_var)
            self._add_var_to_patients_and_donors(chain.pd_pairs, chain.mip_var)
        for cycle in self.cycles:
            self._add_var_to_patients_and_donors(cycle.pd_pairs, cycle.mip_var)

        for person in patients + paired_donors:
            s.Add(s.Sum(person.mip_vars) <= 1)

        for altruist in altruists:
            s.Add(s.Sum(altruist.mip_vars) == 1)

    def _enforce_objective(self, s, z, obj_val, sense):
        if sense=='MAX':
            s.Add(z >= obj_val)
        else:
            s.Add(z <= obj_val)

    def _items_in_optimal_solution(self, items):
        return [item for item in items if item.mip_var.SolutionValue() > 0.5]

    def solve(self, max_solutions):
        patients = self.pool.patients
        paired_donors = self.pool.paired_donors
        altruists = self.pool.altruists

        s = self._create_solver()

        self._create_vars_and_constraints(s, patients, paired_donors, altruists)

        for opt_criterion in opt_criteria[:-1]:
            z, status = self._optimise(s, opt_criterion)
            if status != s.OPTIMAL:
                raise OptimisationException("Solver status was " + str(solve_status))
            self._enforce_objective(s, z, s.Objective().Value(), opt_criterion.sense)

        
        n_solutions = 0
        best_objval_found = -1
        while True:
            z, solve_status = self._optimise(s, opt_criteria[-1])

            if solve_status != s.OPTIMAL:
                raise OptimisationException("Solver status was " + str(solve_status))

            objval = s.Objective().Value()
            if objval+self.EPSILON < best_objval_found:
                break

            optimal_chains = self._items_in_optimal_solution(self.chains)
            optimal_cycles = self._items_in_optimal_solution(self.cycles)
            optimal_altruists = self._items_in_optimal_solution(altruists)
            optimal_vars = [o.mip_var for o in optimal_chains + optimal_cycles + optimal_altruists]

            best_objval_found = objval

            n_solutions += 1
            
            for item in optimal_chains + optimal_cycles + optimal_altruists:
                print str(item)
            print

            if n_solutions==max_solutions:
                return best_objval_found, n_solutions, True

            if len(optimal_vars)==0:
                break

            # Ensure that this optimal set of cycles, chains and unused altruists isn't re-found
            s.Add(s.Sum(optimal_vars) <= len(optimal_vars)-1)


        return best_objval_found, n_solutions, False

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Hierarchical kidney-exchange optimisation")
    parser.add_argument("-f", "--file", help="Input file name", type=str,
                        required=True)
    parser.add_argument("-c", "--criteria",
                        help="A colon-separated list of optimality criteria," +
                             "such as effective:size:3way:backarc:weight",
                        type=str,
                        required=True)
    parser.add_argument("-e", "--cycle",
                        help="Maximum cycle length",
                        type=int,
                        required=True)
    parser.add_argument("-n", "--chain",
                        help="Maximum chain length",
                        type=int,
                        required=True)
    parser.add_argument("-m", "--max",
                        help="Stop after finding this number of solutions",
                        type=int,
                        default=100) 
    args = parser.parse_args()

    opt_criteria = get_criteria(args.criteria)

    with open(args.file) as json_file:
        pool = pool_reader.read(json.load(json_file)["data"])
        pool_optimiser = PoolOptimiser(pool, opt_criteria, args.cycle, args.chain)
        objval, n_solutions, reached_max = pool_optimiser.solve(args.max)
        print "Objective value: {} Number of solutions: {} Reached limit: {}".format(
                objval, n_solutions, "TRUE" if reached_max else "FALSE")
