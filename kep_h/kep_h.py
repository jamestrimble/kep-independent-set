import argparse
import sys
import json
from kep_h_pool import *
from optimality_criteria import *
from kep_h_pool_optimiser import PoolOptimiser
import pool_reader

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

    if args.file.endswith(".json"):
        with open(args.file) as json_file:
            pool = pool_reader.read(json.load(json_file)["data"])
            pool_optimiser = PoolOptimiser(pool, opt_criteria, args.cycle, args.chain)
            objval, n_solutions, reached_max = pool_optimiser.solve(args.max)
            print "Objective value: {} Number of solutions: {} Reached limit: {}".format(
                    objval, n_solutions, "TRUE" if reached_max else "FALSE")
    else:
        print "Input file must be in JSON format"
