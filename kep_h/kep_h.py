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
    parser.add_argument("-o", "--node-order",
                        help="Node order (0=default, 1=random, 2=?)",
                        type=int,
                        default=0)
    parser.add_argument("-r", "--reduce-nodes",
                        help="Remove some nodes that can't be part of a solution",
                        action='store_true') 
    parser.add_argument("-i", "--invert-edges",
                        help="Create complement graph",
                        action='store_true') 
    args = parser.parse_args()

    opt_criteria = get_criteria(args.criteria)

    if args.file.endswith(".json"):
        with open(args.file) as json_file:
            pool = pool_reader.read(json.load(json_file)["data"])
            pool_optimiser = PoolOptimiser(pool, opt_criteria, args.cycle, args.chain)
            pool_optimiser.solve(args.invert_edges, args.reduce_nodes, args.node_order)
    else:
        print "Input file must be in JSON format"
