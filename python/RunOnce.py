"""
A script to do a standard run of a yamlfile, print a table of parameters,
and if asked save the sim to a pickle file.

Usage:
RunOnce.py yamlfile

Optional argument:
    --pkl_out='output_pickle_filename.pkl'

Run default set of optical, then bolo calcs
on instrument (and etc) described in yamlfile.

"""

import numpy as np
import yaml
import pickle
import argparse
from jbolo_funcs import *

# Parse the arguments
parser = argparse.ArgumentParser(description='Run baseline model for a given yaml file.')
parser.add_argument('expt_yaml', type=str, help='yaml file describing experiment')
parser.add_argument(
    '--pkl_out',
    default=None,
    type=str,
    help='Output pickle filename'
)
args = parser.parse_args()

# Load the yaml file into the sim dictionary, where we'll also put all the results.
sim = yaml.safe_load(open(args.expt_yaml))

# Run optical calculations all the way up through the bolometer absorption, ie P_optical.
run_optics(sim)

# Run bolometer calculations given the P_optical for each channel.
# If you don't run_optics first, make sure you load values into sim['results'][channel][P_optical].
#
run_bolos(sim)

# Print the table of NET/etc table
print_full_table(sim)

# Write it all to a pickle file
#  (fwiw, the "with open" syntax automatically cleans up and closes the file)
if (args.pkl_out != None):
    with open(args.pkl_out,'wb') as f:
        pickle.dump( sim, f )

# If you need to read the pickle file back in, do:
#with open(filename, 'rb') as f:
#    sim2 = pickle.load(f)
