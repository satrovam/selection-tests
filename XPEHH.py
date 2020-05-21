#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 15:02:52 2020

@author: mari
"""
# load needed libraries
import allel
import zarr
import argparse
import glob
import os
from datetime import datetime
import numpy as np
import pandas as pd
import sys

# load modules
import logger
import preprocessing
import xp_utils


################## PARSING ARGUMETNS ####################################################

# create ArgumentParser object
parser = argparse.ArgumentParser(description='  \
                                 \nExample: ')

# add positional argumimport pandas as pdent ZARR DIRR
parser.add_argument("zarr_dir", help='Path to directory with zarr directories for chromosomes of interest')
# add positional argument RESULT DIR
parser.add_argument("result_dir", help='Path to dir with results')
# add positional argument PANEL_FILE_PATH
parser.add_argument("panel_file_path", help='Path to panel file that contains information \
                    about each individual (population, subpopulation, sex) file can be downloaded from 1000Genomes ftp')
# add optional argument
parser.add_argument('-chromlist', help='List of chrmosomes to use for the computation; default - all chromosomes',
  nargs="*",  # 0 or more values expected => creates a list
  type=str,
  default=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',\
           '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'],  # default if nothing is provided
)

parser.add_argument('-exceptchrom', help='List of chrmosomes to exclude from the computation; default - empty list',
  nargs="*",  # 0 or more values expected => creates a list
  type=str,
  default=[],  # default if nothing is provided
)

args = parser.parse_args()

chroms = [chrom for chrom in args.chromlist if chrom not in args.exceptchrom]

# get all the paths of zarrs given the chromosome names
zarr_paths = []
for chrom in chroms:
    zarr_path = glob.glob(args.zarr_dir +'*chr' + chrom + '[!0123456789]*.zarr')
    if len(zarr_path) != 0:
        zarr_paths.append(zarr_path[0])


######################### CHECK WITH USER + INITIATE LOG ####################


# check with the user that everything is ok
print("YOU ARE ABOUT TO COMPUTE THE XPEHH STATISTIC ON:")
for zarr_path in zarr_paths:
    print(zarr_path)
print("\nRESULTS WILL BE SAVED IN DIRECTORY:")
print(args.result_dir)

# if everything ok, proceed with the computation
proceed = input("PROCEED? [y/n]: ")
if proceed != 'y':
    sys.exit()


# check if result dir is ended with a slash
if args.result_dir[-1] != "/":
    args.result_dir = args.result_dir + "/"
    
# create result directory if it doesn't already exist
if not os.path.exists(args.result_dir):
    os.makedirs(args.result_dir)
    
# create log file, write header with time stamp, script name and parameters that were used
logger.create_log(args, zarr_paths)

# load the csv description of samples in populations
panel = pd.read_csv(args.panel_file_path, sep='\t', usecols=['sample', 'pop', 'super_pop'])


########################### COMPUTE  ######################################### 

# create an array with all possible combinations of 2 populations
pop_pairs = xp_utils.create_pop_pairs(panel)
      
 
for zarr_path in zarr_paths:
    root_name = zarr_path.split("/")[-1].split(".") # get the name of zarr
    root_name.pop() # pop the extension at the end (.zarr)
    root_name = ".".join(root_name)  
    result_path = args.result_dir + root_name + ".xpehh.tsv"
    
    # UPDATE LOG
    logger.update_log("loading data from: " + zarr_path)
        
    # acess the zarr
    callset = zarr.open_group(zarr_path, mode='r')
    
    # get the filtered gt array and positions of filtered variants
    gt, positions = preprocessing.filter_by_AF(callset, 0.05)

    logger.update_log("Dimensions of genotype data (variants, samples, ploidy): " + " ".join(map(str,gt.shape)))
    
    
    # check whether the ordering of samples in zarr callset is equal to the order in csv
    samples = callset['samples'][:]
    if np.all(samples == panel['sample'].values):
        logger.update_log("Order of samples ok")
    else:
        logger.update_log("Order of samples in panel file does not match order of samples in given zarr. \
              It is possible that you are using wrong panel file path \
              e.g. from different phase than you variant data comes from different phase than your data")
        continue
    
    
    df = pd.DataFrame({"variant_pos" : positions})
    # get the name of chromosome from zarr file
    chrom = root_name.split("chr")[1].split("_")[0]
    # insert chromosome name as the fist column
    df.insert(0, 'chromosome', chrom)
    
    for pair in pop_pairs:
           
        ht_pop1 = xp_utils.get_haplotypes(gt, panel, pair[0])
        ht_pop2 = xp_utils.get_haplotypes(gt, panel, pair[1])
        
        
    
        logger.update_log("computing XPEHH for pair " + pair[0] + " " + pair[1])
        
        
        logger.update_log("dimensions of haplotype data for pop " + pair[0] + ": " + " ".join(map(str,ht_pop1.shape)))
        logger.update_log("dimensions of haplotype data for pop " + pair[1] + ": " + " ".join(map(str,ht_pop2.shape)))
        logger.update_log("dimensions of positions: " + str(len(positions)))
        
        
        result = allel.xpehh(h1=ht_pop1, h2=ht_pop2, pos=positions, map_pos=None, min_ehh=0.05,
                    include_edges=False, gap_scale=20000, max_gap=200000,
                    is_accessible=None, use_threads=True)
        
        df["xpehh"] = result
        
###################### WRITE RESULTS INTO TSV #######################################
        
        pops_res_dir = args.result_dir + pair[0] + '_' + pair[1] + '/'
        if not os.path.exists(pops_res_dir):
            os.makedirs(pops_res_dir)
            
        result_path = pops_res_dir + root_name + ".xpehh.tsv"
        df.to_csv(result_path, index=False, sep='\t')
         
        # UPDATE LOG
        logger.update_log(datetime.now().strftime('%H:%M:%S') + " - Computation completed")
        logger.update_log(" Resuts saved into: " + result_path)
        
        
        
        

        
        
