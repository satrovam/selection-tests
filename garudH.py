# load needed libraries
import allel
import zarr
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
import time

# load modules
import logger
import window_utilities

################## PARSING ARGUMETNS ####################################################

# create ArgumentParser object
parser = argparse.ArgumentParser(description='Compute the H1, H12, H123 and H2/H1 tests \
                                 for detecting signatures of soft sweeps, as defined in Garud et al. (2015) \
                                 on one or all chromosomes in zarr directory. \
                                 \nExample: python3 garudH.py data/zarr/ results/garudH/ -size 30000 -step 3000 -exceptchrom X Y')

# add positional argumimport pandas as pdent ZARR DIRR
parser.add_argument("zarr_dir", help='Path to directory with zarr directories for chromosomes of interest')
# add positional argument RESULT DIR
parser.add_argument("result_dir", help='Path to dir with results')
# add optional argument WINDOW SIZE
parser.add_argument('-size', help='Size of window in bp, default = 10000')
# add optional argument WINDOW STEP
parser.add_argument('-step', help='The number of variants between start positions of windows.\
                                    If not given, defaults to the window size, i.e., windows will not overlap')

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

if args.size:
    size = int(args.size)
else:
    size = 10000
if args.step:
    step = int(args.step)
else:
    step = None

chroms = [chrom for chrom in args.chromlist if chrom not in args.exceptchrom]

# get all the paths of zarrs given the chromosome names
zarr_paths = []
for chrom in chroms:
    zarr_path = glob.glob(args.zarr_dir +'*chr' + chrom + '[!0123456789]*.zarr')
    if len(zarr_path) != 0:
        zarr_paths.append(zarr_path[0])
    
######################### CHECK WITH USER + INITIATE LOG #####################################


# check with the user that everything is ok
print("YOU ARE ABOUT TO COMPUTE THE H1, H12, H123 AND H2/H1 STATISTICS ON:")
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
        
        
########################### COMPUTE GARUD'S H ##########################################3           
# start timer
t0 = time.time()

for zarr_path in zarr_paths:
    root_name = zarr_path.split("/")[-1].split(".") # get the name of zarr
    root_name.pop() # pop the extension at the end (.zarr)
    root_name = ".".join(root_name)  
    # update log and print to stdout
    logger.update_log("loading data from: " + zarr_path)    
    
    # acess the zarr
    callset = zarr.open_group(zarr_path, mode='r')
    
    # get the filtered gt array and positions of filtered variants
    gt, positions = preprocessing.filter_by_AF(callset, 0.05)
    
    # UPDATE LOG
    logger.update_log("\nDimensions of filtered genotype data \ (variants, samples, ploidy): " + " ".join(map(str,gt.shape))) 
    
    # if callset too big, load the genotype as chunked array
    #gt = allel.GenotypeChunkedArray(gt_zarr)
    #print("Dimensions of genotype data (variants, samples, ploidy): ", gt.shape)
    
    # haplotype - needed for garud h
    ht = gt.to_haplotypes()
    
    logger.update_log("Computing Garud's H over variant windows")
    h1, h12, h123, h2_h1 = allel.moving_garud_h(ht, size=size, start=0, stop=None, step=step)
      
    # wndows generator     
    windows = window_utilities.index_windows(positions, size=size, step=step, start=0, stop=None)
    
    # declare numpy arrays that will store start end end positions of windows
    window_start = np.array([], dtype=np.uint32) 
    window_end = np.array([], dtype=np.uint32)
        
    # generate the window positions using windows generator
    for i, j in windows:
        window_start = np.append(window_start, positions[i])
        window_end = np.append(window_end, positions[j-1])
          
    #check lengths
    if len(window_start) != len(h1) or len(window_end) != len(h1):
        logger.update_log("The number of windows differs from the number of computed \
              values for windows! Terminating computation for this chromosome")  
        break
    
    # stop timer
    t1 = time.time()
    logger.update_log('Computation completed in ' + (str(round(t1-t0,2)) + 's'))
    
###################### WRITE RESULTS INTO TSV #######################################
    # get the name of the chromosome
    chrom = root_name.split("chr")[1].split("_")[0]
    
    # H1  
    df = pd.DataFrame({"chromosome" : chrom, "window_start" : window_start, \
                       "window_end" : window_end, "h1" : h1})
    
    result_path = args.result_dir + 'garudH1/' + root_name + ".garudH1.tsv"
    df.to_csv(result_path, index=False, sep='\t')
    logger.update_log("Resuts saved into: " + result_path)
    
    # H12  
    df = pd.DataFrame({"chromosome" : chrom, "window_start" : window_start, \
                       "window_end" : window_end, "h12" : h12})
    
    result_path = args.result_dir + 'garudH12/'+ root_name + ".garudH12.tsv"
    df.to_csv(result_path, index=False, sep='\t') 
    logger.update_log("Resuts saved into: " + result_path)
    
    # H123  
    df = pd.DataFrame({"chromosome" : chrom, "window_start" : window_start, \
                       "window_end" : window_end, "h123" : h123})
    
    result_path = args.result_dir + 'garudH123/' + root_name + ".garudH123.tsv"
    df.to_csv(result_path, index=False, sep='\t')  
    logger.update_log("Resuts saved into: " + result_path)
    
    # H2_H1   
    df = pd.DataFrame({"chromosome" : chrom, "window_start_bp" : window_start, \
                       "window_end_bp" : window_end, "h2_h1" : h2_h1})
    
    result_path = args.result_dir + 'garudH2_H1/' + root_name + ".garudH2_H1.tsv"
    df.to_csv(result_path, index=False, sep='\t') 
    
    logger.update_log("Resuts saved into: " + result_path)
