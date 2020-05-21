#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 19:01:58 2020

@author: mari
"""

import numpy as np
import logger

def create_pop_pairs(panel):
    
    populations = np.unique(panel['pop'].values)
    
    pop_pairs = []
    for i in range(len(populations)):
        pop1 = populations[i]
        for j in range(i + 1, len(populations)):
            pop2 = populations[j]
            pop_pairs.append((pop1, pop2))
    
    return pop_pairs

def get_haplotypes(gt_array, panel, pop):
    
    logger.update_log("Dimensions of genotype data (variants, samples, ploidy): " + " ".join(map(str,gt_array.shape)))    
    # get the indices of samples which belong to given population
    indices_pop = panel.index[panel['pop'] == pop]
    # get genotype data belonging only to given population
    gt_pop = gt_array.take(indices_pop, axis=1)
    return gt_pop.to_haplotypes()

def get_pop_triplets():
    ''' This function returns a 3D lists, so in list[i][j][p]:
        i = superpopulation with j triplets
        j = triplet
        p = population in triplet
    '''
    
    # Too tired to code, here are list of pops grouped in superpops
    pops_lists = [['ESN', 'MSL', 'GWD', 'ACB', 'YRI', 'ASW', 'LWK'],
                  ['PEL', 'PUR', 'MXL', 'CLM'],
                  ['CDX', 'KHV', 'CHS', 'JPT', 'CHB'],
                  ['CEU', 'GBR', 'TSI', 'IBS', 'FIN'],
                  ['STU', 'BEB', 'PJL', 'ITU', 'GIH' ]]
    
    
    # create pairs within each superpopulation
    pop_pairs_list = []
    
    for populations in pops_lists:
        pop_pairs = []
        for i in range(len(populations)):
            pop1 = populations[i]
            for j in range(i + 1, len(populations)):
                pop2 = populations[j]
                pop_pairs.append([pop1, pop2])
        # append list of pop pairs within superpopulation to list of other pop_pairs lists
        pop_pairs_list.append(pop_pairs)        
    
    # add outgroup to each pair
    outgroup1 = 'FIN' # give everyone except europeans FIN as outgroup
    outgroup2 = 'CHB' # give europeans Bejing Chinese as outgroup

    for i in range(len(pop_pairs_list)):
        for j in range(len(pop_pairs_list[i])):
            if i == 3: # if we are in european lists
                pop_pairs_list[i][j].append(outgroup2)
            else:
                pop_pairs_list[i][j].append(outgroup1)
                
    #### FINISH IT
    return pop_pairs_list

def get_pop_alllel_counts(gt, panel, pop):
    ''' Returns allel counts for given population
    '''
    # get the indices of samples (individuals) which belong to pop
    indices_pop = panel.index[panel['pop'] == pop]        
    # get genotype data belonging only to pop 
    gt_pop = gt.take(indices_pop, axis=1)
    #get the allel counts for population (input for pbs)
    logger.update_log("Computing allel counts for population " + pop)
    ac = gt_pop.count_alleles()
    return ac
    
            
    
    
