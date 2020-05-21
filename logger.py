#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 17:33:44 2020

@author: mari
"""
import sys
from datetime import datetime

def create_log(script_args, zarr_dir_paths):
    global log_file
    time = datetime.now().strftime('%d-%m-%Y')
    log_file = script_args.result_dir + time +".log"
    with open(log_file, 'w') as log:
        log.write("script: " + sys.argv[0])
        log.write("\ntime: " + time)
        log.write("\nparameters:\n")
        for arg in vars(script_args):
            log.write(arg + '\t' + str(getattr(script_args, arg)))
            log.write('\n')
            
        log.write("\nCOMPUTING TEST 0N:")
        for zarr_path in zarr_dir_paths:
            log.write("\n" + zarr_path)

def update_log(message):
    # write into log file
    with open(log_file, 'a') as log:
        log.write("\n" + message)
    # print message to stdout
    print(message)