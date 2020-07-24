#!/usr/bin/env python3
import benchmarks_lib as bl
from sys import argv
import pandas as pd
import os

if __name__ == '__main__':
    # read setup file
    try:
        setup = pd.read_csv(argv[1], sep = '\s+')
        benchmarks_lib_location = argv[2]   
        mkrules = argv[3]
    except Exception as e:
        print(e)
        print(f'Usage: {argv[0]} setup_file benchmark_lib_dir mkrules')
        exit(1)

    for row in setup.itertuples():
        newdir = bl.dirname(row)

        # creating
        os.makedirs(newdir,exist_ok = True)
        # copying files
        bl.copyall(benchmarks_lib_location,newdir)
        # adjusting them
        bl.adjust_files(newdir,row)
        # pointing the makefile to the right MkRules
        bl.adjust_makefile(newdir,mkrules)
        #with cd(newdir):
        #    proc = subprocess.run("make", capture_output=True)
        #    print(proc.stdout)
        #    print(proc.stderr)
            



 


