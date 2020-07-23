#!/usr/bin/env python3
'''
This is supposed to be a slimmed version of 
benchsetup.py.
'''
from sys import argv,exit
import pandas as pd
import os
import shutil
import re
from contextlib import contextmanager
import subprocess

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

  
def sed(infilename,outfilename,pardict):
    '''
    Replace in the file `infilename` the values from pardict that have been
    prefixed with 'SED'
    '''

    print(f"Reading {infilename}")
    with open(infilename,'r') as f:
        text = f.read()

    def replace_keep_spaces(match,v):
        groups = match.groups()
        pre = groups[0]
        post = groups[1]
        return pre + str(v) + post


    print(f"In {infilename}:")
    for k,v in pardict.items():
        text,n = re.subn(r'(\W)SED'f'{k}'r'(\W)',
                         lambda match: replace_keep_spaces(match,v),
                         text)
        print(f"Found SED{k} {n} times")

    print(f"Writing {infilename}")
    with open(outfilename,'w') as f:
        f.write(text)

def dirname(row):
    return '_'.join([f'{k}{v}' for k,v in row._asdict().items()])

def copyall(benchmarks_lib_location,destination):
    for filename in [ 'MkFlags.tmpl',
                      'run.tmpl.slurm',
                      'benchmark_params.tmpl.F90',
                      'Makefile.tmpl',
                      'benchmark_congrad.F90',
                      'benchmark_full_md.F90',
                      'benchmark_qmrherm_1.F90',
                      'benchmark_qmrherm_1_sp.F90']:
        shutil.copy(os.path.join(benchmarks_lib_location,filename),
                    os.path.join(destination,filename))

def adjust_files(directory,row):
    for filename in ['MkFlags.tmpl', 
                     'benchmark_params.tmpl.F90',
                     'run.tmpl.slurm']:
        infilename = os.path.join(directory,filename)
        outfilename = infilename.replace('.tmpl','')
        nprocs = (row.NP_X *
                  row.NP_Y *
                  row.NP_T *
                  row.NP_THIRD)

        pardict = dict(row._asdict())
        pardict.update({'NPROCS':nprocs})
        sed(infilename,outfilename,pardict)

def adjust_makefile(directory,topdir):
    filename = 'Makefile.tmpl'
    infilename = os.path.join(directory,filename)
    outfilename = infilename.replace('.tmpl','')
    nup = len(directory.split('/'))
    topdir_pos = os.path.join(*(['..']*nup + [topdir]))
 
    pardict = {'TOPDIR':topdir_pos}
    sed(infilename,outfilename,pardict)


if __name__ == '__main__':
    # read setup file
    try:
        setup = pd.read_csv(argv[1], sep = '\s+')
        benchmarks_lib_location = argv[2]   
        mkrules = argv[3]
    except:
        print(f'Usage: {argv[0]} setup_file benchmark_lib_dir mkrules')
        exit(1)

    for row in setup.itertuples():
        newdir = dirname(row)

        # creating
        os.makedirs(newdir,exist_ok = True)
        # copying files
        copyall(benchmarks_lib_location,newdir)
        # adjusting them
        adjust_files(newdir,row)
        # pointing the makefile to the right MkRules
        adjust_makefile(newdir,mkrules)
        #with cd(newdir):
        #    proc = subprocess.run("make", capture_output=True)
        #    print(proc.stdout)
        #    print(proc.stderr)
            



 


