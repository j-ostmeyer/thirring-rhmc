'''
This is supposed to be a slimmed version of 
benchsetup.py.
'''
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


benchmarks_sources = [ 'benchmark_congrad.F90',
                       'benchmark_full_md.F90',
                       'benchmark_qmrherm_1.F90',
                       'benchmark_qmrherm_1_sp.F90']

run_tmpl = 'run.tmpl.slurm'

def copyall(benchmarks_lib_location,destination):
    for filename in [ 'MkFlags.tmpl',
                      run_tmpl,
                      'benchmark_params.tmpl.F90',
                      'Makefile.tmpl'] + benchmarks_sources:
        shutil.copy(os.path.join(benchmarks_lib_location,filename),
                    os.path.join(destination,filename))

def adjust_files(directory,row):
    for filename in ['MkFlags.tmpl', 
                     'benchmark_params.tmpl.F90',
                     run_tmpl]:
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

def get_benchmark_run_order(benchmark_lib_location):
    '''
    We look at the order of executables in the run template.
    This is also the order of the timing in the slurm output files.
    '''
    run_tmpl_path = os.path.join(benchmark_lib_location, run_tmpl)
    with open(run_tmpl_path,'r') as f:
        benchmarks_ordered = [ benchmark 
                           for line in f.readlines()
                           for bsource in benchmarks_sources
                           if re.search(r'\W'+
                                       (benchmark:=bsource.replace('.F90',''))+
                                       r'\W',line)]

    return benchmarks_ordered

def parse_timing_file(filename):

    with open(filename,'r') as f:
        timings = [ float(line.split()[-1])
                    for line in f.readlines()
                    if 'Time per iteration:' in line ]

    return timings


