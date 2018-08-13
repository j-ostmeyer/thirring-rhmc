#!/usr/bin/python
'''
Fake mode: activated adding the word 'fake' to the argument list, uses MPI=fake
           in MkFlags.
Detailed mode: activated adding the word 'detailed' to the argument list,
               instead of exploring only the hypercubic-symmetric partitionings, it
               tries all the partitionings that give a different value of total ranks.
               Among all the partitionings that give the same value of the total ranks,
               the one having the least variance for the set (NP_X,NP_Y,NP_T) is chosen.

NOTE0: A file called 'benchsetup_parameters.py' must be present, and must contain the 
       definition of the parameters necessary to the script to work. A possible content is:
#####
ranks_per_node = 32 # 
maxnranks = 16
ksizes = [4,6,8,10,12,16]
####
       This is just to keep parameters out of the version control system.

NOTE1: modes must be specified in all three phases, that is 'createdir', 'prepare' and 'run'.

NOTE2: This script has many flaws. But it was built "on the way" without a precise idea
       of the final destination. E.g. : violation of the interface segregation principle.

NOTE3: The 'run' function must not be used on login nodes. 

'''
from os import system,listdir,walk
from sys import argv
import os
import glob
from contextlib import contextmanager
import math
from benchsetup_parameters import * # to keep parameters out ov VCS

setFake = False
detailed_mode = False
@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def create_work_in_dir(func):
    def work_in_dir(size,directory,possible_subsizes):
        with cd(directory):
            #os.system('make clean')
            print directory
            
            for divs,ssize in possible_subsizes:
                func(divs,ssize)
    return work_in_dir

def cycle(func):
     ksizes = [4,6,8,10,12,16]
     for size in ksizes:
        if setFake:
            directory = 'benchmarks_fake'+str(size)
        else:
            directory = 'benchmarks'+str(size)
        
        possible_subsizes = []

        if detailed_mode:
            import numpy as np
            from itertools import combinations_with_replacement as cwp
            possible_divs = []
            for div in range(1,size):
                if size%div is 0:
                    possible_divs.append(div)
            best_choice_for_given_nrank = {}
            for nps in cwp(possible_divs,3):
                totranks = nps[0]*nps[1]*nps[2]
                var = np.var(nps) 
                if totranks in best_choice_for_given_nrank:
                    oldvar,oldnps = best_choice_for_given_nrank[totranks]
                    if oldvar > var:
                        best_choice_for_given_nrank[totranks] = \
                             var,nps
                else:
                    best_choice_for_given_nrank[totranks] = var,nps
            
            for nranks,(var,nps) in best_choice_for_given_nrank.items():
                divs = tuple(np.array(size/np.array(nps),dtype=int))
                possible_subsizes.append((nps,divs))

        else:
            for divs in range(1,size):
                if size%divs is 0:
                    possible_subsizes.append((divs,size/divs))
        
        print size,possible_subsizes
        func(size,directory,possible_subsizes)


def createdir(size,directory,possible_subsizes):
    with open('benchmarks/benchmark_params.tmpl.F90','r') as f:
        benchparam_tmpl_text = f.read()

    benchparam_text = benchparam_tmpl_text.replace('SEDKSIZE',str(size))
    script='''
mkdir -p {directory}
cd {directory} 
ln -s ../benchmarks/benchmark_congrad.F90 ./
ln -s ../benchmarks/benchmark_qmrherm_1.F90 ./
cp ../benchmarks/MkFlags.tmpl ./MkFlags
cp ../benchmarks/Makefile ./

'''.format(directory = directory)
    if setFake:
        fakecorrection = "\nsed -i 's/MPI=yes/MPI=fake/' ./MkFlags\n"
        script += fakecorrection

    scriptname = 'script'
    with open(scriptname,'w') as f:
        f.write(script)
    os.system('bash ./'+scriptname)
    
    with open(os.path.join(directory,'benchmark_params.F90'),'w') as f:
         f.write(benchparam_text)


def prepare(divs,ssize):
    if detailed_mode:
        divX,divY,divT = divs
        script='''
make clean
make -j8 benchmark_congrad benchmark_qmrherm_1 NP_X={divX} NP_Y={divY} NP_T={divT}
NEWDIR={divX}x{divY}x{divT}
mkdir -p $NEWDIR
cp benchmark_congrad benchmark_qmrherm_1 $NEWDIR
'''.format(divX = divX,divY = divY, divT = divT)

    else:
        script='''
make clean
make -j8 benchmark_congrad benchmark_qmrherm_1 NP_X={div} NP_Y={div} NP_T={div}
NEWDIR={div}x{div}x{div}
mkdir -p $NEWDIR
cp benchmark_congrad benchmark_qmrherm_1 $NEWDIR
'''.format(div = divs)

    scriptname = 'script'
    print script
    with open(scriptname,'w') as f:
        f.write(script)
    os.system('bash ./'+scriptname)
    

def get_nranks(divs):
    if detailed_mode:
        divX,divY,divT = divs
        nranks=divX*divY*divT
    else:
        nranks=divs**3
    
    return nranks

def get_newdirname(divs):
    if detailed_mode:
        divX,divY,divT = divs
        newdirname="{divX}x{divY}x{divT}".format(divX = divX,divY = divY, divT = divT)
    else:
        newdirname="{div}x{div}x{div}".format(div=divs)
    return newdirname


def write_runscripts(divs,ssize):
    newdirname = get_newdirname(divs)
    nranks = get_nranks(divs)
    script='''
/usr/bin/time -o timeqmr -p mpirun -n {nrank} ./benchmark_qmrherm_1 > qmroutput
/usr/bin/time -o timecongrad -p mpirun -n {nrank} ./benchmark_congrad > congradoutput
'''.format(div = divs, nrank = nranks)

    suffix = str(int(math.ceil(float(nranks)/ranks_per_node)))
    scriptname = 'scriptrun'+suffix
    print script
    with open(os.path.join(newdirname,scriptname),'w') as f:
        f.write(script)
 
def run(divs,ssize):
    newdirname = get_newdirname(divs)
    nranks = get_nranks(divs)
    suffix = str(int(math.ceil(float(nranks)/ranks_per_node)))
    scriptname = 'scriptrun'+suffix
    if nranks <= maxnranks:
        with cd(newdirname):
            os.system('bash ./'+scriptname)
           
 
modes = {
    'createdir' : createdir,
    'prepare': create_work_in_dir(prepare),
    'runscripts' : create_work_in_dir(write_runscripts),
    'run' : create_work_in_dir(run)
} 


if __name__ == "__main__" : 
    if 'fake' in argv:
        setFake = True
        argv.remove('fake')
    if 'detailed' in argv:
        detailed_mode = True
        argv.remove('detailed')
 
    cycle(modes[argv[1]])

    


