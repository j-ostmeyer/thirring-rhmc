#!/usr/bin/python
from os import system,listdir,walk
from sys import argv
import os
import glob
from contextlib import contextmanager
import math

ranks_per_node = 32

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def work_in_dir(func,directory,possible_subsizes):
    with cd(directory):
        #os.system('make clean')
        print directory
        
        for divs,ssize in possible_subsizes:
            func(divs,ssize)

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
        directory = 'benchmarks'+str(size)
        possible_subsizes = []
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
    scriptname = 'script'
    with open(scriptname,'w') as f:
        f.write(script)
    os.system('bash ./'+scriptname)
    
    with open(os.path.join(directory,'benchmark_params.F90'),'w') as f:
         f.write(benchparam_text)


def prepare(divs,ssize):
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
    

def run(divs,ssize):
    nranks=divs**3
    suffix = str(int(math.ceil(float(nranks)/ranks_per_node)))
    newdirname="{div}x{div}x{div}".format(div=divs)
    script='''
/usr/bin/time -o timeqmr -p mpirun -n {nrank} ./benchmark_qmrherm_1 > qmroutput
/usr/bin/time -o timecongrad -p mpirun -n {nrank} ./benchmark_congrad > congradoutput
'''.format(div = divs, nrank = nranks)
    scriptname = 'scriptrun'+suffix
    print script
    with open(os.path.join(newdirname,scriptname),'w') as f:
        f.write(script)
    #os.system('bash ./'+scriptname)
 
            
modes = {
    'createdir' : createdir,
    'prepare': create_work_in_dir(prepare),
    'run' : create_work_in_dir(run)
} 


if __name__ == "__main__" : 
    cycle(modes[argv[1]])


    


