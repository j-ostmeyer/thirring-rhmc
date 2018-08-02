#!/usr/bin/python
from os import system,listdir,walk
from sys import argv
import os
import glob
from contextlib import contextmanager

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def cycle(func):
     for directory in glob.glob('benchmarks[0-9]*'):
        size = int(directory.replace('benchmarks',''))
        possible_subsizes = []
        for divs in range(1,size):
            if size%divs is 0:
                possible_subsizes.append((divs,size/divs))
    
        print size,possible_subsizes
        with cd(directory):
            #os.system('make clean')
            print directory
            
            for divs,ssize in possible_subsizes:
                func(divs,ssize)
     


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
    script='''
NEWDIR={div}x{div}x{div}
cd $NEWDIR
/usr/bin/time -o timeqmr -p numactl -m 1 mpirun -n {nrank} ./benchmark_qmrherm_1 > qmroutput
/usr/bin/time -o timecongrad -p numactl -m 1 mpirun -n {nrank} ./benchmark_congrad > congradoutput
'''.format(div = divs, nrank = nranks)
    scriptname = 'scriptrun'
    print script
    with open(scriptname,'w') as f:
        f.write(script)
    os.system('bash ./'+scriptname)
 
            
modes = {
    'prepare': prepare,
    'run' : run
} 


if __name__ == "__main__" : 
    cycle(modes[argv[1]])


    


