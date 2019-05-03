#!/usr/bin/python
'''
This script creates and configures a number of benchmarks, compiling different 
versions of the benchmarks for different geometries. The workflow happens by 
 - first creating the directories, 
 - then preparing the content of each directory,
 - then writing the scripts
 - optionally run them locally.
If scripts are created, then we want to execute them together in a batch script.
See benchmark1n.sh for this.

It must be called with a 'stepname' as argument, in the directory 
above this one. Possible stepnames are
written in the 'modes' dictionary, and should be invoked sequentially.
Read all the notes.

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

NOTE1: modes must be specified in all three phases, that is 'createdir', 'prepare', 
       'writescripts' and 'run'.

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
    '''
    This function creates a function that consists of the argument
    run in a group of directories.
    '''
    def work_in_dir(size,directory,possible_subsizes):
        with cd(directory):
            #os.system('make clean')
            print directory
            
            for divs,ssize in possible_subsizes:
                func(divs,ssize)
    return work_in_dir

def cycle(func,tag):
     ksizes = [4,6,8,10,12,16]
     for size in ksizes:
        if setFake:
            directory = tag + '_benchmarks_fake'+str(size)
        else:
            directory = tag + '_benchmarks'+str(size)
        
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
ln -s ../benchmarks/benchmark_qmrherm_split1.F90 ./
ln -s ../benchmarks/benchmark_qmrherm_split_nodir1.F90 ./
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
make -j8 benchmark_congrad NP_X={divX} NP_Y={divY} NP_T={divT}
make -j8 benchmark_qmrherm_1 NP_X={divX} NP_Y={divY} NP_T={divT}
make -j8 benchmark_qmrherm_split1 NP_X={divX} NP_Y={divY} NP_T={divT}
make -j8 benchmark_qmrherm_split_nodir1 NP_X={divX} NP_Y={divY} NP_T={divT}
NP_X={divX} NP_Y={divY} NP_T={divT}
NEWDIR={divX}x{divY}x{divT}
mkdir -p $NEWDIR
cp benchmark_congrad $NEWDIR
cp benchmark_qmrherm_1 $NEWDIR
cp benchmark_qmrherm_split1 $NEWDIR
cp benchmark_qmrherm_split_nodir1 $NEWDIR

'''.format(divX = divX,divY = divY, divT = divT)

    else:
        script='''
make clean
make -j8 benchmark_congrad NP_X={div} NP_Y={div} NP_T={div}
make -j8 benchmark_qmrherm_1 NP_X={div} NP_Y={div} NP_T={div}
make -j8 benchmark_qmrherm_split1 NP_X={div} NP_Y={div} NP_T={div}
make -j8 benchmark_qmrherm_split_nodir1 NP_X={div} NP_Y={div} NP_T={div}

NEWDIR={div}x{div}x{div}
mkdir -p $NEWDIR
cp benchmark_congrad $NEWDIR
cp benchmark_qmrherm_1 $NEWDIR
cp benchmark_qmrherm_split1 $NEWDIR
cp benchmark_qmrherm_split_nodir1 $NEWDIR

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
/usr/bin/time -o timecongrad -p mpirun -n {nrank} ./benchmark_congrad > congradoutput
/usr/bin/time -o timeqmr -p mpirun -n {nrank} ./benchmark_qmrherm_1 > qmroutput
/usr/bin/time -o timeqmr_split -p mpirun -n {nrank} ./benchmark_qmrherm_split1 > qmrsplitoutput
/usr/bin/time -o timeqmr_split_nodir -p mpirun -n {nrank} ./benchmark_qmrherm_split_nodir1 > qmrsplitnodiroutput
'''.format(div = divs, nrank = nranks)

    suffix = str(int(math.ceil(float(nranks)/ranks_per_node)))
    scriptname = 'scriptrun'+suffix
    print script
    with open(os.path.join(newdirname,scriptname),'w') as f:
        f.write(script)

def write_runscripts_extrae(divs,ssize):
    newdirname = get_newdirname(divs)
    nranks = get_nranks(divs)
    script='''
module purge
module use /home/s.michele.mesiti/modules
module load extrae-gnu-8.1-mpi3.1.1 
TRACE=$EXTRAE_HOME/share/example/MPI/ld-preload/trace.sh 
ln -s $EXTRAE_HOME/share/example/MPI/extrae_fixed.xml ./extrae.xml
for type in congrad qmrherm_1 qmrherm_split1 qmrherm_split_nodir1
do
(
mkdir -p $type
cd $type
mpirun -n {nrank} $TRACE ../benchmark_$type              > output_$type &&
mpirun -n {nrank} ${{EXTRAE_HOME}}/bin/mpimpi2prv -syn -f TRACE.mpits -o trace.prv
)
done 
'''.format(div = divs, nrank = nranks)

    suffix = str(int(math.ceil(float(nranks)/ranks_per_node)))
    scriptname = 'scriptrun'+suffix
    print script
    with open(os.path.join(newdirname,scriptname),'w') as f:
        f.write(script)
 
def write_runscripts_scalasca(divs,ssize):
    newdirname = get_newdirname(divs)
    nranks = get_nranks(divs)
    script='''
module purge
module load scalasca
for type in congrad qmrherm_1 qmrherm_split1 qmrherm_split_nodir1
do
(
mkdir -p $type
cd $type
scalasca -analyse mpirun -n {nrank} ../benchmark_$type > output_$type 
)
done 
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
    'writescripts' : create_work_in_dir(write_runscripts),
    'writescripts_extrae' : create_work_in_dir(write_runscripts_extrae),
    'writescripts_scalasca' : create_work_in_dir(write_runscripts_scalasca),
    'run' : create_work_in_dir(run)
} 


if __name__ == "__main__" : 
    if 'fake' in argv:
        setFake = True
        argv.remove('fake')
    if 'detailed' in argv:
        detailed_mode = True
        argv.remove('detailed')

   
    try:
       tag = argv[1]
    except:
       print("Usage: benchsetup.py tag step")
 
    try:
       cycle(modes[argv[2]],tag)
    except KeyError:
       print("'{}' step name not recognized.".format(argv[1]))
       print("Possible step names:")
       print(modes.keys())
    except IndexError:
       print("Please specify step and tag.")
       print("Possible step names:")
       print(modes.keys())
      

    


