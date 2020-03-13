# coding: utf-8
import os
import itertools as it

# getting all filenames
sources = []
for f in os.listdir('.'):
    if '.f90' in f.lower():
        sources.append(f)

# selecting tests who are working 
status = dict()
for k in sources:
    status[k] = True
    if 'test_qmrherm_' in k :
        status[k] = False
    if 'test_qmrherm_4.F90' in k:
        status[k] = True

failing_tests = [ k for k in status if status[k] is False ]

# collecting module usage for each test
uses = dict() 
for source in sources:
    modules_used = []
    with open(source,'r') as f:
        lines = f.readlines()
        for line in lines:
            exclmidx = line.find('!')
            if exclmidx is not -1 :
                line = line[:exclmidx]
            if 'use ' in line:
                words = line.split()[:2]
                if words[0] == 'use':
                    sanitized = words[1].replace(',','')
                    modules_used.append(sanitized.lower())
    uses[source] = modules_used                


# checking if a combination of two modules causes issues
def check_interactionsN(N):
    combinations = dict() 
    for k,v in uses.items():
        combinations[k] = list(it.combinations(v,N))
    
    # working pair of modules
    working = set()
    for k,v in combinations.items():
        if status[k] :
            for item in v:
                working.add(item)
                
    # possibly not working pair of modules
    potentially_not_working = set()
    for k,v in combinations.items():
        if not status[k] :
            for item in v:
                potentially_not_working.add(item)
                
    not_working = potentially_not_working.difference(working)
    print("problems with combinations of {:d} modules:".format(N))
    print(not_working)
    
    #possibly not working single modules shares between all tests failing
    potentially_not_working_common = set(combinations[failing_tests[0]])
    for k,v in combinations.items():
        if not status[k] :
            potentially_not_working_common.intersection_update(v)
                
    not_working_common = potentially_not_working_common.difference(working)
    
    print("problems with combinations of {:d} modules (common):".format(N))
    print(not_working_common)

check_interactionsN(1)
check_interactionsN(2)
check_interactionsN(3)
