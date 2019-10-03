#!/usr/bin/env python3
import pandas as pd

table = pd.read_table(filepath_or_buffer='bookkeeping.txt',comment='#',header=0,sep=r'\s+')
table = table.set_index(['SLURM_JOB_ID','MeasType'])

for filename in ['fort.200','fort.100','fort.11','control','output']:
    with open(filename,'r') as f_input:
        f_input_lines = f_input.readlines()
        for slurm_job_id in table.index.levels[0]:
            end = table.loc[(slurm_job_id,'end'),filename]
            start = table.loc[(slurm_job_id,'start'),filename]
            if start != end:
                f_out_name = filename + '.' + str(slurm_job_id)
                with open(f_out_name,'w') as f_output:
                    for line in f_input_lines[start:end]:
                        f_output.write(line)

            
    
