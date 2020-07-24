#!/usr/bin/env python3
import benchmarks_lib as bl
from sys import argv, exit
import os
from glob import glob
import pandas as pd 
from tabulate import tabulate

if __name__ == '__main__':
    try:
        benchmarks_lib_location = argv[1]   
        glob_expr = argv[2]
        setup = pd.concat(pd.read_csv(filename, sep = '\s+') 
                          for filename in argv[3:])
    except Exception as e:
        print(e)
        print(f'Usage: {argv[0]} benchmark_lib_dir glob_expr [setup_files]')
        print(f"Example: {argv[0]} ../benchmarks 'slurm*' config1.tsv")
        exit(1)

    timing_data = []

    for row in setup.itertuples():
        newdir = bl.dirname(row)

        timing_file_glob_expr = os.path.join(newdir,glob_expr)
        timing_file_candidates = glob(timing_file_glob_expr)
        assert len(timing_file_candidates) == 1, \
                f"There are more than 1 file matching"\
                f"{timing_file_glob_expr}\n"\
                f"{timing_file_candidates}"

        timing_file = timing_file_candidates[0]   

        timing_dict = dict(zip(bl.get_benchmark_run_order(benchmarks_lib_location),
                               bl.parse_timing_file(timing_file)))
        timing_datum = dict(**dict(row._asdict()),
                            **timing_dict)

        timing_data.append(timing_datum)

    timing_df = pd.DataFrame(data = timing_data)

    timing_df.to_csv("timing.tsv",sep = '\t')
    
    with open("timing_table.txt",'w') as f:
        f.write(
                tabulate(
                    timing_df,
                    headers = 'keys',
                    showindex = False,
                    tablefmt = 'psql'))



        




