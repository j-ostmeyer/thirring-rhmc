#!/usr/bin/env python3
import benchmarks_lib as bl
import os
from glob import glob
import pandas as pd 
from tabulate import tabulate

if __name__ == '__main__':
    try:
        setup = pd.read_csv(argv[1], sep = '\s+')
        glob_expr = argv[2]
    except:
        print(f'Usage: {argv[0]} setup_file glob_expr')
        print(f"Example: {argv[0]} setup_file 'slurm*'")
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

        timing_dict = dict(zip(get_benchmark_run_order(),
                               parse_timing_file(timing_file)))

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
                    tablafmt = 'psql'))



        




