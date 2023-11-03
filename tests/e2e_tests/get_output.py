import os
import argparse
import glob
from datetime import datetime
from utils import get_exp_dH, get_acceptance_rate, open_file

def print_runtime(o_file):
    first_line = o_file.readline().decode()

    try:  
        o_file.seek(-2, os.SEEK_END)
        while o_file.read(1) != b'\n':
            o_file.seek(-2, os.SEEK_CUR)
    except OSError: # catch OSError in case of a one line file 
        o_file.seek(0)

    last_line = o_file.readline().decode()

    try:
        format = '%a %d %b %H:%M:%S %Z %Y'
        start = datetime.strptime(first_line.strip(), format)
        end = datetime.strptime(last_line.strip(), format)
        print("Runtime: " + str(end - start))
    except ValueError as e:
        print("Runtime not found in " + o_file.name)

def main():
    parser = argparse.ArgumentParser(description='Extract output from TEST_OUTPUT_* dir')
    parser.add_argument('output_dir', type=str, nargs=1,
                        help='The relative path to the output directory to extract from.')

    args = parser.parse_args()

    output_file = open_file(os.path.join(args.output_dir[0], "output"), 'r')
    o_filename = glob.glob(args.output_dir[0] + '/ThirringTest.o*')[0].split("/")[-1]
    o_file = open_file(os.path.join(args.output_dir[0], o_filename), 'rb')

    print_runtime(o_file)
    print("Acceptance rate: " + str(get_acceptance_rate(output_file)))
    exp_dH = get_exp_dH(output_file)
    print("exp-dH: " + str(exp_dH[0]) + " +/- " +  str(exp_dH[1]))

    output_file.close()
    o_file.close()

if __name__ == "__main__":
    main()
