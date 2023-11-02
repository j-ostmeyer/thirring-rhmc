import os
import pathlib
import re
import argparse
import glob
import time

def open_file(filename, format):
    current_directory = pathlib.Path(__file__).parent.resolve()
    file_path = os.path.join(current_directory, filename)
    if not os.path.isfile(file_path):
        print(file_path + ": File not found.")
        exit(1)
    return open(file_path, format)

def get_first_integer_from_line(line):
    search = re.search("\d+", line)
    return int(line[search.start():search.end()])

def get_acceptance_rate(output_file):
    total_line = ""
    accepted_line = ""
    for line in output_file:
        if total_line != "":
            accepted_line = line
            break
        if " averages for last" in line:
            total_line = line

    if total_line == "" or accepted_line == "":
        print("Trajectry data not found in " + output_file.name)
        exit(1)

    # Get total number of attempted trajectories
    total_trajectories = get_first_integer_from_line(total_line)

    # Get total accepted trajectories 
    accepted_trajectories = get_first_integer_from_line(accepted_line)

    print("Acceptance rate: " + str(accepted_trajectories / total_trajectories))

def get_exp_dH(output_file):
    for line in output_file:
        if "<exp-dH>=" in line:
            break
    print("exp-dH: " + line.split("<exp-dH>=")[-1].strip())

def get_runtime(o_file):
    first_line = o_file.readline().decode()

    try:  # catch OSError in case of a one line file 
        o_file.seek(-2, os.SEEK_END)
        while o_file.read(1) != b'\n':
            o_file.seek(-2, os.SEEK_CUR)
    except OSError:
        o_file.seek(0)

    last_line = o_file.readline().decode()

    start = time.strptime(first_line, '%a %b %d %H:%M:%S %Z %Y')
    end = time.strptime(last_line, '%a %b %d %H:%M:%S %Z %Y')

    print(start - end)

def main():
    parser = argparse.ArgumentParser(description='Extract output from TEST_OUTPUT_* dir')
    parser.add_argument('output_dir', type=str, nargs=1,
                        help='The relative path to the output directory to extract from.')

    args = parser.parse_args()

    output_file = open_file(os.path.join(args.output_dir[0], "output"), 'r')
    o_filename = glob.glob(args.output_dir[0] + '/ThirringTest.o*')[0].split("/")[-1]
    o_file = open_file(os.path.join(args.output_dir[0], o_filename), 'rb')

    get_runtime(o_file)
    get_acceptance_rate(output_file)
    get_exp_dH(output_file)

    output_file.close()
    o_file.close()

if __name__ == "__main__":
    main()
