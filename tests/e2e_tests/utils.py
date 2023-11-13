import os
import pathlib
import re

def open_file(filename, format = 'r'):
    current_directory = pathlib.Path(__file__).parent.resolve()
    file_path = os.path.join(current_directory, filename)
    if not os.path.isfile(file_path):
        print(file_path + ": File not found.")
        exit(1)
    return open(file_path, format)

def get_first_integer_from_line(line):
    search = re.search("\d+", line)
    return int(line[search.start():search.end()])

def get_exp_dH(output_file):
    for line in output_file:
        if "<exp-dH>=" in line:
            break
    exp_dH_str = line.split("<exp-dH>=")[-1].strip()
    return [float(i) for i in exp_dH_str.split(" +/- ")]
        
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

    return accepted_trajectories / total_trajectories