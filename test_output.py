import sys
import os
import pathlib
import re

output_option = "--output-files"

def help(executable):
    print("Usage: " + executable + " " + output_option + " <filename_list>")
    print("   " + output_option + " <filename_list>: A space separated list of relative patha to output files to be tested.")
    exit()

def print_results(passed, failed):
    if len(passed) > 0:
        print("---------------------------------------------")
        for test in passed:
            print("Passed: " + test["name"])
    if len(failed) > 0:
        print("---------------------------------------------")
        for test in failed:
            print("Failed: " + test["name"] + ", " + test["message"])
    print("---------------------------------------------")
    print(str(len(passed)) + " of " + str(len(passed) + len(failed)) + " tests have passed")

def get_lines_from_file(output_filename):
    current_directory = pathlib.Path(__file__).parent.resolve()
    output_file_path = os.path.join(current_directory, output_filename)
    if not os.path.isfile(output_filename):
        print(output_file_path + ": File not found.")
        exit(1)

    # Test its contents 
    # Find appropriate lines in output file
    output_file = open(output_file_path)
    total_line = ""
    accepted_line = ""
    for line in output_file:
        if total_line != "":
            accepted_line = line
            break
        if " averages for last" in line:
            total_line = line
    
    if total_line == "" or accepted_line == "":
        print("Trajectry data not found in " + output_file_path)
        exit(1)

    return total_line, accepted_line

def get_first_integer_from_line(line):
    search = re.search("\d+", line)
    return int(line[search.start():search.end()])

def test_acceptance(output_filename):
    total_line, accepted_line = get_lines_from_file(output_filename)

    # Get total number of attempted trajectories
    total_trajectories = get_first_integer_from_line(total_line)

    # Get total accepted trajectories 
    accepted_trajectories = get_first_integer_from_line(accepted_line)

    # Check that acceptance is 80% or above 
    assert (accepted_trajectories / total_trajectories) >= 0.8, "Acceptance is < 80%"

def main():
    args = sys.argv[1:]
    if len(args) < 2 or args[0] != output_option:
        help(sys.argv[0])

    # Read output file defined via args
    output_filenames = args[1:]

    passed = []
    failed = []
    for filename in output_filenames:
        try:
            test_acceptance(filename)
            passed.append({"name": filename})
        except AssertionError as e:
            failed.append({"name": filename,
                           "message": str(e)})
    
    print_results(passed, failed)

if __name__ == "__main__":
    main()
