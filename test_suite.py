import sys
import os
import pathlib
import re

output_option = "--output-files"
fort_option = "--fort-files"

def help():
    print("Usage: " + sys.argv[0] + " " + output_option + " <output_filenames> " + fort_option + " <fort.11_filenames>")
    print("   " + output_option + " <output_filenames>: A space separated list of relative paths to output files to be tested.")
    print("   " + fort_option + " <fort.11_filenames: A space separated list of relative paths to fort.11 files to be tested.")
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

def run_test(filenames, test_function):
    passed = []
    failed = []
    for filename in filenames:
        try:
            test_function(filename)
            passed.append({"name": filename})
        except AssertionError as e:
            failed.append({"name": filename,
                           "message": str(e)})
    return passed, failed

def open_file(filename):
    current_directory = pathlib.Path(__file__).parent.resolve()
    file_path = os.path.join(current_directory, filename)
    if not os.path.isfile(filename):
        print(file_path + ": File not found.")
        exit(1)
    return open(file_path)

def get_first_integer_from_line(line):
    search = re.search("\d+", line)
    return int(line[search.start():search.end()])


def test_acceptance(output_filename):
    expected_acceptance = 0.8

    output_file = open_file(output_filename)
    total_line = ""
    accepted_line = ""
    for line in output_file:
        if total_line != "":
            accepted_line = line
            break
        if " averages for last" in line:
            total_line = line

    output_file.close()

    if total_line == "" or accepted_line == "":
        print("Trajectry data not found in " + output_file.name)
        exit(1)

    # Get total number of attempted trajectories
    total_trajectories = get_first_integer_from_line(total_line)

    # Get total accepted trajectories 
    accepted_trajectories = get_first_integer_from_line(accepted_line)

    assert accepted_trajectories / total_trajectories >= expected_acceptance, "Acceptance is < " + str(int(expected_acceptance * 10)) + "%"

def test_fort_file(fort_filename):
    fort_file = open_file(fort_filename)
    expected_sequece = ["1.339680","1.288082","1.324598","1.323031","1.285908","1.275051","1.284303","1.289887",
                        "1.285864","1.321300","1.329907","1.329907","1.345052","1.325124","1.335083","1.373880",
                        "1.319285","1.328140","1.328140","1.325442","1.307628","1.297789","1.297789","1.373745",
                        "1.378049","1.372292","1.376962","1.372832","1.372832","1.368621","1.330136","1.312386",
                        "1.312386","1.315233","1.383984","1.369509","1.287796","1.270133","1.282026","1.316900",
                        "1.332728","1.344033","1.356403","1.366509","1.366509","1.323285","1.319908","1.319908",
                        "1.319908","1.359130","1.364986","1.323512","1.322849","1.338154","1.349150","1.348525",
                        "1.348525","1.273444","1.286146","1.303770","1.313136","1.312661","1.312476","1.320717",
                        "1.292957","1.295259","1.304003","1.268134","1.242731","1.260418","1.258866","1.249406",
                        "1.249406","1.254013","1.292001","1.292001","1.271733","1.277426","1.277426","1.282797",
                        "1.282797","1.281287","1.295108","1.268947","1.238418","1.230964","1.227291","1.229165",
                        "1.246840","1.248483","1.248483","1.269551","1.320455","1.320913","1.321494","1.325113",
                        "1.325113","1.301641","1.301641","1.301641"]
    # Read third column of fort file
    actual_sequence = []
    for line in fort_file:
        actual_sequence.append(line.split(' ')[2])
    fort_file.close()
    # Compare with expected sequence
    assert len(actual_sequence) == len(expected_sequece)
    failed_index = -1
    for i, actual in enumerate(actual_sequence):
        expected = expected_sequece[i]
        if actual != expected:
            failed_index = i
            break
    assert failed_index == -1, "The two sequences do not match for row " + str(failed_index)

def main():
    args = sys.argv[1:]
    run_output_tests = output_option in args
    run_fort_tests = fort_option in args

    if not (run_output_tests or run_fort_tests):
        help()

    # Get list of output files
    output_filenames = []
    if run_output_tests:
        output_filenames = args[args.index(output_option) + 1:]
        if fort_option in output_filenames:
            output_filenames = output_filenames[:output_filenames.index(fort_option)]

    # Get list of fort.11 files
    fort_filenames = []
    if run_fort_tests:
        fort_filenames = args[args.index(fort_option) + 1:]
        if output_option in fort_filenames:
            fort_filenames = fort_filenames[:fort_filenames.index(output_option)]

    if len(fort_filenames) + len(output_filenames) == 0:
        help()

    output_passed, output_failed = run_test(output_filenames, test_acceptance)
    fort_passed, fort_failed = run_test(fort_filenames, test_fort_file)

    print_results(output_passed + fort_passed,
                  output_failed + fort_failed)

if __name__ == "__main__":
    main()
