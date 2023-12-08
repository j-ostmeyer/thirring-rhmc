# 
# This script is intented to act as an end-to-end test for the thirring-rhmc code base.
# The script validates simulation outputs such as the acceptance rate and the <exp-dH>
# These values are extracted from the output files (i.e. output, fort.xxx) of a 
# completed simulation.
# 
# To ensure these tests are reliable, consistent inputs must be used to generate the 
# output files to be tested. Therefore, these inputs should be:
#     -  KSIZE:    4
#     -  KSIZET:   4
#     -  KTHIRD:   8
#     -  ITER2:    100
#     -  NP_X:     4
#     -  NP_Y:     4
#     -  NP_T:     2
#     -  NP_THIRD: 1
# We should also be running using the measure method rather than meson in dwf3d_lib.f90.
from os import path
import argparse
from utils import get_exp_dH, get_acceptance_rate, open_file

def print_results(passed, failed):
    colour_end = "\033[0m"
    green_start = "\033[92m"
    red_start =  '\033[91m'
    if len(passed) > 0:
        print("---------------------------------------------")
        for test in passed:
            print(green_start + "Passed: " + test["function"] + " - " + test["name"] + colour_end)
    if len(failed) > 0:
        print("---------------------------------------------")
        for test in failed:
            print(red_start + "Failed: " + test["function"] + " - " + test["name"] + colour_end + ", " + test["message"])
    print("---------------------------------------------")
    print(str(len(passed)) + " of " + str(len(passed) + len(failed)) + " tests have passed")

def run_test(file, test_function):
    passed = []
    failed = []
    test_dir = file.name.split("/")[-2]
    try:
        test_function(file)
        passed.append({"name": test_dir, 
                       "function": test_function.__name__})
    except AssertionError as e:
        failed.append({"name": test_dir,
                       "function": test_function.__name__,
                       "message": str(e)})
    return passed, failed

def test_acceptance(output_file):
    # Acceptance rate
    expected_acceptance = 0.8
    actual_acceptance = get_acceptance_rate(output_file)
    print("Acceptance rate: " + str(actual_acceptance))
    assert actual_acceptance >= expected_acceptance, "Acceptance " + str(actual_acceptance) + " is < " + str(int(expected_acceptance * 100)) + "%"

def test_exp_dH(output_file):
    # exp-dH
    expected_exp_dH = 1
    actual_exp_dH_val, actual_exp_dH_err = get_exp_dH(output_file)
    print("exp-dH: " + str(actual_exp_dH_val) + " +/- " +  str(actual_exp_dH_err))
    assert abs(expected_exp_dH - actual_exp_dH_val) <= actual_exp_dH_err, "Expected value of " + str(expected_exp_dH) + " is outside exp-dH = " + str(actual_exp_dH_val) + " +/- " + str(actual_exp_dH_err)

def main():
    parser = argparse.ArgumentParser(description='Extract output from TEST_OUTPUT_* dir')
    parser.add_argument('output_dirs', type=str, nargs='+',
                        help='The relative paths to the output directories to test.')
    
    args = parser.parse_args()

    passed = []
    failed = []
    for output_dir in args.output_dirs:
        output_file = open_file(path.join(output_dir, "output"))
        
        acceptance_passed, acceptance_failed = run_test(output_file, test_acceptance)
        exp_passed, exp_failed = run_test(output_file, test_exp_dH)
        
        output_file.close()

        passed += acceptance_passed + exp_passed
        failed += acceptance_failed + exp_failed

    print_results(passed, failed)

if __name__ == "__main__":
    main()
