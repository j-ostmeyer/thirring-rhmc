import argparse
from utils import get_exp_dH, get_acceptance_rate, open_file

output_option = "--output-files"
fort_option = "--fort-11-files"

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

def run_test(filenames, test_function):
    passed = []
    failed = []
    for filename in filenames:
        try:
            test_function(filename)
            passed.append({"name": filename, 
                           "function": test_function.__name__})
        except AssertionError as e:
            failed.append({"name": filename,
                           "function": test_function.__name__,
                           "message": str(e)})
    return passed, failed

def test_acceptance(output_filename):
    output_file = open_file(output_filename)
    # Acceptance rate
    expected_acceptance = 0.8
    actual_acceptance = get_acceptance_rate(output_file)
    output_file.close()
    assert actual_acceptance >= expected_acceptance, "Acceptance is < " + str(int(expected_acceptance * 10)) + "%"

def test_exp_dH(output_filename):
    output_file = open_file(output_filename)
    # exp-dH
    expected_exp_dH = 1
    actual_exp_dH_val, actual_exp_dH_err = get_exp_dH(output_file)
    output_file.close()
    assert abs(expected_exp_dH - actual_exp_dH_val) <= actual_exp_dH_err, "Expected value of " + str(expected_exp_dH) + " is outside exp-dH = " + str(actual_exp_dH_val) + " +/- " + str(actual_exp_dH_err)

def get_fort_11_sequence(filename):
    file = open_file(filename)
    sequence = []
    for line in file:
        sequence.append(line.split(' ')[2])
    file.close()
    return sequence

def test_fort_11_file(fort_filename):    
    # Get expected sequence 
    expected_sequence = get_fort_11_sequence("samples/ref_fort.11")

    # Get actual sequence 
    actual_sequence = get_fort_11_sequence(fort_filename)

    # Compare with expected sequence
    assert len(actual_sequence) == len(expected_sequence), "The length of the actual sequence " + str(len(actual_sequence)) + " does not match that of the expected sequence " + str(len(expected_sequence))
    failed_index = -1
    for i, actual in enumerate(actual_sequence):
        expected = expected_sequence[i]
        if actual != expected:
            failed_index = i
            break
    assert failed_index == -1, "The two sequences do not match for row " + str(failed_index)

def main():
    parser = argparse.ArgumentParser(description='Extract output from TEST_OUTPUT_* dir')
    parser.add_argument(output_option, dest="outputs", type=str, nargs='+', default=[],
                        help='A space separated list of relative paths to the output files to test.')
    parser.add_argument(fort_option, dest="forts", type=str, nargs='+', default=[],
                        help='A space separated list of relative paths to the fort.11 files to test.')
    
    args = parser.parse_args()

    acceptance_passed, acceptance_failed = run_test(args.outputs, test_acceptance)
    exp_passed, exp_failed = run_test(args.outputs, test_exp_dH)
    fort_passed, fort_failed = run_test(args.forts, test_fort_11_file)

    print_results(acceptance_passed + exp_passed + fort_passed,
                  acceptance_failed + exp_failed + fort_failed)

if __name__ == "__main__":
    main()
