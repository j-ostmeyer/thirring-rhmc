# Tests

This directory contains several unit tests to aid in the future development of this code base. Unit tests are an excellent tool to validate any new features or changes made to the code.

The principal is that, a single unit (subroutine, function, etc.) runs in isolation with well-defined inputs. The output of these runs are then validated against well-defined, expected values. 

For example, 
- If we wish to improve the performance of the subroutine `congrad` in [measure_module.F90](../measure_module.F90), whilst ensuring we do not break the currently implemented functionality or lose any accuracy, we can utilise the test [test_congrad.F90](./test_congrad.F90). 
- This test sets up some well-defined inputs (`Phi`, `u`, `X`, etc.) and then calls the subroutine `congrad` which makes changes to `X`. 
- The test then compares the newly calculated value of `X` with a "snapshot" of `X` that is known to be correct. This is stored in the data file [test_congrad.dat](./test_congrad.dat).
- If the two versions of `X` match, the test passes and you should see an output similar to...
<img src="congrad-test-output.png"/>
- We can see from the above that `test_congrad.F90` validates `X` as well as the value of `itercg`.
- So long as this test passes it is likely that any changes we have made to `congrad` in `measure_module.F90` have not broken the subroutine. 

It is important to note that a single unit test is rarely enough to be certain we haven't broken anything. Therefore, it is best to write and run as many unit tests as possible to keep increasing your confidence. Additionally, [regression tests](./e2e_tests/README.md) can give us further confidence by testing the interaction between subroutines.

## Running the Unit Tests

A shell script is provided to allow easy running of the unit tests, [run_tests.sh](./run_tests.sh). This script...
1. **Takes in optional paramaters**
    - `KSIZE`:        default = 12
    - `KSIZET`:       default = 12
    - `KTHIRD`:       default = 24
    - `NP_X`:         default = 1
    - `NP_Y`:         default = 1
    - `NP_T`:         default = 1
    - `NP_THIRD`:     default = 1
    - `TESTS`:        default = $ALL_TESTS
    - `GENERATE`:     default = 0
    - `SKIP_COMPILE`: default = 0
2. **Compiles** the test executables using the provided parameters, if `SKIP_COMPILE != 1`.
3. **Runs the test**
    - If `GENERATE = 1`, new data files are generated for all tests specified through the parameter `TESTS`
    - IF `GENERATE != 1`, the tests specified though `TESTS` are ran.

## Generating new data files

Sometimes, if a significate change to a subroutines implementation or purpose has been made. it makes sense to alter the data that is known to be correct (for example, [test_congrad.dat](./test_congrad.dat)), however, this should not be done lightly. 

To generate new data files, you simply run the test you wish to generate with the `GENERATE` flag switched on. For example, for `test_congrad.F90`,
```sh
GENERATE=1 TESTS=test_congrad ./run_tests.sh
```

If you are happy with the generated data file(s) they should then be committed to git and used for all future testing. 
