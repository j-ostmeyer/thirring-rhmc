# Tests

This directory contains several unit tests to aid in the development of this code base. Unit tests are an excellent tool to validate any new features or changes made to the code.

The principal is that, a single unit (subroutine, function, etc.) runs in isolation with well-defined inputs. The output of these runs are then validated against well-defined, expected values. 

Some unit tests include hard-coded expected outputs, while others are "black-box" tests that should be run once (before any changes are made) to generate test data, then run again to compare the outputs of changed subroutines to the generated test data.

For example, 
- If we wish to improve the performance of the subroutine `congrad` in [measure_module.F90](../measure_module.F90), whilst ensuring we do not break the currently implemented functionality or lose any accuracy, we can utilise the test [test_congrad.F90](./test_congrad.F90). 
- This test sets up some well-defined inputs (`Phi`, `u`, `X`, etc.) and then calls the subroutine `congrad` which makes changes to `X`. 
- Since this test requires a `.dat` file, we must generate that before making any changes. The test will run the subroutine and save its output. We can then make changes to the code.
- Running the test after making some change then compares a newly calculated value of `X` with the previously generated "snapshot" of `X` stored in the data file [test_congrad.dat](./test_congrad.dat).
- If the two versions of `X` match, the test passes and you should see an output similar to...
    ``` 
    Test test_congrad:
        Initialising MPI with grid (NP_X * NP_Y * NP_T * NP_THIRD)           1 *           1 *           1 *           1
            Passed test_congrad itercg is the expected value
            Passed test_congrad sum of x is within error
            Passed test_congrad max x is within error
    ```
- We can see from the above that `test_congrad.F90` validates `X` as well as the value of `itercg`.
- So long as this test passes it is likely that any changes we have made to `congrad` in `measure_module.F90` have not broken the subroutine. 

It is important to note that a single unit test is rarely enough to be certain we haven't broken anything. Therefore, it is best to write and run as many unit tests as possible to keep increasing your confidence. Additionally, [regression tests](./e2e_tests/README.md) can give us further confidence by testing the entire simulation against previously recorded results.

## Running the Unit Tests

A shell script is provided to allow easy running of the unit tests, [run_tests.sh](./run_tests.sh). This script...
1. **Takes in optional parameters**
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
  
   <u>**WARNING:**</u> The problem size (`KSIZE`, `KSIZET` and `KTHIRD`) must remain the same for unit tests with hard-coded values to pass and be valid. Unit tests that generate their own test data in `.dat` files are unaffected (however the parameters must remain consistent between generation and a later test).
2. **Compiles** the test executables using the provided parameters, if `SKIP_COMPILE != 1`.
3. **Runs the test**
    - If `GENERATE = 1`, new data files are generated for all tests specified through the parameter `TESTS`
    - IF `GENERATE != 1`, the tests specified though `TESTS` are ran.

## Generating new data files

To generate new data files, you simply run the test you wish to generate with `GENERATE=1`. For example, for `test_congrad.F90`,
```sh
GENERATE=1 TESTS=test_congrad ./run_tests.sh
```
