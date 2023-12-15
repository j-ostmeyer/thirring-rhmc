# End to End testing

This is intended to facilitate regression testing this code base. 

Currently, the process is made up of submitting a well defined job to a cluster and then validating the outputs using a python script.

To improve these tests, additional output quantities can be extracted and verified by editting the python script [test_suite.py](./test_suite.py).

## E2E process 
1. Build the code and submission scripts:
    - The script `build_test_dir.sh` will update the src code and compile with the parmeters given to it.
    - The paramaters accepted and the default values are given below:
        - KSIZE=4
        - KSIZET=4
        - KTHIRD=8
        - ITER2=100
        - NP_THIRD=1
        - NP_X=4
        - NP_Y=4
        - NP_T=2
    - To change the defaults value for any of the above, pass it as an environemnt variable. For example `NP_X=2 ./build_test_dir.sh`.
    - The script will produce a directory with a name matching `TEST_OUTPUT_<KSIZE>_<KSIZET>_<KTHIRD>_<ITER2>_<NP_X>_<NP_Y>_<NP_T>_<NP_THIRD>` which should contain th following
        - **bulk_rhmc**: The compiled executable.
        - **midout**: The default midout file to be used for the test run (stored in [samples](./samples/)).
        - **con**: The default con file to be used for the run. This is determined by the `KSIE` and `KSIZET` values passed to `build_test_dir.sh`. The sample con files are stored in [samples](./samples/).
        - **remez2, remez2g, remez4 and remez4g**: The remez dependencies required by the executable.
        - **csd_submit.sh**: The submission script required for submitting jobs to CSD *(Work in progress)*.
        - **myriad_submit.sh**: The submission script required for submitting jobs to the UCL clusters Myriad and Kathleen.

2. Submit your job:
     - Using the files generated in step 1, we can now submit a job to one of our chosen clusters (Myriad, Kathleen or CSD).
     - To do this 
        1. change into the new directory produced by `build_test_dir.sh`
        2. Submit the job. 
            - **Myriad or Kathleen**: `qsub myriad_submit.sh`
            - **CSD**: `sbatch csd_submit.sh`

3. Test the outputs:
    - Once the job submitted in step 2 has completed you can then test the output.
    - To test the output, utilise the provided python script `test_suite.py`. This script extracts `acceptance rate` and `exp-dH` values from the `output` file inside every directory passed to it then prints them to the screen and verifies these values against expected ones. 
    - Example usage:
      ```
      python3 test_suite.py TEST_OUTPUT_4_4_8_100_4_4_2_1 TEST_OUTPUT_4_4_8_100_4_2_4_1
      ```
      The above would test the outputs of two jobs, one with default inputs and another with `NP_Y=4` and `NP_y=2`.
