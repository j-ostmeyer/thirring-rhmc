#!/bin/bash

SAMPLES_DIR="../samples"
CONS_DIR="${SAMPLES_DIR}/cons"
MIDOUT_FILE="${SAMPLES_DIR}/midout"
MYRIAD_SUBMIT_FILE="myriad_submit_test.sh"

# A comma separated list of all tests to be ran
TESTS="test_dslash_1"

# get inputs from cli
if [[ ( -z "${KSIZE}" ) || ( -z "${KSIZET}" ) || \
	  ( -z "${KTHIRD}" ) || ( -z "${ITER2}" ) || \
	  ( -z "${NP_X}" ) || ( -z "${NP_Y}" ) || \
	  ( -z "${NP_T}" ) || ( -z "${NP_THIRD}" ) ]]; then
	echo "Please ensure, the env vars KSIZE, KSIZET, KTHIRD, ITER2, NP_X, NP_Y, NP_T and NP_THIRD are defined"
	exit 0
else
	echo "Using config:"
	echo "  KSIZE:    ${KSIZE}"
	echo "  KSIZET:   ${KSIZET}"
	echo "  KTHIRD:   ${KTHIRD}"
	echo "  ITER2:    ${ITER2}"
	echo "  NP_X:     ${NP_X}"
	echo "  NP_Y:     ${NP_Y}"
	echo "  NP_T:     ${NP_T}"
	echo "  NP_THIRD: ${NP_THIRD}"
fi

# TODO: Write above to file

# Create output folder TEST_OUTPUT_<KSIZE>_<KSIZET>_<KTHIRD>_<ITER2>
OUTPUT_DIR="TEST_OUTPUT_${KSIZE}_${KSIZET}_${KTHIRD}_${ITER2}_${NP_X}_${NP_Y}_${NP_T}_${NP_THIRD}"
if [ -d "${OUTPUT_DIR}" ]; then
	echo "Error: A submission for these inputs has already been made"
	exit 0
fi
mkdir $OUTPUT_DIR

# Copy MYRIAD_SUBMIT_FILE into output folder
cp $MYRIAD_SUBMIT_FILE "./${OUTPUT_DIR}/${MYRIAD_SUBMIT_FILE}"
cp $CON_FILE "./${OUTPUT_DIR}/con"
cp $MIDOUT_FILE "./${OUTPUT_DIR}/midout"
cp remez2 remez2g remez4 remez4g "${OUTPUT_DIR}/"

# update params.F90 with user inputted KSIZE* values and setting to start from con 
sed -i "s/#define KSIZE .*/#define KSIZE ${KSIZE}/g" test_params.F90
sed -i "s/#define KSIZET .*/#define KSIZET ${KSIZET}/g" test_params.F90
sed -i "s/#define KTHIRD .*/#define KTHIRD ${KTHIRD}/g" test_params.F90
sed -i "s/integer, parameter :: istart = .*/integer, parameter :: istart = -1/g" test_params.F90
sed -i "s/integer, parameter :: iread = .*/integer, parameter :: iread = 1/g" test_params.F90
# Update MkFlags with NP values 
sed -i "s/NP_X=.*/NP_X=${NP_X}/g" MkFlags
sed -i "s/NP_Y=.*/NP_Y=${NP_Y}/g" MkFlags
sed -i "s/NP_T=.*/NP_T=${NP_T}/g" MkFlags
sed -i "s/NP_THIRD=.*/NP_THIRD=${NP_THIRD}/g" MkFlags
sed -i "s/SITE_RANDOM=no/SITE_RANDOM=yes/g" MkFlags

# Compile
make -f ./MakefileNew

# cp executable into output dir to allow recompiling for a different build without effecting this one
cp "${TESTS}" "${OUTPUT_DIR}/"

# replace iter2 in midout with user input
sed -i "s/<ITER2>/${ITER2}/g" "${OUTPUT_DIR}/midout"

# Run test
NP_TOTAL="$(($NP_X * $NP_Y * $NP_T * $NP_THIRD))"
cd $OUTPUT_DIR
for TEST in $TESTS; do
	mpirun --oversubscribe -n $NP_TOTAL "./${TEST}"
done

# # Update submission script with total NP
# NP_TOTAL="$(($NP_X * $NP_Y * $NP_T * $NP_THIRD))"
# echo "NP_TOTAL: ${NP_TOTAL}"
# echo "Total number of processors: ${NP_TOTAL}"
# sed -i "s/#$ -pe mpi .*/#$ -pe mpi ${NP_TOTAL}/g" "${OUTPUT_DIR}/${MYRIAD_SUBMIT_FILE}"

# sed -i "s/<EXECUTABLE_PATH>/${OUTPUT_DIR}\/${TEST}/g" "${OUTPUT_DIR}/${MYRIAD_SUBMIT_FILE}"

# # submit job
# cd $OUTPUT_DIR
# qsub $MYRIAD_SUBMIT_FILE
