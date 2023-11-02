#!/bin/bash

SAMPLES_DIR="samples"
CONS_DIR="${SAMPLES_DIR}/cons"
MIDOUT_FILE="${SAMPLES_DIR}/midout"
MYRIAD_SUBMIT_FILE="myriad_submit.sh"

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

# Check we have a con file matching the inputs
CON_FILE="${CONS_DIR}/con_${KSIZE}x${KSIZE}x${KSIZET}"
if [ -f "${CON_FILE}" ]; then
	echo "Found matching con file: ${CON_FILE}"
else
	echo "No matching con file found for user inputs"
	exit 0
fi

# Create output folder TEST_OUTPUT_<KSIZE>_<KSIZET>_<KTHIRD>_<ITER2>
OUTPUT_DIR="TEST_OUTPUT_${KSIZE}_${KSIZET}_${KTHIRD}_${ITER2}_${NP_X}_${NP_Y}_${NP_T}_${NP_THIRD}"
if [ -d "${OUTPUT_DIR}" ]; then
	echo "Error: A submission for these inputs has already been made"
	exit 0
fi
mkdir $OUTPUT_DIR

touch "${OUTPUT_DIR}/inputs"
echo "Using config:" >> "${OUTPUT_DIR}/inputs"
echo "  KSIZE:    ${KSIZE}" >> "${OUTPUT_DIR}/inputs"
echo "  KSIZET:   ${KSIZET}" >> "${OUTPUT_DIR}/inputs"
echo "  KTHIRD:   ${KTHIRD}" >> "${OUTPUT_DIR}/inputs"
echo "  ITER2:    ${ITER2}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_X:     ${NP_X}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_Y:     ${NP_Y}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_T:     ${NP_T}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_THIRD: ${NP_THIRD}" >> "${OUTPUT_DIR}/inputs"

# Copy myriad_submit.sh, con, midout and remez* files from samples into output folder
cp $MYRIAD_SUBMIT_FILE "./${OUTPUT_DIR}/${MYRIAD_SUBMIT_FILE}"
cp $CON_FILE "./${OUTPUT_DIR}/con"
cp $MIDOUT_FILE "./${OUTPUT_DIR}/midout"
cp remez2 remez2g remez4 remez4g "${OUTPUT_DIR}/"

# update ../../params.F90 with user inputted KSIZE* values and setting to start from con 
sed -i "s/#define KSIZE .*/#define KSIZE ${KSIZE}/g" ../../params.F90
sed -i "s/#define KSIZET .*/#define KSIZET ${KSIZET}/g" ../../params.F90
sed -i "s/#define KTHIRD .*/#define KTHIRD ${KTHIRD}/g" ../../params.F90
sed -i "s/integer, parameter :: istart = .*/integer, parameter :: istart = -1/g" ../../params.F90
sed -i "s/integer, parameter :: iread = .*/integer, parameter :: iread = 1/g" ../../params.F90

# Update ../../MkFlags with NP values 
sed -i "s/NP_X=.*/NP_X=${NP_X}/g" ../../MkFlags
sed -i "s/NP_Y=.*/NP_Y=${NP_Y}/g" ../../MkFlags
sed -i "s/NP_T=.*/NP_T=${NP_T}/g" ../../MkFlags
sed -i "s/NP_THIRD=.*/NP_THIRD=${NP_THIRD}/g" ../../MkFlags
sed -i "s/SITE_RANDOM=no/SITE_RANDOM=yes/g" ../../MkFlags

# run make
make -f ../../Makefile

# cp executable into output dir to allow recompiling for a different build without effecting this one
cp ../../bulk_rhmc "${OUTPUT_DIR}/bulk_rhmc"

# replace iter2 in midout with user input
sed -i "s/<ITER2>/${ITER2}/g" "${OUTPUT_DIR}/midout"

# Update submission script with total NP
NP_TOTAL="$(($NP_X * $NP_Y * $NP_T * $NP_THIRD))"
echo "NP_TOTAL: ${NP_TOTAL}"
echo "Total number of processors: ${NP_TOTAL}"
sed -i "s/#$ -pe mpi .*/#$ -pe mpi ${NP_TOTAL}/g" "${OUTPUT_DIR}/myriad_submit.sh"
sed -i "s/<EXECUTABLE_PATH>/${OUTPUT_DIR}\/bulk_rhmc/g" "${OUTPUT_DIR}/myriad_submit.sh"

# submit job
cd $OUTPUT_DIR
qsub myriad_submit.sh
