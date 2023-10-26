#!/bin/bash

SAMPLES_DIR="./samples"
CONS_DIR="${SAMPLES_DIR}/cons"
MIDOUT_file="${SAMPLES_DIR}/midout"
MYRIAD_SUBMIT_FILE="myriad_submit.sh"

# get inputs from cli
if [[ ( -z "${KSIZE}" )  || ( -z "${KSIZET}" )   || \
	  ( -z "${KTHIRD}" ) || ( -z "${ITER2}" )    || \
	  ( -z "${NP_X}" )   || ( -z "${NP_Y}" )     || \
	  ( -z "${NP_T}" )   || ( -z "${NP_THIRD}" ) || ]]; then
	echo "Please ensure, the env vars KSIZE, KSIZET, KTHIRD and ITER2 are define"
	exit 0
else
	echo "Using config:"
	echo "  KSIZE:  ${KSIZE}"
	echo "  KSIZET: ${KSIZET}"
	echo "  KTHIRD: ${KTHIRD}"
	echo "  NP_X:  ${NP_X}"
	echo "  NP_Y:  ${NP_Y}"
	echo "  NP_T:  ${NP_T}"
	echo "  NP_THIRD:  ${NP_THIRD}"
fi

# Check we have a con file matching the inputs
CON_FILE="${CONS_DIR}/con_${KSIZE}x${KSIZET}"
if [ -f "${CON_FILE}" ]; then
	echo "Found matching con file: ${CON_FILE}"
else
	echo "No matching con file found for user inputs"
	exit 0
fi

# Create output folder TEST_OUTPUT_<KSIZE>_<KSIZET>_<KTHIRD>_<ITER2>
OUTPUT_DIR="TEST_OUTPUT_${KSIZE}_${KSIZET}_${KTHIRD}_${IETR2}"
mkdir $OUTPUT_DIR

# replace iter2 in midout with user input
sed "s/<ITER2>/${ITER2}/g" "${OUTPUT_DIR}/midout"

# update params.F90 with user inputted KSIZE* values
sed "s/#define KSIZE .*/#define KSIZE ${KSIZE}/g" params.F90
sed "s/#define KSIZET .*/#define KSIZET ${KSIZET}/g" params.F90
sed "s/#define KTHIRD .*/#define KTHIRD ${KTHIRD}/g" params.F90

# Update MkFlags with NP values 
sed "s/NP_X=.*/NP_X=${NP_X}/g" MkFlags
sed "s/NP_Y=.*/NP_Y=${NP_Y}/g" MkFlags
sed "s/NP_T=.*/NP_T=${NP_T}/g" MkFlags
sed "s/NP_THIRD=.*/NP_THIRD=${NP_THIRD}/g" MkFlags

# run make
make

# Update submission script with total NP
NP_TOTAL="$(($NP_X + $NP_Y + $NP_T + $NP_THIRD))"
echo "Total number of processors: ${NP_TOTAL}"
sed "s/#$ -pe mpi .*/#$ -pe mpi ${NP_TOTAL}/g" myriad_submit.sh

# Copy myriad_submit.sh, con and midout files from samples into output folder
cp $MYRIAD_SUBMIT_FILE "${OUTPUT_DIR}/${MYRIAD_SUBMIT_FILE}"
cp $CON_FILE "${OUTPUT_DIR}/con"
cp $MIDOUT_FILE "${OUTPUT_DIR}/midout"

# submit job
cd $OUTPUT_DIR
qsub myriad_submit.sh
