#!/bin/bash

SOURCE_DIR="../.."
SAMPLES_DIR="samples"
CONS_DIR="${SAMPLES_DIR}/cons"
MIDOUT_FILE="${SAMPLES_DIR}/midout"
SUBMIT_FILE="csd_submit.sh"

# E2E test inputs
KSIZE="4"
KSIZET="4"
KTHIRD="8"
ITER2="100"
NP_THIRD="1"
USE_MEASURE="1"
NP_X="${NP_X:-4}"
NP_Y="${NP_Y:-4}"
NP_T="${NP_T:-2}"

# Check we have a con file matching the inputs
CON_FILE="${CONS_DIR}/con_${KSIZE}x${KSIZE}x${KSIZET}"
if [ -f "${CON_FILE}" ]; then
	echo "Found matching con file: ${CON_FILE}"
else
	echo "No matching con file found for user inputs"
	exit 0
fi

# Create output folder
OUTPUT_DIR="TEST_OUTPUT_${KSIZE}_${KSIZET}_${KTHIRD}_${ITER2}_${NP_X}_${NP_Y}_${NP_T}_${NP_THIRD}_${USE_MEASURE}"
if [ -d "${OUTPUT_DIR}" ]; then
	echo "Error: A submission for these inputs has already been made"
	exit 0
fi
mkdir $OUTPUT_DIR

touch "${OUTPUT_DIR}/inputs"
echo "Using config:" >> "${OUTPUT_DIR}/inputs"
echo "  KSIZE:       ${KSIZE}" >> "${OUTPUT_DIR}/inputs"
echo "  KSIZET:      ${KSIZET}" >> "${OUTPUT_DIR}/inputs"
echo "  KTHIRD:      ${KTHIRD}" >> "${OUTPUT_DIR}/inputs"
echo "  ITER2:       ${ITER2}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_X:        ${NP_X}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_Y:        ${NP_Y}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_T:        ${NP_T}" >> "${OUTPUT_DIR}/inputs"
echo "  NP_THIRD:    ${NP_THIRD}" >> "${OUTPUT_DIR}/inputs"
echo "  USE_MEASURE: ${USE_MEASURE}" >> "${OUTPUT_DIR}/inputs"

# Copy csd_submit.sh, con, midout and remez* files from samples into output folder
cp $SUBMIT_FILE "./${OUTPUT_DIR}/${SUBMIT_FILE}"
cp $CON_FILE "./${OUTPUT_DIR}/con"
cp $MIDOUT_FILE "./${OUTPUT_DIR}/midout"
cp "${SOURCE_DIR}/remez2" "${SOURCE_DIR}/remez2g" "${SOURCE_DIR}/remez4" "${SOURCE_DIR}/remez4g" "${OUTPUT_DIR}/"

# update params.F90 with user inputted KSIZE* values and setting to start from con 
sed -i "s/#define KSIZE .*/#define KSIZE ${KSIZE}/g" "${SOURCE_DIR}/params.F90"
sed -i "s/#define KSIZET .*/#define KSIZET ${KSIZET}/g" "${SOURCE_DIR}/params.F90"
sed -i "s/#define KTHIRD .*/#define KTHIRD ${KTHIRD}/g" "${SOURCE_DIR}/params.F90"
sed -i "s/integer, parameter :: istart = .*/integer, parameter :: istart = -1/g" "${SOURCE_DIR}/params.F90"
sed -i "s/integer, parameter :: iread = .*/integer, parameter :: iread = 1/g" "${SOURCE_DIR}/params.F90"

# Switch between meson and measure
measure_line='call measure(pbp, respbp, ancgm, am, imass, isweep + isweep_total_start)'
meson_line='call meson(rescgm,itercg,ancgm,am,imass, isweep + isweep_total_start)'
if [ $(($USE_MEASURE)) -eq 1 ]; then
	sed -E -i "s|^(\s*!*)\s*$measure_line|\1$measure_line|" "${SOURCE_DIR}/dwf3d_lib.F90"
	sed -E -i "s|^(\s*)(!*)\s*$meson_line|\1!\2$meson_line|" "${SOURCE_DIR}/dwf3d_lib.F90"
else
	sed -E -i "s|^(\s*!*)\s*$meson_line|\1$meson_line|" "${SOURCE_DIR}/dwf3d_lib.F90"
	sed -E -i "s|^(\s*)(!*)\s*$measure_line|\1!\2$measure_line|" "${SOURCE_DIR}/dwf3d_lib.F90"
fi

# Update MkFlags with NP values 
sed -i "s/NP_X=.*/NP_X=${NP_X}/g" "${SOURCE_DIR}/MkFlags"
sed -i "s/NP_Y=.*/NP_Y=${NP_Y}/g" "${SOURCE_DIR}/MkFlags"
sed -i "s/NP_T=.*/NP_T=${NP_T}/g" "${SOURCE_DIR}/MkFlags"
sed -i "s/NP_THIRD=.*/NP_THIRD=${NP_THIRD}/g" "${SOURCE_DIR}/MkFlags"
sed -i "s/SITE_RANDOM=no/SITE_RANDOM=yes/g" "${SOURCE_DIR}/MkFlags"

# run make
pushd "${SOURCE_DIR}/"
make
popd

# cp executable into output dir to allow recompiling for a different build without effecting this one
cp "${SOURCE_DIR}/bulk_rhmc" "${OUTPUT_DIR}/bulk_rhmc"

# replace iter2 in midout with user input
sed -i "s/<ITER2>/${ITER2}/g" "${OUTPUT_DIR}/midout"

# Update submission script with total NP
NP_TOTAL="$(($NP_X * $NP_Y * $NP_T * $NP_THIRD))"
echo "NP_TOTAL: ${NP_TOTAL}"
echo "Total number of processors: ${NP_TOTAL}"
sed -i "s/#SBATCH --ntasks=.*/#SBATCH --ntasks=${NP_TOTAL}/g" "${OUTPUT_DIR}/${SUBMIT_FILE}"
sed -i "s/<EXECUTABLE_PATH>/${OUTPUT_DIR}\/bulk_rhmc/g" "${OUTPUT_DIR}/${SUBMIT_FILE}"

if [ $(($NP_TOTAL)) -gt 56 ]
then
   	echo "Error: Cannot submit job with more than 56 processors to csd without manually updating nodes."
	exit 0
fi

# submit job
cd $OUTPUT_DIR
sbatch $SUBMIT_FILE
