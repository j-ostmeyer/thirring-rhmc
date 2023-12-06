#!/bin/bash

# A space separated list of all tests to be ran
ALL_TESTS="test_derivs test_meson test_measure test_congrad test_dslash_shamir test_dslashd_shamir test_dslashd_reqs_shamir test_dslash_wilson test_dslashd_wilson test_dslashd_reqs_wilson"

# Read optional user inputs
KSIZE="${KSIZE:-12}"
KSIZET="${KSIZET:-12}"
KTHIRD="${KTHIRD:-24}"
NP_X="${NP_X:-1}"
NP_Y="${NP_Y:-1}"
NP_T="${NP_T:-1}"
NP_THIRD="${NP_THIRD:-1}"
TESTS="${TESTS:-$ALL_TESTS}"
GENERATE="${GENERATE:-0}"
SKIP_COMPILE="${SKIP_COMPILE:-0}"

echo "Using config:"
echo "  KSIZE:        ${KSIZE}"
echo "  KSIZET:       ${KSIZET}"
echo "  KTHIRD:       ${KTHIRD}"
echo "  NP_X:         ${NP_X}"
echo "  NP_Y:         ${NP_Y}"
echo "  NP_T:         ${NP_T}"
echo "  NP_THIRD:     ${NP_THIRD}"
echo "  TESTS:        ${TESTS}"
echo "  GENERATE:     ${GENERATE}"
echo "  SKIP_COMPILE: ${SKIP_COMPILE}"

# Turn on generate if set by user
if [ $GENERATE -ne 0 ]; then
	for TEST in ${TESTS//,/ }
	do
		sed -i "s/logical :: generate = .false./logical :: generate = .true./g" "test_utils.F90"
	done
else 
    for TEST in ${TESTS//,/ }
	do
		sed -i "s/logical :: generate = .true./logical :: generate = .false./g" "test_utils.F90"
	done
fi

if [ $SKIP_COMPILE -ne 1 ]; then
	# update params.F90 with user inputted KSIZE* values and setting to start from con 
	sed -i "s/#define KSIZE .*/#define KSIZE ${KSIZE}/g" test_params.F90
	sed -i "s/#define KSIZET .*/#define KSIZET ${KSIZET}/g" test_params.F90
	sed -i "s/#define KTHIRD .*/#define KTHIRD ${KTHIRD}/g" test_params.F90
	# Not clear if this is needed but will result in reading in a con file
	# sed -i "s/integer, parameter :: istart = .*/integer, parameter :: istart = -1/g" test_params.F90
	# sed -i "s/integer, parameter :: iread = .*/integer, parameter :: iread = 1/g" test_params.F90

	# Update MkFlags with NP values 
	# TODO: only edit if desired value is not already present
	sed -i "s/NP_X=.*/NP_X=${NP_X}/g" MkFlags
	sed -i "s/NP_Y=.*/NP_Y=${NP_Y}/g" MkFlags
	sed -i "s/NP_T=.*/NP_T=${NP_T}/g" MkFlags
	sed -i "s/NP_THIRD=.*/NP_THIRD=${NP_THIRD}/g" MkFlags
	sed -i "s/SITE_RANDOM=no/SITE_RANDOM=yes/g" MkFlags

	# Compile
	make "TESTS=${TESTS}"
fi

# Run test
NP_TOTAL="$(($NP_X * $NP_Y * $NP_T * $NP_THIRD))"
for TEST in ${TESTS//,/ }
do
	if [ -f "${TEST}" ]; then
		echo "Test ${TEST}:"
		mpirun -n $NP_TOTAL "./${TEST}"
	fi
done
