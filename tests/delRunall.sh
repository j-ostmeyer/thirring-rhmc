source ./MkFlags

MPI_RUNNER="mpirun --oversubscribe -n $((${NP_X} * ${NP_Y} * ${NP_T} * ${NP_THIRD}))"

TESTS=( test_dslash_1 test_dslash_3 test_dslash_5 \
        test_dslashd_1 test_dslashd_3 test_dslashd_5 \
        test_derivs test_qmrherm_0 test_qmrherm_1 test_qmrherm_2 test_qmrherm_3 \
        test_qmrherm_4 test_gauss0 test_gaussp test_measure test_hamilton \
        test_force test_force_regression test_congrad test_congrad_1 )

for test in "${TESTS[@]}"; do
  echo "Testing " $test
  $MPI_RUNNER $test
  echo ""
done
