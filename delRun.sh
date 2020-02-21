source ./MkFlags

MPI_RUNNER="mpirun --oversubscribe -n $((${NP_X} * ${NP_Y} * ${NP_T} * ${NP_THIRD}))"

$MPI_RUNNER bulk_rhmc
