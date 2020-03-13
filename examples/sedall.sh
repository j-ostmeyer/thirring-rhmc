sed -i 's/use mpi_f08/use mpi/g' *.F90 */*.F90
sed -i 's/use pmpi_f08/use pmpi/g' */*.F90
sed -i 's/type(MPI_Datatype)/integer/g' *.F90 */*.F90
sed -i 's/type(MPI_Comm)/integer/g' *.F90 */*.F90
sed -i 's/type(MPI_Request)/integer/g' *.F90 */*.F90
sed -i 's/type(MPI_request)/integer/g' *.F90 */*.F90
sed -i 's/type(MPI_Status)/integer/g' *.F90 */*.F90
sed -i 's/type(MPI_File)/integer/g' *.F90 */*.F90
