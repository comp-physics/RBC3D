Code requires:
MPI environment with Intel compiler (need mpif90 and ifort)
Load torque with this usually
PETSc 3.0.0 lite (details below)
Intel MKL 10.0 (tarball can be found in ~/packages) which is then to be installed using install.sh in home dir
FFTW (best to use module file here)
Lapack, Blas, and Lapack 95 (install via tar in ~/packages)
Sherepack, netcdf which can just copy from ~_3d_cell_mylib apparently

To compile PetSc (need version 3.0.0):
._config_configure.py --with-cc=mpicc --with-fc=mpif90 --with-blas-lib=_home_bryngel2/packages_BLAS-mpi_blas_mpi.a --with-lapack-lib=_home_bryngel2/packages_lapack-3.6.0-mpi_liblapack.a --with-valgrind-dir=_srv_local_data1_work_ostoich1_valgrind

For GOLUB i used —with-x=0

For new Intel compiler we need in the .bash_profile
export LD_PRELOAD=_usr_local_intel_parallel_studio_xe_2018/compilers_and_libraries_linux_mkl_lib_intel64/libmkl_def.so:_usr_local_intel_parallel_studio_xe_2018/compilers_and_libraries_linux_mkl_lib_intel64/libmkl_avx2.so:_usr_local_intel_parallel_studio_xe_2018/compilers_and_libraries_linux_mkl_lib_intel64/libmkl_core.so:_usr_local_intel_parallel_studio_xe_2018/compilers_and_libraries_linux_mkl_lib_intel64/libmkl_intel_thread.so

module load xpacc-mpich2-intel
module load torque

To install FFTW (never worked on xpacc, ended up using module files):
 $ ./configure CC=icc
--enable-threads \
--enable-mpi MPICC=mpicc \
--prefix=_usr_local/fftw-3.3.3-single
 $ make -j8
 $ make install
