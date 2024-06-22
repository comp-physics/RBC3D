include ../../Makefile.in

# replace $(MKL_LIB) with $(BLAS_LIB) $(LAPACK_LIB) if using lapack install
LIB = $(COMMON_LIB) $(SPHPK_LIB) $(FFTW_LIB) $(PETSC_LIB) $(NETCDF_LIB) $(LAPACK_LIB) $(STATIC)
EXECUTABLES = tube initcond

all : $(EXECUTABLES)

lib : $(wildcard $(WORK_DIR)/common/*.F90)
	make -C $(WORK_DIR)/common

$(EXECUTABLES) : % : %.o $(wildcard $(COMMON_DIR)/*.F90)
	make lib
	$(FC) $(LDFLAGS) -o $@ $< $(LIB)

clean :
	-rm -f $(EXECUTABLES) *.o *.mod *.a core

# Dependency
.depend : $(wildcard *.F90)
	$(MAKEDEPEND_BIN) $^ > .depend

include .depend
