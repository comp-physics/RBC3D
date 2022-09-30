include ../Makefile.in

# Object files
OBJS = $(patsubst %.F90,%.o,$(wildcard Mod*.F90))

all : libcommon.a

libcommon.a : $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

clean ::
	rm -f *.o *.mod *.a core*

# Dependency
depend .depend : $(wildcard *.F90)
	makedepf90 *.F90 > .depend

include .depend
