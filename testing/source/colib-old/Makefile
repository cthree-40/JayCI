# Makefile for building colib.a archive
OBJECTS := blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f colib6.f colib7.f \
		colib8.f colib9.f colib10.f
FC := gfortran
BLASLIBRARY := $(BLAS_LIB)
AR := ar
RANL := ranlib

colib.a: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS)
	$(AR) colib.a *.o $(BLAS_LIB)
	$(RANL) colib.a

%.o:%.f
	$(FC) $(FFLAGS) -o $@ $< $(BLAS_LIB)
