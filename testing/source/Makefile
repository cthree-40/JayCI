## Makefile for JayCI
#
## By Christopher L. Malbon
## Yarkony Group
## The Johns Hopkins University
## 
# Currently tested with Intel Fortran Compilers ONLY.
# Currently tested on Linux ONLY.
# Need LAPACK and BLAS. Tested for Intel MKL ONLY.
#
#
#######################################
# Build environment
JOBDIR	:= ..

#SOURCE	:= $(JOBDIR)/source2
SOURCE	:= .
COLIBDIR:= $(SOURCE)/colib
UNIXDIR	:= $(SOURCE)/UNIX
BIN	:= $(JOBDIR)/bin

# Set compiler
# We build on edison
ifeq ($(NERSC_HOST),edison)
	COMPLR	:= ftn -i8
	COLIBFC := ftn -i8
	CC	:= cc
	AR	:= ar rv
	RANL	:= ranlib
# Or we are on Hyperion
else
	COMPLR	:= gfortran -fdefault-integer-8
	COLIBFC	:= gfortran -fdefault-integer-8 
	CC	:= gcc
	AR	:= ar rv
	RANL	:= ranlib
endif

PREFIX		:= sh -c

# Set Libraries
# We build on edison
ifeq ($(NERSC_HOST),edison)
	LIBS1	:=
else
	# we are on hyperion or helios
	LIBS1	:= -L/usr/lib64 -llapack -lblas
endif

LIBS2		:= $(SOURCE)/colib.a

ifneq ($(findstring gfortran,$(COMPLR)),)
#GNU compilers
	CPOPS	:= -c -frecord-marker=4 -O3
	CLIBOPS	:= -c -frecord-marker=4 -O3
	#DEBUG	:= -g -fcheck=all -g -fbacktrace 
	CDEF	:= -c -DFLUSH -DINT64 -DEXTNAME
else
#Intel compilers
	CPOPS	:= -auto -c -assume byterecl -O3
#	DEBUG	:= -debug -check all -check bounds -check uninit -gen-interfaces -warn interfaces -g -traceback
	LKOPT	:= -auto -lpthread
	VTUNE	:= -g -O3
endif

# Executables ------------------------------------------------------------------
EXEC1		:= $(BIN)/jayci_exp.x
EXEC2		:= $(BIN)/jayci.x
TESTEX		:= $(BIN)/test.x

# Clean up ---------------------------------------------------------------------
REMOVE	:= rm -f

#######################################
# Objects

COLIBO	:= blaswrapper.o colib1.o colib2.o colib3.o colib4.o colib5.o colib6.o colib7.o \
			colib8.o colib9.o colib10.o
COLIBF	:= blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f colib6.f colib7.f \
			colib8.f colib9.f colib10.f
UNIXO	:= fdate.o falloc.o fwtime.o hostnm.o flushstdout.o fsize.o
UNIXC	:= fdate.c falloc.c fwtime.c hostnm.c flushstdout.c fsize.c
OBJS1		:= common_util.o truncation1.o
OBJS2		:= $(OBJS1) construct.o  action_util.1.o acthv.o david_util.o initialguess.o davidson.o
OBJSTEST	:= $(OBJS2) test_module.o test.o
DRIV1		:= $(OBJS1) jayci_exp.o
DRIV2		:= $(OBJS2) jayci.o

#######################################
# Preprocessor flags

PREPROC	:= -DDEBUGGING
PROFILE	:= -prof-gen -prof-dir$(JOBDIR)/Profile_Info
PROFILE2	:= -prof-use -ipo -prof-dir$(JOBDIR)/Profile_Info

#######################################
# BUILD
driver1:	$(DRIV1) | $(BIN)
	@echo " Building driver1..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(EXEC1) $(DRIV1) $(LIBS1) $(LIBS2)"
	@echo " ------------------------------------------------------------------------"


driver2:$(DRIV2) | $(BIN)
	@echo " Building driver2..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(EXEC2) $(DRIV2) $(LIBS1) $(LIBS2) $(DEBUG) $(VTUNE) "
	@echo " ------------------------------------------------------------------------"

test:	$(OBJSTEST) | $(BIN)
	@echo " Building test..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(TESTEX) $(OBJSTEST) $(LIBS1) $(LIBS2) $(DEBUG) $(PREPROC)"
	@echo " ------------------------------------------------------------------------"

colib:	$(COLIBO) $(UNIXO)
	chmod +x colib.a
	$(AR) colib.a $(COLIBO) $(UNIXO)
	$(RANL) colib.a

%.o:$(SOURCE)/$(UNIXDIR)/%.c
	@echo " Building $< "
	$(CC) -c $(CDEF) -o $@ $<

%.o:$(SOURCE)/$(COLIBDIR)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COLIBFC) -o $@ $< $(CPOPS) $(VTUNE) "

%.o:$(SOURCE)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -c -o $@ $< $(CPOPS) $(VTUNE)"

%.o:$(SOURCE)/%.f90
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -c -cpp -o $@ $< $(CPOPS) $(DEBUG) $(PREPROC) $(VTUNE)"

clean1:
	$(RM) $(DRIV1) ./*__genmod.f90./*.o ./*.mod 
	@echo " Finished cleaning."


clean2:
	$(RM) $(DRIV2) ./*.o ./*.mod ./*__genmod.f90
	@echo " Finished cleaning."

cleancolib:
	$(RM) ./blaswrapper.o ./colib/colib*.o ./colib*.o ./fwtime.o ./hostnm.o ./falloc.o ./fdate.o ./fsize.o ./flushstdout.o
	@echo " Finished cleaning."
