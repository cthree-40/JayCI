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
# Or we are on Hyperion
else
	LIBS1	:= -L/usr/lib64 -llapack64 -lblas
endif

LIBS2		:= $(SOURCE)/colib.a

ifneq ($(findstring gfortran,$(COMPLR)),)
#GNU compilers
	CPOPS	:= -c -frecord-marker=4 -O3
	CLIBOPS	:= -c -frecord-marker=4 -O3
	DEBUG	:= -g -fcheck=all -g -fbacktrace
	CDEF	:= -c -DFLUSH -DINT64 -DEXTNAME
else
#Intel compilers
	CPOPS	:= -auto -c -assume byterecl -O3
	DEBUG	:= -debug -check all -check bounds -check uninit -gen-interfaces -warn interfaces -g -traceback
	LKOPT	:= -auto -lpthread
	VTUNE	:= -g -O3
endif

# Executables
EXEC1		:= $(BIN)/jayCI1.x
EXEC2		:= $(BIN)/jayCI2.x
TESTEX		:= $(BIN)/test.x

# Clean up
REMOVE	:= rm -f

#######################################
# Objects

COLIBO	:= blaswrapper.o colib1.o colib2.o colib3.o colib4.o colib5.o colib6.o colib7.o \
			colib8.o colib9.o colib10.o
COLIBF	:= blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f colib6.f colib7.f \
			colib8.f colib9.f colib10.f
UNIXO	:= fdate.o falloc.o fwtime.o hostnm.o flushstdout.o fsize.o
UNIXC	:= fdate.c falloc.c fwtime.c hostnm.c flushstdout.c fsize.c
OBJS1		:= common_util.o truncation.o citrunc.o
OBJS2		:= $(OBJS1) construct.o action_util.1.o acthv.o david_util.o initial_guess.o Davidson.o
OBJSTEST	:= common_util.o construct.o action_util.1.o acthv.o test.o
DRIV1		:= $(OBJS1) jayCI1.o
DRIV2		:= $(OBJS2) jayCI2.o

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
	@echo " Cleaning... "
	$(RM) $(DRIV1) ./*genmod.f90 ./*.o ./*.mod


driver2:$(DRIV2) | $(BIN)
	@echo " Building driver2..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(EXEC2) $(DRIV2) $(LIBS1) $(LIBS2) $(VTUNE) "
	@echo " ------------------------------------------------------------------------"

test:	$(OBJSTEST) | $(BIN)
	@echo " Building test..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(TESTEX) $(OBJSTEST) $(LIBS1) $(LIBS2) $(DEBUG)"
	@echo " ------------------------------------------------------------------------"
	@echo " Cleaning... "
	$(RM) $(OBJSTEST) ./*__genmod.f90 ./*.o ./*.mod

colib:	$(COLIBO) $(UNIXO)
	$(AR) colib.a $(COLIBO) $(UNIXO)
	$(RANL) colib.a

%.o:$(SOURCE)/$(UNIXDIR)/%.c
	@echo " Building $< "
	$(CC) $(CDEF) -o $@ $<

%.o:$(SOURCE)/$(COLIBDIR)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COLIBFC) -o $@ $< $(CPOPS) $(VTUNE) $(DEBUG)"
	chmod +x colib.a

%.o:$(SOURCE)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -o $@ $< $(CPOPS) $(VTUNE)"

%.o:$(SOURCE)/%.f90
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -cpp -o $@ $< $(CPOPS) $(DEBUG) $(PREPROC) $(VTUNE)"

clean1:
	$(RM) $(DRIV1) ./*__genmod.f90./*.o ./*.mod 
	@echo " Finished cleaning."


clean2:
	$(RM) $(DRIV2) ./*.o ./*.mod ./*__genmod.f90
	@echo " Finished cleaning."

cleancolib:
	$(RM) ./blaswrapper.o ./colib/colib*.o ./colib*.o ./fwtime.o ./hostnm.o ./falloc.o ./fdate.o ./fsize.o ./flushstdout.o
	@echo " Finished cleaning."
