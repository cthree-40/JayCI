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
BIN		:= $(JOBDIR)/bin

# Set compiler
# We build on edison
ifeq ($(NERSC_HOST),edison)
	COMPLR	:= ftn -i8
# Or we are on Hyperion
else
	COMPLR	:= ifort -i8
endif

PREFIX		:= sh -c

# Set Libraries
# We build on edison
ifeq ($(NERSC_HOST),edison)
	LIBS1	:= 
# Or we are on Hyperion
else
	LIBS1	:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif

LIBS2		:= $(SOURCE)/colib.a $(SOURCE)/blaswrapper.a
CPOPS		:= -auto -c -assume byterecl -O3 
DEBUG		:= -debug -check all -check bounds -check uninit -gen-interfaces -warn interfaces -g -traceback
LKOPT		:= -auto -lpthread
VTUNE		:= -g -O3

# Executables
EXEC1		:= $(BIN)/jayCI1.x
EXEC2		:= $(BIN)/jayCI2.x
TESTEX		:= $(BIN)/test.x

# Clean up
REMOVE	:= rm -f

#######################################
# Objects

OBJS1		:= common_util.o truncation.o citrunc.o
OBJS2		:= $(OBJS1) construct.o action_util.1.o acthv.o david_util.o initial_guess.o Davidson.o \
			iwfmt.o
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
	$(RM) $(DRIV1) ./genmod.f90 ./*.o ./*.mod


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
	$(RM) $(OBJSTEST) ./*genmod.f90 ./*.o ./*.mod

%.o:$(SOURCE)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -o $@ $< $(CPOPS) $(VTUNE)"

%.o:$(SOURCE)/%.f90
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -cpp -o $@ $< $(CPOPS) $(DEBUG) $(PREPROC) $(VTUNE)"

clean1:
	$(RM) $(DRIV1) ./*genmod.f90./*.o ./*.mod 
	@echo " Finished cleaning."


clean2:
	$(RM) $(DRIV2) ./*.o ./*.mod ./*genmod.f90
	@echo " Finished cleaning."
