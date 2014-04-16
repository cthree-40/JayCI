## Makefile for Determinant CI
#
## By Christopher L. Malbon
## Yarkony Group
## The Johns Hopkins University
## 
#

#######################################
# Build environment
JOBDIR	:= $(shell pwd)
SOURCE	:= $(JOBDIR)/source
BIN		:= $(JOBDIR)/bin
COMPLR	:= ifort -i8
PREFIX	:= sh -c
LIBS1		:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIBS2		:= $(SOURCE)/colib.a $(SOURCE)/blaswrapper.a
CPOPS		:= -auto -c -assume byterecl -O3
DEBUG		:= -debug -gen-interfaces -warn interfaces -g -traceback -fp-stack-check -check bounds -check all
LKOPT		:= -auto -lpthread

# Executables
EXEC1		:= $(BIN)/driver1.x
EXEC2		:= $(BIN)/driver2.x

# Clean up
REMOVE	:= rm -f

#######################################
# Objects

OBJS1		:= detci1.o detci2.o detci5.o truncation.o citrunc.o iwfmt.o possex1.o cannon.o \
						singrepinfo.o doublerepinfo.o construct.o hvdiag.o actionutil.o acthv.o
OBJS2		:= $(OBJS1) 
OBJS3		:= $(OBJS2) cannonreal.o lowdiags.o prediag.o orthogroutines.o david_util.o davidson.o

DRIV1		:= $(OBJS2) driver1.o
DRIV2		:= $(OBJS3) driver2.o

#######################################
# Preprocessor flags

PREPROC	:= 
PROFILE	:= -prof-gen -prof-dir$(JOBDIR)

#######################################
# BUILD
driver1:	$(DRIV1) | $(BIN)
	@echo " Building driver1..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(EXEC1) $(DRIV1) $(LIBS1) $(LIBS2) "
	@echo " ------------------------------------------------------------------------"
	@echo " Cleaning... "
	$(RM) $(DRIV1) ./genmod.f90 ./*.o ./*.mod


driver2:$(DRIV2) | $(BIN)
	@echo " Building driver1..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building with $(PREPROC)..."
	@echo " ------------------------------------------------------------------------"
	@echo " Building $@..."
	@echo " --- "
	$(PREFIX) "$(COMPLR) -o $(EXEC2) $(DRIV2) $(LIBS1) $(LIBS2) $(DEBUG)"
	@echo " ------------------------------------------------------------------------"

%.o:$(SOURCE)/%.f
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -o $@ $< $(CPOPS)"

%.o:$(SOURCE)/%.f90
	@echo " Building $< "
	$(PREFIX) "$(COMPLR) -cpp -o $@ $< $(CPOPS) $(PREPROC) $(DEBUG)"

clean1:
	$(RM) $(DRIV1) ./*genmod.f90./*.o ./*.mod 
	@echo " Finished cleaning."


clean2:
	$(RM) $(DRIV2) ./*.o ./*.mod ./*genmod.f90
	@echo " Finished cleaning."
