# Makefile for driver1
#
NAME			:= $(shell hostname)
JDIR			:= $(shell pwd)
SOURCE		:= $(JDIR)/src

# Set any preprocessing flags
PREPROC		:=

# Set any Debugging flags
DEBUGFLG	:= -g -traceback -fp-stack-check -check bounds

# Set any compiler options
CPOPT			:= -auto -c -assume byterecl -O3

# Set compiler
COMPILER	:= ifort -i8
PREFIX		:= sh -c

# Set libraries
LIBS			:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIBS2			:= $(JDIR)/colib.a $(JDIR)/blaswrapper.a

# Linking options
LKOPT			:= -auto -lpthread

# Objects for building
OBJECTS		:= detci1.o detci2.o detci5.o iwfmt.o possex1.o cannon.o singrepinfo.o doublerepinfo.o \
             hvdiag.o actionutil.o acthv.o orthogroutines.o david_util.o davidson.o cannonreal.o \
             lowdiags.o \
             construct.o prediag.o driver2.o

# Executable
EXEC1			:= driver.x

# Cleaning
RM				:= rm -f

#Building
driver:$(OBJECTS)
	@echo "Bulding on $(NAME) with: $(PREPROC)"
	@echo "Building $@"
	$(PREFIX) "$(COMPILER) -o $(EXEC1) $(OBJECTS) $(LIBS) $(LIBS2) $(DEGUBFLG)"
	@echo "Finished building $@"

%.o:$(SOURCE)/%.f90
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -cpp -o $@ $< $(CPOPT) $(DEBUGFLG)"
	@echo "Finished building $<"

%.o:$(SOURCE)/%.f
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -cpp -o $@ $< $(CPOPT) $(DEBUGFLG)"
	@echo "Finished building $<"

clean:
	$(RM) $(OBJECTS) ./*.o ./*.mod
	@echo "Finished cleaning."
