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

# Linking options
LKOPT			:= -auto -lpthread

# Objects for building
OBJECTS		:= detci1.o detci2.o truncation.o citrunc.o driver1.o

# Executable
EXEC1			:= driver.x

# Cleaning
RM				:= rm -f

#Building
driver:$(OBJECTS)
	@echo "Bulding on $(NAME) with: $(PREPROC)"
	@echo "Building $@"
	$(PREFIX) "$(COMPILER) -o $(EXEC1) $(OBJECTS) $(LIBS)"
	@echo "Finished building $@"

%.o:$(SOURCE)/%.f90
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -cpp -o $@ $< $(CPOPT)"
	@echo "Finished building $<"

clean:
	$(RM) $(OBJECTS) ./*.o ./*.mod ./*.dets ./*.ref ./*.out ./*.list
	@echo "Finished cleaning."
