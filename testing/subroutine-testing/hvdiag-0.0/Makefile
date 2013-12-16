# MAKEFILE for tester
PREPROC2	= -DHVDIAG

OBJS2			= detci1.o detci2.o detci5.o truncation.o citrunc.o hvdiag.o iwfmt.o


OBJS			= $(OBJS2) test.o

COMPILER	= ifort -i8
PREFIX		= sh -c
LIBS			= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIB2			= ./colib.a ./blaswrapper.a
CPOPT			= -auto -c -assume byterecl -g -traceback -fp-stack-check -check bounds
LKOPT			= -auto -lpthread
EXEC			= test.x

RM				= rm -f

# BUILD
test:$(OBJS)
	@echo "Building tester."
	@echo "Building $@"
	$(PREFIX) "$(COMPILER) -o $(EXEC) $(OBJS) $(LIB2) $(LIBS)"

%.o:./%.f
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -o $@ $< $(CPOPT)"

%.o:./%.f90
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -cpp -o $@ $< $(CPOPT) $(PREPROC2)"
  
clean:
	$(RM) $(OBJS) ./*.o ./*.mod ./*.dets ./*.list ./*.ref
	@echo "Finished cleaning."  

