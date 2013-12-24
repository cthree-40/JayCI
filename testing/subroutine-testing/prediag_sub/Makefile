# MAKEFILE for tester
PREPROC2	=-DPREDIAG -DHVDIAG -DHVACTION -DEXPBUILD 
SOURCE		= ./src
JDIR			= ./
OBJS2			= detci1.o detci2.o detci5.o truncation.o citrunc.o hvdiag.o iwfmt.o \
						possex1.o cannon.o singrepinfo.o doublerepinfo.o actionutil.o \
						acthv.o
OBJS3			= $(OBJS2) construct.o

OBJS			= $(OBJS3) cannonreal.o lowdiags.o prediag.o test.o

COMPILER	= ifort -i8
PREFIX		= sh -c
LIBS			= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm
LIB2			= ./colib.a ./blaswrapper.a
CPOPT			= -auto -c -assume byterecl -g -traceback -fp-stack-check -check bounds -check all
LKOPT			= -auto -lpthread
EXEC			= test.x

RM				= rm -f

# BUILD
test:$(OBJS)
	@echo "Building tester."
	@echo " "
	@echo "-------------------------------------------------"
	@echo " "
	@echo " Building with $(PREPROC2)"
	@echo " "
	@echo "-------------------------------------------------"
	@echo " "
	@echo "Building $@"
	$(PREFIX) "$(COMPILER) -o $(EXEC) $(OBJS) $(LIB2) $(LIBS)"

%.o:$(SOURCE)/%.f
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -o $@ $< $(CPOPT)"

%.o:$(SOURCE)/%.f90
	@echo "Building $<"
	$(PREFIX) "$(COMPILER) -cpp -o $@ $< $(CPOPT) $(PREPROC2)"
  
clean:
	$(RM) $(OBJS) ./problem* ./debug.acthv ./*.o ./*.mod ./*.dets ./*.list ./*.ref ./alpha-single.cont ./diag.cont ./output_acthv.1 ./alpha.cont
	@echo "Finished cleaning."  

