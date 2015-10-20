######################################################################
# Makefile for jayci
#
#
######################################################################

# Compilers ----------------------------------------------------------
ifeq ($(NERSC_HOST),edison)
	CC := cc
	FC := ftn
else
	CC := gcc 
	FC := gfortran 
endif

# Options ------------------------------------------------------------
FCOPS  :=  
CCDEF  := -DFLUSH -DINT64 -DEXTNAME
CPOPS  := 

# FLAGS --------------------------------------------------------------
FFLAGS := -fdefault-integer-8 -frecord-marker=4 -O3
CFLAGS := -march=native -funroll-loops -ffast-math -fomit-frame-pointer \
	  -fstrict-aliasing

# Debugging ----------------------------------------------------------
DEBUG  := -g -DDEBUGGING
FDEBUG := -fbacktrace
CDEBUG :=

# Libraries ----------------------------------------------------------
ifeq ($(NERSC_HOST),edison)
	MATHLIBS := $(BLAS_LIB)
else
	MATHLIBS := -L/usr/lib64 -llapack -lblas -lm
endif
COLIBLIB := ./colib.a

# Executables --------------------------------------------------------
TESTX := test.x
JEXPX := jayci_exp.x
JYCIX := jayci.x

# Objects ------------------------------------------------------------
OBJS  := arrayutil.o bitutil.o binary.o combinatorial.o \
	 straddress.o moindex.o readmoints.o readnamelist.o ioutil.o \
	 abecalc.o citruncate.o binarystr.o dlamch_fcn.o mathutil.o \
	 action_util.o prediagfcns.o
TESTO := $(OBJS) test.o
JEXPO := $(OBJS) jayci_exp.o
JYCIO := $(OBJS) jayci.o

# Build --------------------------------------------------------------
test: $(TESTO)
	$(CC) -o $(TESTX) $(TESTO) $(MATHLIBS) $(COLIBLIB) -lgfortran

jayci_exp: $(JEXPO)
	$(CC) -o $(JEXPX) $(JEXPO) $(MATHLIBS) $(COLIBLIB) -lgfortran

jayci: $(JYCIO)
	@echo "-----------------------------------------------------------"
	@echo "=>Building: jayci.x executable."
	@echo "  "
	@echo " MATHLIBS = $(MATHLIBS)"
	@echo " COLIBLIB = $(COLIBLIB)"
	@echo " DEBUG    = $(DEBUG)"
	@echo ""
	$(CC) -o $(JYCIX) $(JYCIO) $(MATHLIBS) $(COLIBLIB) $(DEBUG) $(CFLAGS) -lgfortran
	@echo " Build complete."
	@echo ""

# Clean --------------------------------------------------------------
clean:
	rm -rf $(OBJS) test.o jayci_exp.o jayci.o

# Rules --------------------------------------------------------------
%.o:%.c
	@echo "->Building: $@ "
	$(CC) -c $< $(DEBUG) $(CFLAGS)
	@echo ""

dlamch_fcn.o:dlamch_fcn.f90
	@echo "->Building: $@ "
	$(FC) -c dlamch_fcn.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

readmoints.o:readmoints.f90
	@echo "->Building: $@ "
	$(FC) -c readmoints.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

readnamelist.o:readnamelist.f90
	@echo "->Building: $@ "
	$(FC) -c readnamelist.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

ioutil.o:ioutil.c
	@echo "->Building: $@ "
	$(CC) -c ioutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

mathutil.o:mathutil.c
	@echo "->Building: $@ "
	$(CC) -c mathutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""
