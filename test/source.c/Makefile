######################################################################
# Makefile for jayci
#
#
######################################################################

# Compilers ----------------------------------------------------------
CC := gcc
FC := gfortran -fdefault-integer-8

# Options ------------------------------------------------------------
FCOPS  := -c -frecord-marker=4 -O3
CCDEF  := -c -DFLUSH -DINT64 -DEXTNAME
CPOPS  := 

# Debugging ----------------------------------------------------------
DEBUG  := -g -DDEBUGGING
FDEBUG := -fbacktrace
CDEBUG :=

# Libraries ----------------------------------------------------------
MATHLIBS := -L/usr/lib64 -llapack64 -lblas -lm
COLIBLIB := ./colib.a

# Executables --------------------------------------------------------
TESTX := test.x
JEXPX := jayci_exp.x
JYCIX := jayci.x

# Objects ------------------------------------------------------------
OBJS  := bitutil.o binary.o initarray.o combinatorial.o straddress.o moindex.o \
	 readmoints.o readnamelist.o ioutil.o abecalc.o citruncate.o \
	 binarystr.o dlamch_fcn.o mathutil.o action_util.o prediagfcns.o
TESTO := $(OBJS) test.o
JEXPO := $(OBJS) jayci_exp.o
JYCIO := $(OBJS) jayci.o

# Build --------------------------------------------------------------
test: $(TESTO)
	$(CC) -o $(TESTX) $(TESTO) $(MATHLIBS) $(COLIBLIB) -lgfortran

jayci_exp: $(JEXPO)
	$(CC) -o $(JEXPX) $(JEXPO) $(MATHLIBS) $(COLIBLIB) -lgfortran

jayci: $(JYCIO)
	$(CC) -o $(JYCIX) $(JYCIO) $(MATHLIBS) $(COLIBLIB) $(DEBUG) -lgfortran

# Clean --------------------------------------------------------------
clean:
	rm -rf $(OBJS) test.o jayci_exp.o jayci.o

# Rules --------------------------------------------------------------
%.o:%.c
	$(CC) -c $< $(DEBUG)

dlamch_fcn.o:dlamch_fcn.f90
	$(FC) -c dlamch_fcn.f90 $(DEBUG) $(FDEBUG)

readmoints.o:readmoints.f90
	$(FC) -c readmoints.f90 $(DEBUG) $(FDEBUG)

readnamelist.o:readnamelist.f90
	$(FC) -c readnamelist.f90 $(DEBUG) $(FDEBUG)

ioutil.o:ioutil.c
	$(CC) -c ioutil.c $(DEBUG) -Wno-implicit-function-declaration

mathutil.o:mathutil.c
	$(CC) -c mathutil.c $(DEBUG) -Wno-implicit-function-declaration
