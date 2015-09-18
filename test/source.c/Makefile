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

# Debugging ----------------------------------------------------------
DEBUG  := -g
FDEBUG := -fbacktrace
CDEBUG :=
# Libraries ----------------------------------------------------------
MATHLIBS := -L/usr/lib64 -llapack64 -lblas -lm
COLIBLIB := ./colib.a

# Executables --------------------------------------------------------
TESTX := test.x

# Objects ------------------------------------------------------------
OBJS  := binary.o combinatorial.o straddress.o moindex.o readmoints.o ioutil.o
TESTO := $(OBJS) test.o

# Build --------------------------------------------------------------
test: $(TESTO)
	$(CC) -o $(TESTX) $(TESTO) $(MATHLIBS) $(COLIBLIB) -lgfortran

# Clean --------------------------------------------------------------
clean:
	rm ./*.o

# Rules --------------------------------------------------------------
%.o:%.c
	$(CC) -c $< $(DEBUG)

readmoints.o:readmoints.f90
	$(FC) -c readmoints.f90 $(DEBUG) $(FDEBUG)
