######################################################################
# Makefile for jayci
#
#
######################################################################

UNIXDIR := ./UNIX
COLIBDIR:= ./colib

# Compilers ----------------------------------------------------------
ifeq ($(NERSC_HOST),edison)
	CC := cc -fopenmp -I ./include 
	FC := ftn
	AR := ar rv
	RANL := ranlib
else
	ifeq ($(NERSC_HOST),cori)
		CC := gcc -fopenmp -I ./include
		FC := gfortran
		AR := ar rv
		RANL:= ranlib
	else
	# Assuming standard Linux
		CC := gcc -fopenmp -I ./include
		FC := gfortran
		AR := ar rv
		RANL := ranlib
	endif
endif

# Options ------------------------------------------------------------
FCOPS  :=  
CCDEF  := -DFLUSH -DINT64 -DEXTNAME
CPOPS  := 

# FLAGS --------------------------------------------------------------
ifeq ($(findstring gfortran,$(FC)),)
# intel compilers
	FFLAGS := -i8 -auto -assume byterecl -O3
	CFLAGS := -Wall -std=c11 -march=native -funroll-all-loops \
		  -fomit-frame-pointer -fstrict-aliasing -O3
	VTUNE := -g -dynamic
else
# gnu compilers
	FFLAGS := -fdefault-integer-8 -frecord-marker=4 -O3
	CFLAGS := -Wall -std=c11 -march=native -funroll-loops -ffast-math \
	          -fomit-frame-pointer -fstrict-aliasing -O3

endif
# Debugging ----------------------------------------------------------
ifeq ($(findstring gfortran,$(FC)),)
# intel compilers
	DEBUG  := -g -DDEBUGGING $(VTUNE)
	FDEBUG := -traceback
else
# gnu compilers
	DEBUG  := -g -DDEBUGGING
	FDEBUG := -fbacktrace
endif

# Libraries ----------------------------------------------------------
ifeq ($(NERSC_HOST),edison)
	MATHLIBS := 
else
	ifeq ($(findstring BSD,$(shell uname -s)),)
	# building on linux
		MATHLIBS := -L/usr/lib64 -llapack -lblas -lm -lgfortran
	else
		MATHLIBS := -L/usr/local/lib -llapack -lblas \
				-L/usr/local/lib/gcc48 -lgfortran -lm
	endif
endif
COLIBLIB := ./colib.a

# Executables --------------------------------------------------------
TESTX := test.x
JEXPX := jayci_exp.x
JYCIX := jayci.x

# Objects ------------------------------------------------------------
OBJS  := errorlib.o cmdline.o iminmax.o arrayutil.o bitutil.o binary.o \
	 progheader.o combinatorial.o \
	 straddress.o moindex.o readmoints.o readnamelist.o ioutil.o \
	 abecalc.o citruncate.o binarystr.o dlamch_fcn.o mathutil.o \
	 action_util.o prediagfcns.o

COLIBO := blaswrapper.o colib1.o colib2.o colib3.o colib4.o colib5.o \
	  colib6.o colib7.o colib8.o colib9.o colib10.o
COLIBF := blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f \
          colib6.f colib7.f colib8.f colib9.f colib10.f
UNIXO := fdate.o falloc.o fwtime.o hostnm.o flushstdout.o fsize.o
UNIXC := fdate.c falloc.c fwtime.c hostnm.c flushstdout.c fsize.c
TESTO := $(OBJS) test.o
JEXPO := $(OBJS) jayci_exp.o
JYCIO := $(OBJS) jayci.o


# Build --------------------------------------------------------------
colib: $(COLIBO) $(UNIXO)
	$(AR) colib.a $(COLIBO) $(UNIXO)
	$(RANL) colib.a 
	chmod +x colib.a

test: $(TESTO)
	$(CC) -o $(TESTX) $(TESTO) $(MATHLIBS) $(COLIBLIB) 

jayci_exp: $(JEXPO)
	$(CC) -o $(JEXPX) $(JEXPO) $(MATHLIBS) $(COLIBLIB) 

jayci: $(JYCIO)
	@echo "-----------------------------------------------------------"
	@echo "=>Building: jayci.x executable."
	@echo "  "
	@echo " MATHLIBS = $(MATHLIBS)"
	@echo " COLIBLIB = $(COLIBLIB)"
	@echo " DEBUG    = $(DEBUG)"
	@echo ""
	$(CC) -o $(JYCIX) $(JYCIO) $(MATHLIBS) $(COLIBLIB) $(DEBUG) $(CFLAGS)
	@echo " Build complete."
	@echo ""

# Clean --------------------------------------------------------------
clean:
	rm -rf $(OBJS) test.o jayci_exp.o jayci.o $(COLIBO) $(UNIXO)

# Rules --------------------------------------------------------------
%.o: $(UNIXDIR)/%.c
	@echo "->Building: $@"
	$(CC) -c $(CCDEF) -o $@ $< $(CFLAGS) -Wno-implicit-function-declaration

%.o: $(COLIBDIR)/%.f
	@echo "->Building $@"
	$(FC) -c -o $@ $< $(CPOPS) $(DEBUG) $(FFLAGS)

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
