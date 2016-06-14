######################################################################
# Makefile for jayci
# --------------------------------------------------------------------
# By Christopher L Malbon
# Yarkony Group
# Department of Chemistry
# The Johns Hopkins University
#
# Date		Version
# -----------	-------
# 2016-06-10:     1.0.0
# 2016-06-13:     1.0.1
######################################################################

# Release version
JAYCIVER := 1.0.1

# Get OS name and version
UNAME	:= $(shell uname -a)
OS	:= $(word 1,$(UNAME))
OSV	:= $(word 3,$(UNAME))
ARC	:= $(shell uname -m)

$(info Operating system: $(OS))
$(info OS Version: $(OSV))

# Set up directories
bindir	:= bin
libdir	:= lib
tstdir	:= test
srcdir	:= source
JDIR	:= $(shell pwd)
BDIR	:= $(JDIR)/$(bindir)
LDIR	:= $(JDIR)/$(libdir)
SDIR	:= $(JDIR)/$(srcdir)
TDIR	:= $(JDIR)/$(tstdir)
IDIR	:= $(SDIR)/include
UNIXDIR := $(SDIR)/UNIX
COLIBDIR:= $(SDIR)/colib

# Set up compilers 
# In this section we assume 3 possibilities: 
#  1) We are building on NERSC machine edison.
#  2) We are building on the NERSC machine cori.
#  3) We are building on a standard Linux machine.
CORI:=cori
EDISON:=edison
ifneq ($(filter cori edison,$(NERSC_HOST)),)
	CC := cc -fopenmp -I $(IDIR) 
	FC := ftn
	AR := ar rv
	RANL := ranlib
else
	CC := gcc -fopenmp -I $(IDIR)
	FC := gfortran
	AR := ar rv
	RANL := ranlib
endif

# Compiler options 
FCOPS  :=  
CCDEF  := -DFLUSH -DINT64 -DEXTNAME
CPOPS  := 

# Compiler flags
# These are compiler dependent. Changing FFLAGS from
# the default value is not recommended. 
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

# Debugging flags
ifeq ($(findstring gfortran,$(FC)),)
# intel compilers
	DEBUG  := -g -DDEBUGGING $(VTUNE)
	FDEBUG := -traceback
else
# gnu compilers
	DEBUG  := -g -DDEBUGGING
	FDEBUG := -fbacktrace
endif

# Libraries
ifneq ($(filter cori edison,$(NERSC_HOST)),)
	MATHLIBS := 
else
	MATHLIBS := -L/usr/lib64 -llapack -lblas -lm -lgfortran -lgomp
endif

COLIBLIB := $(LDIR)/colib.a

# Objects for jayci
OBJS := errorlib.o \
	allocate_mem.o \
	cmdline.o \
	iminmax.o \
	arrayutil.o \
	bitutil.o \
	binary.o \
	progheader.o \
	combinatorial.o \
	straddress.o \
	moindex.o \
	readmoints.o \
	readnamelist.o \
	ioutil.o \
	abecalc.o \
	citruncate.o \
	binarystr.o \
	dlamch_fcn.o \
	mathutil.o \
	cimapping.o \
	action_util.o \
	prediagfcns.o \
	genbindet.o \
	det2string.o \
	initguess_sbd.o \
	davidson.o \
	execute_ci_calculation.o \
	run_jayci.o
# Objects for colib library 
COLIBO:=blaswrapper.o colib1.o colib2.o colib3.o colib4.o colib5.o colib6.o \
	colib7.o colib8.o colib9.o colib10.o
# Source files for colib lirary
COLIBF:=blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f colib6.f \
	colib7.f colib8.f colib9.f colib10.f 
UNIXO := fdate.o falloc.o fwtime.o hostnm.o flushstdout.o fsize.o
UNIXC := fdate.c falloc.c fwtime.c hostnm.c flushstdout.c fsize.c
TESTO := $(OBJS) test.o
JEXPO := $(OBJS) jayci_exp.o
JYCIO := $(OBJS) jayci.o

JEXPOBJS := $(addprefix $(SDIR)/,$(JEXPO))
JYCIOBJS := $(addprefix $(SDIR)/,$(JYCIO))
TESTOBJS := $(addprefix $(SDIR)/,$(TESTO))
COLIBOBJS:= $(addprefix $(COLIBDIR)/,$(COLIBO))
COLIBSRCF:= $(addprefix $(COLIBDIR)/,$(COLIBF))
UNIXOBJS := $(addprefix $(UNIXDIR)/,$(UNIXO))
UNIXSRCC := $(addprefix $(UNIXDIR)/,$(UNIXC))

# Set up executable names
JCIEXE := $(BDIR)/jayci-$(JAYCIVER)-$(OS)-$(ARC)
JXPEXE := $(BDIR)/jayci_exp-$(JAYCIVER)-$(OS)-$(ARC)
TESTEXE:= $(TDIR)/test.x
COLIBF := $(LDIR)/colib-$(JAYCIVER)-$(OS)-$(ARC).a

CDS := cd $(SDIR)
RM  := rm -rf

# Build --------------------------------------------------------------
all: colib jayci_exp jayci
	@echo "Finished building jayci."
	@echo ""

colib: $(COLIBOBJS) $(UNIXOBJS) | $(LDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   COLIB LIBRARY "
	@echo " Program Version:        $(JAYCIVER)"
	@echo " Archiver (AR):		$(AR)"
	@echo " Library file name:	$(COLIBF)"
	@echo "------------------------------------------------------"
	$(CDS); $(AR) $(COLIBF) $(COLIBOBJS) $(UNIXOBJS)
	$(RANL) $(COLIBF)
	chmod +x $(COLIBF)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new library"
	ln -sf $(COLIBF) $(LDIR)/colib.a
	@echo "------------------------------------------------------"
	@echo " Finished build."
	@echo ""

test: $(TESTOBJS) | $(TDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   JAYCI TESTING PROGRAM "
	@echo " Program version:	$(JAYCIVER)"
	@echo " Test program:		$(TESTEXE)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options:	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDS); $(CC) -o $(TESTEXE) $(TESTOBJS) $(MATHLIBS) $(COLIBLIB)
	@echo "------------------------------------------------------"
	@echo " Finished build."
	@echo ""

jayci_exp: $(JEXPOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   JAYCI_EXP PROGRAM "
	@echo " Program version:	$(JAYCIVER)"
	@echo " Test program:		$(TESTEXE)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options:	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDS); $(CC) -o $(JXPEXE) $(JEXPOBJS) $(MATHLIBS) $(COLIBLIB) 
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new binary"
	ln -sf $(JXPEXE) $(BDIR)/jayci_exp
	@echo "------------------------------------------------------"
	@echo " Finished build."
	@echo ""

jayci: $(JYCIOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   JAYCI PROGRAM "
	@echo " Program version:	$(JAYCIVER)"
	@echo " Test program:		$(TESTEXE)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options:	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDS); $(CC) -o $(JCIEXE) $(JYCIOBJS) $(MATHLIBS) $(COLIBLIB) $(DEBUG) $(CFLAGS)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new binary"
	ln -sf $(JCIEXE) $(BDIR)/jayci
	@echo "------------------------------------------------------"
	@echo " Finished build."	
	@echo ""

# Clean --------------------------------------------------------------
clean:
	rm -rf $(JYCIOBJS) 
	rm -rf $(SDIR)test.o $(SDIR)/jayci_exp.o $(SDIR)/jayci.o 
	rm -rf $(COLIBOBJS) $(UNIXOBJS)

# Rules --------------------------------------------------------------
$(UNIXDIR)/%.o: $(UNIXDIR)/%.c
	$(CC) -c $(CCDEF) -o $@ $< $(CFLAGS) -Wno-implicit-function-declaration

$(COLIBDIR)/%.o: $(COLIBDIR)/%.f
	$(FC) -c -o $@ $< $(CPOPS) $(DEBUG) $(FFLAGS)

$(SDIR)/%.o:$(SDIR)/%.c
	$(CC) -c -o $@ $< $(DEBUG) $(CFLAGS)
	@echo ""

$(SDIR)/dlamch_fcn.o:$(SDIR)/dlamch_fcn.f90
	$(FC) -c -o $(SDIR)/dlamch_fcn.o $(SDIR)/dlamch_fcn.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(SDIR)/readmoints.o:$(SDIR)/readmoints.f90
	$(FC) -c -o $(SDIR)/readmoints.o $(SDIR)/readmoints.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(SDIR)/readnamelist.o:$(SDIR)/readnamelist.f90
	$(FC) -c -o $(SDIR)/readnamelist.o $(SDIR)/readnamelist.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(SDIR)/ioutil.o:$(SDIR)/ioutil.c
	$(CC) -c -o $(SDIR)/ioutil.o $(SDIR)/ioutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(SDIR)/mathutil.o:$(SDIR)/mathutil.c
	$(CC) -c -o $(SDIR)/mathutil.o $(SDIR)/mathutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(BDIR) $(LDIR):
	@echo "Creating directory $@"
	mkdir -p $@
