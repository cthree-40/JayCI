######################################################################
# Makefile for jayci, pjayci, and dycicalc
# --------------------------------------------------------------------
# By Christopher L Malbon
# Yarkony Group
# Department of Chemistry
# The Johns Hopkins University
#
# JAYCI ----------------
# Date		Version
# -----------	-------
# 2016-06-10:     1.0.0
# 2016-06-13:     1.0.1
# 2016-06-14:     1.0.2
# 2016-06-19:     1.1.0
# 2018-01-09:     2.0.0
# 2018-01-26:     2.0.1
# ----------------------
#
# PJAYCI ---------------
# Date          Version
# -----------   -------
# 
# ----------------------
######################################################################

# Release version
JAYCIVER := 2.0.1
PJAYCIVER:= 1.2.0
DYCICALCVER:= 1.0.0
PDYCICALCVER:= 1.0.0

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
mpisrc	:= mpi_source
JDIR	:= $(shell pwd)
BDIR	:= $(JDIR)/$(bindir)
LDIR	:= $(JDIR)/$(libdir)
SDIR	:= $(JDIR)/$(srcdir)
MPISDIR := $(JDIR)/$(mpisrc)
DYCISDIR:= $(JDIR)/$(dycisrc)
TDIR	:= $(JDIR)/$(tstdir)
IDIR	:= $(SDIR)/include
MPIIDIR := $(MPISDIR)/include
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
	CC := cc -qopenmp -I $(IDIR) 
	MPICC := cc -qopenmp -I $(MPIIDIR) -DBIGMEM
	FC := ftn
	MPIFC := ftn -DBIGMEM
	AR := ar rv
	RANL := ranlib
else
	ifneq ($(findstring intel,$(COMPILERS)),)
		CC := icc -qopenmp -I${MKLROOT}/include -I $(IDIR) 
		MPICC := mpicc -qopenmp -I${MKLROOT}/include -I $(MPIIDIR)
		#FC := ifort -I${MKLROOT}/include
		FC := gfortran
		MPIFC := mpif90
		AR := ar rv
		RANL := ranlib
	else
		CC := gcc -fopenmp -I $(IDIR)
		MPICC := mpicc -fopenmp -I $(MPIIDIR)
		FC := gfortran
		MPIFC := mpif90
		AR := ar rv
		RANL := ranlib
	endif
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
	CFLAGS := -Wall -std=c11 -funroll-all-loops \
		  -fomit-frame-pointer -fstrict-aliasing -O3
	VTUNE := -g
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
	GALIBS   := -L/usr/local/lib -lga -larmci
endif

COLIBLIB := $(LDIR)/colib.a

# Objects for jayci
OBJS := 	timestamp.o \
		errorlib.o \
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
		readmocoef.o \
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
		initguess_roldv.o \
		davidson.o \
		write_wavefunction.o \
		read_wavefunction.o \
		execute_ci_calculation.o \
		run_jayci.o

# Objects for mpi_jayci
MPIOBJS :=	timestamp.o \
	   	errorlib.o \
		arrayutil.o \
		dlamch_fcn.o \
		mathutil.o \
		allocate_mem.o \
		readmoints.o \
		readnamelist.o \
		readmocoef.o \
		ioutil.o \
		mpi_utilities.o \
		abecalc.o \
		combinatorial.o \
		moindex.o \
		straddress.o \
		iminmax.o \
		bitutil.o \
		binary.o \
		binarystr.o \
		citruncate.o \
		action_util.o \
                pdavidson.o \
		execute_pjayci.o	

DYCIOBJS := 	errorlib.o \
		arrayutil.o \
		allocate_mem.o \
		abecalc.o \
		combinatorial.o \
		straddress.o \
		iminmax.o \
		bitutil.o \
		binary.o \
		binarystr.o \
		dlamch_fcn.o \
		mathutil.o \
		moindex.o \
		readmoints.o \
		readnamelist.o \
		readmocoef.o \
		ioutil.o \
		buildao.o \
		atomic_orbitals.o \
		read_wavefunction.o \
		write_wavefunction.o \
		citruncate.o \
		progheader.o \
		cimapping.o \
		action_util.o \
		dysoncomp.o \
		run_dycicalc.o
#		execute_dycicalc.o \
	    	run_dycicalc.o

PDYCIOBJS :=	timestamp.o \
		errorlib.o \
		arrayutil.o \
		allocate_mem.o \
		abecalc.o \
		straddress.o \
		iminmax.o \
		mpi_utilities.o \
		moindex.o \
		readmoints.o \
		readnamelist.o \
		bitutil.o \
		binary.o \
		binarystr.o \
		combinatorial.o \
		citruncate.o \
                ioutil.o \
		action_util.o \
		dysoncomp.o \
		run_pdycicalc.o

# Objects for AO evaluation test
AOTESTOBJS	:=	errorlib.o \
			mathutil.o \
			buildao.o \
			atomic_orbitals.o

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
MPIJYCIO := $(MPIOBJS) pjayci.o
DYCIO := $(DYCIOBJS) dycicalc.o
PDYCIO:= $(PDYCIOBJS) pdycicalc.o
AOTESTO := $(DYCIOBJS) test_aorbitals.o

JEXPOBJS := $(addprefix $(SDIR)/,$(JEXPO))
JYCIOBJS := $(addprefix $(SDIR)/,$(JYCIO))
PJYCIOBJS:= $(addprefix $(MPISDIR)/,$(MPIJYCIO))
DYCIOBJS := $(addprefix $(SDIR)/,$(DYCIO))
PDYCIOBJS:= $(addprefix $(MPISDIR)/,$(PDYCIO))
AOTESTOBJS:= $(addprefix $(SDIR)/,$(AOTESTO))
TESTOBJS := $(addprefix $(SDIR)/,$(TESTO))
COLIBOBJS:= $(addprefix $(COLIBDIR)/,$(COLIBO))
COLIBSRCF:= $(addprefix $(COLIBDIR)/,$(COLIBF))
UNIXOBJS := $(addprefix $(UNIXDIR)/,$(UNIXO))
UNIXSRCC := $(addprefix $(UNIXDIR)/,$(UNIXC))

# Set up executable names
JCIEXE := $(BDIR)/jayci-$(JAYCIVER)-$(OS)-$(ARC)
PJCIEXE:= $(BDIR)/pjayci-$(PJAYCIVER)-$(OS)-$(ARC)
JXPEXE := $(BDIR)/jayci_exp-$(JAYCIVER)-$(OS)-$(ARC)
DYCIEXE:= $(BDIR)/dycicalc-$(DYCICALCVER)-$(OS)-$(ARC)
PDYCIEXE:=$(BDIR)/pdycicalc-$(PDYCICALCVER)-$(OS)-$(ARC)
TESTEXE:= $(TDIR)/test.x
ATESTEXE:= $(BDIR)/testao.x
COLIBX := $(LDIR)/colib-$(JAYCIVER)-$(OS)-$(ARC).a
CDS := cd $(SDIR)
CDPS:= cd $(MPISDIR)
RM  := rm -rf

# Build --------------------------------------------------------------
all: colib jayci_exp jayci pjayci dycicalc pdycicalc
	@echo "Finished building jayci."
	@echo ""

colib: $(COLIBOBJS) $(UNIXOBJS) | $(LDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   COLIB LIBRARY "
	@echo " Program Version:        $(JAYCIVER)"
	@echo " Archiver (AR):		$(AR)"
	@echo " Library file name:	$(COLIBX)"
	@echo "------------------------------------------------------"
	$(CDS); $(AR) $(COLIBX) $(COLIBOBJS) $(UNIXOBJS)
	$(RANL) $(COLIBX)
	chmod +x $(COLIBX)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new library"
	ln -sf $(COLIBX) $(LDIR)/colib.a
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

testao: $(AOTESTOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "  TEST AO Program "
	$(CDS); $(CC) -o $(ATESTEXE) $(AOTESTOBJS) $(MATHLIBS) $(COLIBLIB)
	@echo "------------------------------------------------------"
	@echo " Finished"
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

pjayci: $(PJYCIOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   PJAYCI PROGRAM "
	@echo " Program version:	$(PJAYCIVER)"
	@echo " Test program:		$(TESTEXE)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options:	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDPS); $(MPICC) -o $(PJCIEXE) $(PJYCIOBJS) $(MATHLIBS) $(COLIBLIB) $(GALIBS) $(DEBUG) $(CFLAGS)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new binary"
	ln -sf $(PJCIEXE) $(BDIR)/pjayci
	@echo "------------------------------------------------------"
	@echo " Finished build."	
	@echo ""

dycicalc: $(DYCIOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   DYCICALC PROGRAM "
	@echo " Program version:	$(DYCICALCVER)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options: 	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDS); $(CC) -o $(DYCIEXE) $(DYCIOBJS) $(MATHLIBS) $(COLIBLIB) $(DEBUG) $(CFLAGS)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new binary"
	ln -sf $(DYCIEXE) $(BDIR)/dycicalc
	@echo "------------------------------------------------------"
	@echo " Finished build."
	@echo ""

pdycicalc: $(PDYCIOBJS) | $(BDIR)
	@echo ""
	@echo "------------------------------------------------------"
	@echo "   PARALLEL DYCICALC PROGRAM "
	@echo " Program version:	$(PDYCICALCVER)"
	@echo " BLAS/LAPACK Lib:	$(MATHLIBS)"
	@echo " COLIB library:		$(COLIBLIB)"
	@echo " Debug flags:		$(DEBUG)"
	@echo " C Compiler options: 	$(CFLAGS)"
	@echo " F90 Compiler options:	$(FFLAGS)"
	@echo "------------------------------------------------------"
	$(CDS); $(MPICC) -o $(PDYCIEXE) $(PDYCIOBJS) $(MATHLIBS) $(COLIBLIB) $(GALIBS) $(DEBUG) $(CFLAGS)
	@echo "------------------------------------------------------"
	@echo " Creating symbolic link to new binary"
	ln -sf $(PDYCIEXE) $(BDIR)/pdycicalc
	@echo "------------------------------------------------------"
	@echo " Finished build."
	@echo ""

# Clean --------------------------------------------------------------
clean:
	rm -rf $(JYCIOBJS) $(PJYCIOBJS) $(DYCIOBJS) $(PDYCIOBJS)
	rm -rf $(SDIR)/test.o $(MPISDIR)/test.o $(SDIR)/jayci_exp.o $(MPISDIR)/pjayci.o $(SDIR)/jayci.o $(MPISDIR)/pdycalc.o
	rm -rf $(SDIR)/dycicalc.o

deepclean:
	rm -rf $(JYCIOBJS) $(PJYCIOBJS) $(DYCIOBJS) $(PDYCIOBJS)
	rm -rf $(SDIR)/test.o $(MPISDIR)/test.o $(SDIR)/jayci_exp.o $(MPISDIR)/pjayci.o $(SDIR)/jayci.o $(MPISDIR)/pdycalc.o
	rm -rf $(SDIR)/dycicalc.o
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

$(SDIR)/readmocoef.o:$(SDIR)/readmocoef.f90
	$(FC) -c -o $(SDIR)/readmocoef.o $(SDIR)/readmocoef.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(SDIR)/ioutil.o:$(SDIR)/ioutil.c
	$(CC) -c -o $(SDIR)/ioutil.o $(SDIR)/ioutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(SDIR)/mathutil.o:$(SDIR)/mathutil.c
	$(CC) -c -o $(SDIR)/mathutil.o $(SDIR)/mathutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(MPISDIR)/%.o:$(MPISDIR)/%.c
	$(MPICC) -c -o $@ $< $(DEBUG) $(CFLAGS)
	@echo ""

$(MPISDIR)/readnamelist.o:$(MPISDIR)/readnamelist.f90
	$(MPIFC) -c -o $(MPISDIR)/readnamelist.o $(MPISDIR)/readnamelist.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(MPISDIR)/readmoints.o:$(MPISDIR)/readmoints.f90
	$(MPIFC) -c -o $(MPISDIR)/readmoints.o $(MPISDIR)/readmoints.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(MPISDIR)/readmocoef.o:$(MPISDIR)/readmocoef.f90
	$(MPIFC) -c -o $(MPISDIR)/readmocoef.o $(MPISDIR)/readmocoef.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(MPISDIR)/ioutil.o:$(MPISDIR)/ioutil.c
	$(MPICC) -c -o $(MPISDIR)/ioutil.o $(MPISDIR)/ioutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(MPISDIR)/dlamch_fcn.o:$(MPISDIR)/dlamch_fcn.f90
	$(FC) -c -o $(MPISDIR)/dlamch_fcn.o $(MPISDIR)/dlamch_fcn.f90 $(DEBUG) $(FDEBUG) $(FFLAGS)
	@echo ""

$(MPISDIR)/mathutil.o:$(MPISDIR)/mathutil.c
	$(MPICC) -c -o $(MPISDIR)/mathutil.o $(MPISDIR)/mathutil.c $(DEBUG) $(CFLAGS) -Wno-implicit-function-declaration
	@echo ""

$(BDIR) $(LDIR):
	@echo "Creating directory $@"
	mkdir -p $@
