# $Id$
#
# Makefile for compiling and linking LISFLOOD

# ************************************************************
# DEFAULT CONFIGURATION
# Default options for simple scalar compilation using g++
# for bluecrystal/dylan/any linux installation.
# ************************************************************

# Use "g++" to compile source files
CC = g++
# Use "g++" to link object files
LD = g++
# Compiler flags
CFLAGS = -w -g -O3
# -openmp-profile
# Linker flags
LDFLAGS = -w -g -O3


# ************************************************************
# OMP PARALLEL CONFIGURATION
# Options for openmp parallel compilation using icpc
# for bluecrystal.
# Compiled static to include openmp libraries for increased
# portability.
# Use edit omp variable to access openmp version
# 0=scalar
# 1=openmp
# ************************************************************

# Use "icpc" to compile source files
CCIN = g++
# Use "icpc" to link object files
LDIN = g++
# Compiler flags
# don't worry about MSVC++ pragmas, or '\r\n' characters
# the other disabled warnings require some cleanup
DISABLED_WARNINGS = -Wno-unknown-pragmas -Wno-multichar -Wno-write-strings -Wno-maybe-uninitialized -Wno-unused-but-set-variable

OMPCFLAGS = -O3 -g  -fopenmp -Dunix -Wall $(DISABLED_WARNINGS)
# Linker flags
OMPLDFLAGS = -O3 -g -fopenmp -Wall -Wno-unknown-pragmas $(DISABLED_WARNINGS)
# -openmp-profile


# ************************************************************
# Remove files
RM = rm -f
# ************************************************************

# ************************************************************
# Declare the main object files to link
# Only change these when core files are
# added (or renamed) to the LISFLOOD model
MAINOBJS = lisflood.o
LIBOBJS = pars.o \
	input.o \
	iterateq.o \
	boundary.o \
	fp_flow.o \
	fp_acc.o \
	fp_trent.o \
	ch_flow.o \
	util.o \
	output.o \
	chkpnt.o \
	por_flow.o \
	weir_flow.o \
	sgc.o \
	infevap.o \
	initialize.o \
	finalize.o \
	update.o \
	lib_bmi.o

# ************************************************************
# Program executable filename
EXE = lisflood_lin
LIB = liblisflood.so

# ************************************************************
# Openmp switch (if=0 then scalar; if=1 then openmp)
omp = 1

# Top-level rule to compile everything
all: $(EXE)

ifeq ($(omp)$,0)
# Linkage rule
$(EXE): $(MAINOBJS) # $(LIB)
	$(LD) $(LDFLAGS) $(MAINOBJS) -L. -llisflood -o $(EXE)

$(LIB): $(LIBOBJS)
	$(LD) $(LDFLAGS) $(LIBOBJS) -fPIC -shared -o $(LIB)
# Compile main object files
$(MAINOBJS) $(LIBOBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
endif

ifeq ($(omp)$,1)
# Linkage rule
$(EXE): $(MAINOBJS)  $(LIB)
	$(LDIN) $(OMPLDFLAGS) $(MAINOBJS) -L . -llisflood -o $(EXE)

$(LIB): $(LIBOBJS)
	$(LD) $(OMPLDFLAGS) $(LIBOBJS) -fPIC -shared -o $(LIB)

# Compile main object files
$(MAINOBJS): %.o: %.cpp
	$(CCIN) $(OMPCFLAGS) -c $< -o $@
$(LIBOBJS): %.o: %.cpp
	$(CCIN) $(OMPCFLAGS) -fPIC -c $< -o $@
endif

# ************************************************************
# Remove object files after compilation and linking
clean:
	$(RM) $(MAINOBJS) $(LIBOBJS) $(EXE) $(LIB)
