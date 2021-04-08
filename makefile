# FLAGLET package
# Copyright (C) 2012
# Boris Leistedt & Jason McEwen
# ======================================== #

# Directory for FLAG (required)
FLAGDIR = ${FLAG}
# Directory for S2LET (required)
S2LETDIR = ${S2LET}
# Directory for SO3 (required)
SO3DIR = ${SO3}
# Directory for SSHT (required)
SSHTDIR	= ${SSHT}
# Directory for FFTW (required)
FFTWDIR	= ${FFTW}
# Directory for MATLAB (optional)
MLAB	=  ${MATLAB}
# Directory for DOXYGEN (optional)
DOXYGEN_PATH = doxygen

UNAME 	:= $(shell uname)

# Compilers and options for C
CC	= gcc
OPT	= -Wall -g -fPIC -fopenmp -DFLAGLET_VERSION=\"0.0.1\" -DFLAGLET_BUILD=\"`git rev-parse HEAD`\"

# Config for dynamic library
ifeq ($(UNAME), Linux)
  DYLIBEXT = so
  DYLIBCMD = cc -flat_namespace -undefined suppress
endif
ifeq ($(UNAME), Darwin)
  DYLIBEXT = dylib
  DYLIBCMD = g++ -flat_namespace -dynamiclib -undefined suppress
endif

# ======================================== #

# === MATLAB ===
MLABINC	= ${MLAB}/extern/include
MLABLIB	= ${MLAB}/extern/lib
# --------------------
ifeq ($(UNAME), Linux)
  MEXEXT	= mexa64
endif
ifeq ($(UNAME), Darwin)
  MEXEXT	= mexmaci64
endif
# --------------------
MEX 		= ${MLAB}/bin/mex
MEXFLAGS	= -cxx

# === FLAGLET ===
FLAGLETDIR = .
FLAGLETLIB = $(FLAGLETDIR)/lib
FLAGLETINC = $(FLAGLETDIR)/include
FLAGLETBIN = $(FLAGLETDIR)/bin
FLAGLETLIBNM= flaglet
FLAGLETSRC = $(FLAGLETDIR)/src/main/c
FLAGLETOBJ = $(FLAGLETSRC)
FLAGLETTESTSRC = $(FLAGLETDIR)/src/test/c
FLAGLETTESTOBJ = $(FLAGLETTESTSRC)

# === FLAG ===
FLAGLIB   = $(FLAGDIR)/lib/c
FLAGINC   = $(FLAGDIR)/include/c
FLAGLIBNM = flag

# === S2LET ===
S2LETLIB = $(S2LETDIR)/lib/c
S2LETINC = $(S2LETDIR)/include/c
S2LETLIBNM = s2let

# === SO3 ===
SO3LIB = $(SO3DIR)/lib/c
SO3INC = $(SO3DIR)/include/c
SO3LIBNM = so3

# === SSHT ===
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBNM = ssht

# === FFTW ===
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3
# FFTWOMPLIBNM = fftw3_omp

# ======================================== #

FLAGLETSRCMAT	= $(FLAGLETDIR)/src/main/matlab
FLAGLETOBJMAT  = $(FLAGLETSRCMAT)
FLAGLETOBJMEX  = $(FLAGLETSRCMAT)

vpath %.c $(FLAGLETSRC)
vpath %.c $(FLAGLETTESTSRC)
vpath %.h $(FLAGLETINC)
vpath %_mex.c $(FLAGLETSRCMAT)

LDFLAGS = -L$(FLAGLETLIB) -l$(FLAGLETLIBNM) -L$(FLAGLIB) -l$(FLAGLIBNM) -L$(S2LETLIB) -l$(S2LETLIBNM) -L$(SO3LIB) -l$(SO3LIBNM) -L$(FFTWLIB) -l$(FFTWLIBNMM) -L$(SSHTLIB) -l$(SSHTLIBNM) -lm -lc

LDFLAGSMEX = -L$(FLAGLETLIB) -l$(FLAGLETLIBNM) -L$(FLAGLIB) -l$(FLAGLIBNM) -L$(S2LETLIB) -l$(S2LETLIBNM) -L$(SO3LIB) -l$(SO3LIBNM) -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNMM) -L$(SSHTLIB) -l$(SSHTLIBNM)

FFLAGS  = -I$(FLAGLETINC) -I$(FFTWINC) -I$(SSHTINC) -I$(FLAGINC) -I$(S2LETINC) -I$(SO3INC)

FLAGLETOBJS = $(FLAGLETOBJ)/flaglet_transform.o \
			  $(FLAGLETOBJ)/flaglet_tiling.o \
			  $(FLAGLETOBJ)/flaglet_axisym.o

FLAGLETOBJSMAT = $(FLAGLETOBJMAT)/flaglet_tiling_mex.o \
				 $(FLAGLETOBJMAT)/flaglet_analysis_mex.o \
				 $(FLAGLETOBJMAT)/flaglet_synthesis_mex.o

FLAGLETOBJSMEX = $(FLAGLETOBJMEX)/flaglet_tiling_mex.$(MEXEXT) \
				 $(FLAGLETOBJMAT)/flaglet_analysis_mex.$(MEXEXT) \
				 $(FLAGLETOBJMAT)/flaglet_synthesis_mex.$(MEXEXT)

$(FLAGLETOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGLETTESTOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGLETOBJMAT)/%_mex.o: %_mex.c $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC}

$(FLAGLETOBJMEX)/%_mex.$(MEXEXT): $(FLAGLETOBJMAT)/%_mex.o $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a
	$(MEX) $< -output $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test tidy

.PHONY: matlab
matlab: $(FLAGLETOBJSMEX)

.PHONY: all
all: lib matlab test about tidy

.PHONY: lib
lib: $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a
$(FLAGLETLIB)/lib$(FLAGLETLIBNM).a: $(FLAGLETOBJS)
	ar -r $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a $(FLAGLETOBJS)

.PHONY: test
test: lib $(FLAGLETBIN)/flaglet_test
$(FLAGLETBIN)/flaglet_test: $(FLAGLETTESTOBJ)/flaglet_test.o $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a
	$(CC) $(OPT) $< -o $(FLAGLETBIN)/flaglet_test $(LDFLAGS)

.PHONY: about
about: $(FLAGLETBIN)/flaglet_about
$(FLAGLETBIN)/flaglet_about: $(FLAGLETOBJ)/flaglet_about.o
	$(CC) $(OPT) $< -o $(FLAGLETBIN)/flaglet_about

.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(FLAGLETDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(FLAGLETDIR)/doc/c/*

.PHONY: clean
clean:	tidy
	rm -f $(FLAGLETLIB)/lib$(FLAGLETLIBNM).a
	rm -f $(FLAGLETOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(FLAGLETBIN)/flaglet_test
	rm -f $(FLAGLETBIN)/flaglet_about

.PHONY: tidy
tidy:
	rm -f $(FLAGLETOBJ)/*.o
	rm -f $(FLAGLETTESTOBJ)/*.o
	rm -f $(FLAGLETOBJMEX)/*.o
	rm -f *~
