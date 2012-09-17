# ======================================== #

# Directory for FLAG
FLAGDIR	= ${FLAG}
# Directory for S2LET
S2LETDIR = ${S2LET}
# Directory for SSHT
SSHTDIR	= ${SSHT}
# Directory for FFTW
FFTWDIR	= ${FFTW}
# Directory for GSL
GSLDIR	= ${GSL}
# Directory for MATLAB
MLAB	=  /Applications/MATLAB_R2011b.app
# Directory for DOXYGEN
DOXYGEN_PATH=/Applications/Doxygen.app/Contents/Resources/doxygen

# Compiler and options
CC	= gcc
OPT	= -Wall -O3 -g -DFLAGLET_VERSION=\"1.0\" -DFLAGLET_BUILD=\"`svnversion -n .`\"
UNAME := $(shell uname)

# ======================================== #

# === MATLAB ===
ifeq ($(UNAME), Linux)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexa64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif
ifeq ($(UNAME), Darwin)
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib
  MEXEXT	= mexmaci64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif

# === FLAGLET ===
FLAGLETDIR = .
FLAGLETLIB = $(FLAGLETDIR)/lib
FLAGLETINC = $(FLAGLETDIR)/include
FLAGLETBIN = $(FLAGLETDIR)/bin
FLAGLETLIBN= flaglet
FLAGLETSRC = $(FLAGLETDIR)/src/main/c
FLAGLETOBJ = $(FLAGLETSRC)

# === S2LET ===
S2LETLIB = $(S2LETDIR)/lib
S2LETINC = $(S2LETDIR)/include
S2LETLIBN= s2let

# === FLAG ===
FLAGLIB = $(FLAGDIR)/lib
FLAGINC = $(FLAGDIR)/include
FLAGLIBN= flag

# === SSHT ===
SSHTLIB	= $(SSHTDIR)/lib/c
SSHTINC	= $(SSHTDIR)/include/c
SSHTLIBN= ssht

# === FFTW ===
FFTWINC	    = $(FFTWDIR)/include
FFTWLIB     = $(FFTWDIR)/lib
FFTWLIBNM   = fftw3

# ======================================== #

FLAGLETSRCMAT	= $(FLAGLETDIR)/src/main/matlab
FLAGLETOBJMAT  = $(FLAGLETSRCMAT)
FLAGLETOBJMEX  = $(FLAGLETSRCMAT)

vpath %.c $(FLAGLETSRC)
vpath %.h $(FLAGLETINC)
vpath %_mex.c $(FLAGLETSRCMAT)

LDFLAGS = -L$(FLAGLETLIB) -l$(FLAGLETLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -lm -lc

LDFLAGSMEX = -L$(FLAGLETLIB) -l$(FLAGLETLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN)

FFLAGS  = -I$(FLAGLETINC) -I$(FFTWINC) -I$(SSHTINC) -I$(FLAGINC) -I$(S2LETINC) 

FLAGLETOBJS= $(FLAGLETOBJ)/flaglet_axisym.o	\
	$(FLAGLETOBJ)/flaglet_tiling.o

FLAGLETOBJSMAT = $(FLAGLETOBJMAT)/flaglet_axisym_tiling_mex.o	\
	$(FLAGLETOBJMAT)/flaglet_axisym_analysis_mex.o	\
	$(FLAGLETOBJMAT)/flaglet_axisym_synthesis_mex.o

FLAGLETOBJSMEX = $(FLAGLETOBJMEX)/flaglet_axisym_tiling_mex.$(MEXEXT)	\
	$(FLAGLETOBJMEX)/flaglet_axisym_analysis_mex.$(MEXEXT)	\
	$(FLAGLETOBJMEX)/flaglet_axisym_synthesis_mex.$(MEXEXT)

$(FLAGLETOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(FLAGLETOBJMAT)/%_mex.o: %_mex.c $(FLAGLETLIB)/lib$(FLAGLETLIBN).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(FLAGLETOBJMEX)/%_mex.$(MEXEXT): $(FLAGLETOBJMAT)/%_mex.o $(FLAGLETLIB)/lib$(FLAGLETLIBN).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test about tidy

.PHONY: matlab
matlab: $(FLAGLETOBJSMEX)

.PHONY: all
all: lib matlab doc test about tidy

.PHONY: lib
lib: $(FLAGLETLIB)/lib$(FLAGLETLIBN).a
$(FLAGLETLIB)/lib$(FLAGLETLIBN).a: $(FLAGLETOBJS)
	ar -r $(FLAGLETLIB)/lib$(FLAGLETLIBN).a $(FLAGLETOBJS)

.PHONY: test
test: lib $(FLAGLETBIN)/flaglet_test
$(FLAGLETBIN)/flaglet_test: $(FLAGLETOBJ)/flaglet_test.o $(FLAGLETLIB)/lib$(FLAGLETLIBN).a
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
clean:	tidy cleandoc
	rm -f $(FLAGLETLIB)/lib$(FLAGLETLIBN).a
	rm -f $(FLAGLETOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(FLAGLETBIN)/flaglet_test
	rm -f $(FLAGLETBIN)/flaglet_about

.PHONY: tidy
tidy:
	rm -f $(FLAGLETOBJ)/*.o
	rm -f $(FLAGLETOBJMEX)/*.o
	rm -f *~ 

