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
OPT	= -Wall -O3 -g -DB3LET_VERSION=\"1.0\" -DB3LET_BUILD=\"`svnversion -n .`\"
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

# === B3LET ===
B3LETDIR = .
B3LETLIB = $(B3LETDIR)/lib/c
B3LETINC = $(B3LETDIR)/include/c
B3LETBIN = $(B3LETDIR)/bin/c
B3LETLIBN= b3let
B3LETSRC = $(B3LETDIR)/src/c
B3LETOBJ = $(B3LETSRC)

# === S2LET ===
S2LETLIB = $(S2LETDIR)/lib/c
S2LETINC = $(S2LETDIR)/include/c
S2LETLIBN= s2let

# === FLAG ===
FLAGLIB = $(FLAGDIR)/lib/c
FLAGINC = $(FLAGDIR)/include/c
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

B3LETSRCMAT	= $(B3LETDIR)/src/matlab
B3LETOBJMAT  = $(B3LETSRCMAT)
B3LETOBJMEX  = $(B3LETSRCMAT)

vpath %.c $(B3LETSRC)
vpath %.h $(B3LETSRC)
vpath %_mex.c $(B3LETSRCMAT)

LDFLAGS = -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -L$(B3LETLIB) -l$(B3LETLIBN) -lm

LDFLAGSMEX = -I/usr/local/include -L$(FFTWLIB) -l$(FFTWLIBNM) -L$(SSHTLIB) -l$(SSHTLIBN) -L$(FLAGLIB) -l$(FLAGLIBN) -L$(S2LETLIB) -l$(S2LETLIBN) -L$(B3LETLIB) -l$(B3LETLIBN)

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC) -I$(FLAGINC) -I$(S2LETINC) -I$(B3LETINC)

B3LETOBJS= $(B3LETOBJ)/b3let_axisym.o	\
	$(B3LETOBJ)/b3let_tilling.o

B3LETOBJSMAT = $(B3LETOBJMAT)/b3let_axisym_tilling_mex.o	\
	$(B3LETOBJMAT)/b3let_axisym_analysis_mex.o

B3LETOBJSMEX = $(B3LETOBJMEX)/b3let_axisym_tilling_mex.$(MEXEXT)	\
	$(B3LETOBJMEX)/b3let_axisym_analysis_mex.$(MEXEXT)

$(B3LETOBJ)/%.o: %.c
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

$(B3LETOBJMAT)/%_mex.o: %_mex.c $(B3LETLIB)/lib$(B3LETLIBN).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC} 

$(B3LETOBJMEX)/%_mex.$(MEXEXT): $(B3LETOBJMAT)/%_mex.o $(B3LETLIB)/lib$(B3LETLIBN).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

# ======================================== #

.PHONY: default
default: lib test about tidy

.PHONY: matlab
matlab: $(B3LETOBJSMEX)

.PHONY: all
all: lib matlab doc test about tidy

.PHONY: lib
lib: $(B3LETLIB)/lib$(B3LETLIBN).a
$(B3LETLIB)/lib$(B3LETLIBN).a: $(B3LETOBJS)
	ar -r $(B3LETLIB)/lib$(B3LETLIBN).a $(B3LETOBJS)

.PHONY: test
lib: lib $(B3LETBIN)/b3let_test
$(B3LETBIN)/b3let_test: $(B3LETOBJ)/b3let_test.o $(B3LETLIB)/lib$(B3LETLIBN).a
	$(CC) $(OPT) $< -o $(B3LETBIN)/b3let_test $(LDFLAGS)
	$(B3LETBIN)/b3let_test

.PHONY: about
about: $(B3LETBIN)/b3let_about
$(B3LETBIN)/b3let_about: $(B3LETOBJ)/b3let_about.o 
	$(CC) $(OPT) $< -o $(B3LETBIN)/b3let_about
	$(B3LETBIN)/b3let_about

.PHONY: doc
doc:
	$(DOXYGEN_PATH) $(B3LETDIR)/src/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(B3LETDIR)/doc/html/*

.PHONY: clean
clean:	tidy cleandoc
	rm -f $(B3LETLIB)/lib$(B3LETLIBN).a
	rm -f $(B3LETOBJMEX)/*_mex.$(MEXEXT)
	rm -f $(B3LETBIN)/b3let_test
	rm -f $(B3LETBIN)/b3let_about

.PHONY: tidy
tidy:
	rm -f $(B3LETOBJ)/*.o
	rm -f *~ 

