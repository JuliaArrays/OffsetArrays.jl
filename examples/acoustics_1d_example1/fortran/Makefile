
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  

# ---------------------------------
# List of sources for this program:
# ---------------------------------

MODULES = \

SOURCES = \
  qinit.f90 \
  setprob.f90 \
  ./rp1_acoustics.f90 \
  ./setaux.f90 \
  ./bc1.f \
  ./b4step1.f90 \
  ./driver.f90 \
  ./claw1ez.f \
  ./claw1.f \
  ./copyq1.f \
  ./inlinelimiter.f90 \
  ./opendatafile.f \
  ./out1.f \
  ./src1.f90 \
  ./step1.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

