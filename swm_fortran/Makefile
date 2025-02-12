# Define the compiler
FC = gfortran

# Source files for each executable
EXEC_SRCS = swm_fortran.F90

# Common source files used by all executables
COMMON_SRCS =

FFLAGS = -O2 -fdefault-real-8# -cpp
#FFLAGS = -Wall -g

ifdef GPU
    FFLAGS += gpu=cc70 -Minfo=accel -Mnofma
endif

# Define the executable names based on source files
EXECS = $(basename $(EXEC_SRCS))

# Libraries to link to
LDLIBS = -lm

# Default target
all: $(EXECS)

# Implicit rule to compile the executables
$(EXECS):
	$(FC) $(FFLAGS) -o $(EXECS) $(EXEC_SRCS)

clean:
	rm -f $(EXECS) *.o