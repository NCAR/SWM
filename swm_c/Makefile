# Define the compiler
CC = nvc
#CC = gcc

# Source files for each executable
EXEC_SRCS = shallow_swap.c shallow_swap.acc.c shallow_swap.acc.Tile.c

# Common source files used by all executables
COMMON_SRCS = wtime.c 

CFLAGS = -O2
#CFLAGS = -Wall -g

ifdef GPU
    CFLAGS += gpu=cc70 -Minfo=accel -Mnofma
endif

# Define the executable names based on source files
EXECS = $(basename $(EXEC_SRCS))

# Libraries to link to
LDLIBS = -lm

# Default target
all: $(EXECS)

# Implicit rule to compile the executables
$(EXECS): $(COMMON_SRCS)

clean:
	rm -f $(EXECS) *.o
