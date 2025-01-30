# Define the compiler
#CC = nvc
CC = gcc

# Source files for each executable
EXEC_SRCS = shallow_swap.c shallow_swap.acc.c shallow_swap.acc.Tile.c

# Common source files used by all executeables
COMMON_SRCS = wtime.c 

CFLAGS = -O2
#CFLAGS = -Wall -g

ifdef GPU
    CFLAGS += gpu=cc70 -Minfo=accel -Mnofma
endif

# Define the executable names based on source files
EXECS = $(basename $(EXEC_SRCS))

# Define the libraries to link
LDLIBS = -lm

# Default target
all: $(EXECS)

# Implicit rule to compile the executables
$(EXECS): $(COMMON_SRCS)

# Rule to clean the project
clean:
	rm -f $(EXECS) *.o
