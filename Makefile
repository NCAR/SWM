CC = nvc
CFLAGS = -O2
ifdef GPU
    CFLAGS += gpu=cc70 -Minfo=accel -Mnofma
endif
SRC = shallow_swap.acc.c wtime.c
OBJ = $(SRC:.c=.o)
EXEC = SWM_acc

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(EXEC)