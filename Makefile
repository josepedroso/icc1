DEPS = utils.h

CC = gcc

CFLAGS = -Wall -I/usr/include

LFLAGS = -lm -L/usr/lib 

all: invmat

invmat: main.o matriz.o LU.o utils.o dados.o
	$(CC) -o invmat main.o matriz.o LU.o utils.o dados.o $(LFLAGS)

matriz.c: matriz.h
LU.c: LU.h
utils.c: utils.h
dados.c: dados.h

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

clean:
	-rm -f *~ *\# \#* *.o invmat

purge: clean
	-rm -f invmat
