CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	./$<

main: main.o rkstep.o

rkstep.o: rkstep.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONEY:clean
clean:
	$(RM) main *.txt *.o *.png
