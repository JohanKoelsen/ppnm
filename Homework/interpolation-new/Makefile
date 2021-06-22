CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


.PHONEY: default

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o interpolation.o

interpolation.o: interpolation.c
	$(CC) $(CFLAGS) -c -o $@ $<



.PHONEY:clean
clean:
	$(RM) main *.txt *.o *.png
