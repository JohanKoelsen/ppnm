CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS  += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: dummy
	./dummy > out.txt

dummy.o: dummy.c
	$(CC) $(CFLAGS) -c $< -o $@

dummy: dummy.o
	$(CC) $(LDFLAGS) dummy.o -o dummy $(LDLIBS)

clean:
	$(RM) *.o out* dummy
