CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o qnewton.o downhill.o
main.c: qnewton.c
downhill.o: downhill.c

.PHONEY: clean


clean:
	$(RM) *.o out* main
