CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt plot.png
	cat out.txt

plot.png: plots.pgi points.txt functions.txt
	pyxplot $<


out.txt: main
	./main > out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o ann.o qnewton.o
main.c: ann.c
ann.c: qnewton.c

.PHONEY: clean


clean:
	$(RM) *.o out* main
