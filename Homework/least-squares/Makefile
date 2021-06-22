CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt fit.png
	cat out.txt

fit.png: fit.pgi exp_vals.txt linfit.txt
	pyxplot $<

out.txt: main
	./main > out.txt


main: main.o lsfit.o gs.o

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main.c: lsfit.c

lsfit.c: gs.c


clean:
	$(RM) *.o out* main
