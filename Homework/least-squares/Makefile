CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: out.txt fit.png
	cat out.txt
fit.png: exp_vals.txt linfit.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "time";\
	set ylabel "activity";\
	set title "exp values and fit of 224Ra decay";\
	plot \
	"exp_vals.txt" using 1:2:3 with yerrorbars title "exp values with errors"\
	,"linfit.txt" using 1:2 with line linetype 1 color black title "linear fit"\
	,"linfit.txt" using 1:3 with line linetype 2 color red title "upper uncertainties"\
	,"linfit.txt" using 1:4 with line linetype 2 color red title "lower uncertainties"\
	'| pyxplot

out.txt: main
	./main > out.txt


main: main.o lsfit.o gs.o

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main.c: lsfit.c

lsfit.c: gs.c


clean:
	$(RM) *.o out* main
