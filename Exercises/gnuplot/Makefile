CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: erf.pyxplot.png gam.pyxplot.png

erf.pyxplot.png: data.txt tabulated.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "myerf"\
	,"tabulated.data.txt" using 1:2 with points title "tabulated data"\
	'| pyxplot

gam.pyxplot.png: data.txt gam-tabulated.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:5 with line title "tgamma from math.h"\
	,"$<" using 1:6 with line title "gsl-sf-gamma"\
	,"$<" using 1:7 with line title "mygam"\
	,"gam-tabulated.data.txt" using 1:5 with points title "gam-tabulated data"\
	'| pyxplot
 
data.txt: main
	./$< > $@

tabulated.data.txt: main
	./$< > $@

gam-tabulated.data.txt: main
	./$< > $@

main: main.o myerf.o mygam.o
.PHONEY: clean
clean:
	$(RM) *.o main *.png data.txt
