CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default

default: exp.png exp_wrong.png rapport.pdf
rapport.pdf: main.tex exp.png exp_wrong.png
	pdflatex $<
	pdflatex $<

exp.png: data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "exp(x)";\
	set title "Exp(x)";\
	plot \
	"$<" using 1:2 with line title "exp(x)"\
	,"$<" using 1:3 with line title "exp from math.h"\
	'| pyxplot

exp_wrong.png: data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set tics out;\
	set xlabel "x";\
	set ylabel "exp(x)";\
	set title "Exp(x) without conditions";\
	plot \
	"$<" using 1:4 with line title "exp(x)"\
	,"$<" using 1:3 with line title "exp from math.h"\
	'| pyxplot



data.txt: Latex
	./$< > $@

Latex: Latex.o
.PHONEY: clean
clean:
	$(RM) *.o Latex *.png data.txt
