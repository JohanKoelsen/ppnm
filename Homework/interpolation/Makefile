CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


.PHONEY: default

default: out.txt linplot.png integplot.png qsplineplot.png
	cat out.txt

qsplineplot.png: qsplineplot.pgi functions.txt qspline.txt
	pyxplot $<

integplot.png: integplot.pgi integration_points.txt
	pyxplot $<

linplot.png: linplot.pgi points.txt functions.txt
	pyxplot $<

out.txt: main
	./main > out.txt

main: main.o interpolation.o

interpolation.o: interpolation.c
	$(CC) $(CFLAGS) -c -o $@ $<



.PHONEY:clean
clean:
	$(RM) main *.txt *.o *.png
