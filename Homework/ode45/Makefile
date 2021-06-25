CC = gcc
CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


.PHONEY: default

default: out.txt harmonics.png sir1.png sir2.png sir3.png threebody.png
	cat out.txt

threebody.png: threebody.pgi threebody.txt
	pyxplot $<

harmonics.png: harmonics.pgi harmonics.txt
	pyxplot $<

sir1.png: sir1.pgi SIR1.txt
	pyxplot $<

sir2.png: sir2.pgi SIR2.txt
	pyxplot $<

sir3.png: sir3.pgi SIR3.txt
	pyxplot $<

out.txt: main
	./main > out.txt

main: main.o rkstep.o gravity.o

rkstep.o: rkstep.c
	$(CC) $(CFLAGS) -c -o $@ $<
gravity.o: gravity.c
	$(CC) $(CFLAGS) -c -o $@ $<


.PHONEY:clean
clean:
	$(RM) main *.txt *.o *.png
