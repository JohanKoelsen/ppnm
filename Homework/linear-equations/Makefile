CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: QR_decomp
	./QR_decomp > out.txt

QR_decomp.o: QR_decomp.c
	$(CC) $(CFLAGS) -c $< -o $@

QR_decomp: QR_decomp.o gs.o
QR_decomp.c: gs.c
clean:
	$(RM) *.o out* GS_decomp QR_decomp
