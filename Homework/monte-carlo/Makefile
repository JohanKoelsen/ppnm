CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt error.png
	cat out.txt

error.png: error.pgi error.txt
	pyxplot $<

out.txt: main
	./main > out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o mc.o
main.c: mc.c
clean:
	$(RM) *.o out* main
