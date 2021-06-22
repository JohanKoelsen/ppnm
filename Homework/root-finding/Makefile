  CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: out.txt hydrogen.png
	cat out.txt
out.txt: main
	./main > out.txt

hydrogen.png: hydrogen.pgi hydrogen.txt
	pyxplot $<


main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o roots.o gs.o rkstep.o
main.c: roots.c
roots.c: gs.c rkstep.c
clean:
	$(RM) *.o out* main
