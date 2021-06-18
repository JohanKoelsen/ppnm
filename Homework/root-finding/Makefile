  CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

default: out.txt Hydrogen.png
	#cat out.txt

out.txt: main
	./main #> out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o roots.o gs.o rkstep.o
main.c: roots.c
roots.c: gs.c rkstep.c
clean:
	$(RM) *.o out* main
