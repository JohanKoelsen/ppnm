CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o Jd.o
main.c: Jd.c
clean:
	$(RM) *.o out* main
