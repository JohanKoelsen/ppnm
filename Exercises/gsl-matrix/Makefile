CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

default: out.txt
	cat out.txt

out.txt: matrix
	./matrix > out.txt

matrix.o: matrix.c
	$(CC) $(CFLAGS) -c $< -o $@


matrix: matrix.o
	$(CC) $(LDFLAGS) matrix.o -o matrix $(LDLIBS)

clean:
	$(RM) *.o out* matrix

