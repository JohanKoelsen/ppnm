CC = gcc
CFLAGS =-std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math.o: math.c
	$(CC) $(CFLAGS) -c $< -o $@


math: math.o
	$(CC) $(LDFLAGS) math.o -o math $(LDLIBS)

clean:
	$(RM) *.o out* math
