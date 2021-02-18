CC = gcc
CFLAGS =-std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon epsilon_tau
	./epsilon > out.txt
	./epsilon_tau >> out.txt

epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c $< -o $@

epsilon_tau.o: epsilon_tau.c
	$(CC) $(CFLAGS) -c $< -o $@
epsilon_tau: epsilon_tau.o
	$(CC) $(LDFLAGS) epsilon_tau.o -o epsilon_tau $(LDLIBS)


epsilon: epsilon.o
	$(CC) $(LDFLAGS) epsilon.o -o epsilon $(LDLIBS)

clean:
	$(RM) *.o out* epsilon
