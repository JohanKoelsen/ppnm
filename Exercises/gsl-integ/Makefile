CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@


main: main.o
	$(CC) $(LDFLAGS) main.o -o main $(LDLIBS)

clean:
	$(RM) *.o out* mainy

