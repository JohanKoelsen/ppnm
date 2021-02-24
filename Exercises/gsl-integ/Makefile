CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

default: out.txt
	cat out.txt

out.txt: gsl-integ
	./gsl-integ > out.txt

gsl-integ.o: gsl-integ.c
	$(CC) $(CFLAGS) -c $< -o $@


gsl-integ: gsl-integ.o
	$(CC) $(LDFLAGS) gsl-integ.o -o gsl-integ $(LDLIBS)

clean:
	$(RM) *.o out* gsl-integ
