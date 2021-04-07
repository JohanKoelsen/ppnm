CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: interpolation
	./interpolation > out.txt

interpolation.o: interpolation.c
	$(CC) $(CFLAGS) -c $< -o $@


interpolation: interpolation.o linear_spline.o

clean:
	$(RM) *.o out* interpolation
