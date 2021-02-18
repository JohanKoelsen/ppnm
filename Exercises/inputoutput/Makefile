CC = gcc
CFLAGS =-std=gnu99 -O -Wall
LDFLAGS =
LDLIBS = -lm

default: cmdline stdios loadtext
	touch out.txt
	rm out.txt
	./cmdline $(param1) $(param2) > out.txt
	cat out.txt



cmdline.o: cmdline.c
	$(CC) $(CFLAGS) -c $< -o $@
stdios.o: stdios.c
	$(CC) $(CFLAGS) -c $< -o $@

loadtext.o: loadtext.c
	$(CC) $(CFLAGS) -c $< -o $@

cmdline: cmdline.o
	$(CC) $(LDFLAGS) cmdline.o -o cmdline $(LDLIBS)

stdios: stdios.o
	$(CC) $(LDFLAGS) stdios.o -o stdios $(LDLIBS)

loadtext: loadtext.o
	$(CC) $(LDFLAGS) loadtext.o -o loadtext $(LDLIBS)

clean:
	$(RM) *.o out* cmdline stdios
