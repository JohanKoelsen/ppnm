CC = gcc

LDLIBS = -lm

default: out.txt
	cat out.txt
out.txt: hello
	./hello >out.txt
hello: hello.o
	$(CC) $(CFLAGS) -o hello hello.o
hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c
clean:
	$(RM) hello.o hello out.txt
text:
	echo $(LDLIBS)
	ech $(CC)
