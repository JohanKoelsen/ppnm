CC = gcc
default: out.txt
	cat out.txt
out.txt: hello
	./hello >out.txt
hello: hello.o
	$(CC) $(CFLAGS) -c hello.c
clean:
	$(RM) hello.0 hello out.txt
text:
	echo $(LDLIBS)
	ech $(CC)
