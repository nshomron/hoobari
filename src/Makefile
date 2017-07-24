all:phred

phred:phred.c
	gcc -c -Wall -fpic phred.c -lm    
	gcc -shared -o libphred.so phred.o

clean:
	rm phred.o libphred.so
