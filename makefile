
all:

#module load gcc/5.2.0 

#gcc -E main.c -o main.o
#	gcc -E drift.c -o drift.o
#	gcc -c main.c
#	gcc -c drift.c
#	gcc -o main.o drift.o -o pebble
	gcc -mcmodel=medium *.c -o pebble -lstdc++

clean:
	rm -rf *.o pebble
