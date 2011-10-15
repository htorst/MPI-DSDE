MPICXX=mpicxx

all: dsde.o

dsde.o: dsde.cpp dsde.h
	$(MPICXX) -c dsde.cpp -o dsde.o -I.

clean:
	rm -f dsde.o
