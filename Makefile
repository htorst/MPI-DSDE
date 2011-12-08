MPICXX=mpicxx
CFLAGS=-g 

all: dsde.o main

dsde.o: dsde.cpp *.h
	$(MPICXX) $(CFLAGS) -c dsde.cpp -o dsde.o -I.

main: main.cpp dsde.o *h
	$(MPICXX) $(CFLAGS) $< dsde.o -o $@

clean:
	rm -f dsde.o
