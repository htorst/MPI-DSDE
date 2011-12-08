MPICXX=mpicxx
CFLAGS=-g 

all: libdsde.a main

libdsde.a: dsde.o dsde_exchange_alltoall.o dsde_exchange_reduce_scatter.o
	ar r $@ dsde.o dsde_exchange_alltoall.o dsde_exchange_reduce_scatter.o
	ranlib $@

%.o: %.cpp *.h
	$(MPICXX) $(CFLAGS) -c $< -o $@ -I.

main: main.cpp dsde.o *h libdsde.a
	$(MPICXX) $(CFLAGS) $< libdsde.a -o $@

clean:
	rm -f *.o main
