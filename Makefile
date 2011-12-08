LIBNBC_HOME=$(HOME)/root
MPICXX=mpicxx
CFLAGS=-g -I$(LIBNBC_HOME)/include
LDFLAGS=-L$(LIBNBC_HOME)/lib -lnbc

all: libdsde.a main

libdsde.a: dsde.o dsde_exchange_alltoall.o dsde_exchange_reduce_scatter.o dsde_exchange_accumulate.o dsde_exchange_ibarrier.o
	ar r $@ dsde.o dsde_exchange_alltoall.o dsde_exchange_reduce_scatter.o dsde_exchange_accumulate.o dsde_exchange_ibarrier.o
	ranlib $@

%.o: %.cpp *.h
	$(MPICXX) $(CFLAGS) -c $< -o $@ -I.

main: main.cpp dsde.o *h libdsde.a
	$(MPICXX) $(CFLAGS) $< libdsde.a -o $@ $(LDFLAGS)

clean:
	rm -f *.o main
