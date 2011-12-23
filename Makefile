LIBNBC_HOME=$(HOME)/root
MPICXX=mpicxx
CFLAGS=-g -O0 -I$(LIBNBC_HOME)/include
#LDFLAGS=-L$(LIBNBC_HOME)/lib -lnbc
LDFLAGS=

#  dsde_exchange_ibarrier.o
OBJS= \
  dsde.o \
  dsde_exchange_brucks_inline.o \
  dsde_exchangev_alltoall.o \
  dsde_exchangev_reduce_scatter.o \
  dsde_exchangev_accumulate.o \
  dsde_exchangev_brucks.o \
  dsde_reduce_scatter_block_brucks.o \
  dsde_reduce_scatter_block_hbrucks.o

all: libdsde.a main

libdsde.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

%.o: %.cpp *.h
	$(MPICXX) $(CFLAGS) -c $< -o $@ -I.

main: main.cpp dsde.o *h libdsde.a
	$(MPICXX) $(CFLAGS) $< libdsde.a -o $@ $(LDFLAGS)

clean:
	rm -f *.o main libdsde.*
