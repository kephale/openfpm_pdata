include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH) -I/Users/Antonio/HDF5/include 

agsim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -L/usr/local/lib -lparmetis -lmetis

all: agsim

.PHONY: clean all

clean:
	rm -f *.o *~ core agsim

