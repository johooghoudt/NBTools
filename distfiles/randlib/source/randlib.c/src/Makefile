CC = gcc
OBJS = linpack.o com.o randlib.o

.PHONY: all clean

all: librandlibc.a

librandlibc.a: $(OBJS)
	ar cr librandlibc.a $(OBJS)
	ranlib librandlibc.a

.c.o:
	$(CC) -O2 -c $<

clean:
	$(RM) librandlibc.a *.o
