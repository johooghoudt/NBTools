.PHONY: all clean

CXX = g++
TARGETS = lib/librandlibc.a lib/libcpplapack_plus.a lib/libNBHmm.a include/cpplapack_plus.h include/nbhmm.h

all: $(TARGETS)

clean:
	$(RM) $(TARGETS) *.o

lib/librandlibc.a: 
	$(MAKE) -C distfiles/randlib/source/randlib.c/src/
	cp distfiles/randlib/source/randlib.c/src/librandlibc.a $@

lib/libcpplapack_plus.a: cpplapack_plus.o
	ar cr $@ $<
	ranlib $@

lib/libNBHmm.a: nbhmm.o
	ar cr $@ $<
	ranlib $@

cpplapack_plus.o: cpplapack_plus/cpplapack_plus.cpp
	$(CXX) -c -I include -I include/cpplapack $< -o $@

nbhmm.o: NBHmm/nbhmm.cpp include/cpplapack_plus.h
	$(CXX) -c -I include -I include/cpplapack $< -o $@

include/cpplapack_plus.h: cpplapack_plus/cpplapack_plus.h
	cp $< $@

include/nbhmm.h: NBHmm/nbhmm.h
	cp $< $@
