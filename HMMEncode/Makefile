CXX = g++
CXXFLAGS += -I../include -I../include/cpplapack
LDFLAGS += -lblas -llapack
LIBS += ../lib/librandlibc.a ../lib/libcpplapack_plus.a ../lib/libNBHmm.a

TARGET = encode
.PHONY: clean check-syntax

encode: encode.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) *.o $(TARGET)

check-syntax:
	$(CXX) $(CXXFLAGS) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)


