CC=g++

CFLAGS=-c
PFLAGS=-pg
LFLAGS=-L ../lib/lib
IFLAGS=-I ../lib/include
LIBS=../lib/lib/libsndfile.a
FLAGS=-static
SOURCES=convolveFFT.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=convolve

all: convolve

$(EXECUTABLE): $(OBJECTS)
	$(CC) -O2 $(PFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) -O2 $(CFLAGS) $(PFLAGS) $(LFLAGS) $(IFLAGS) $(FLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)
