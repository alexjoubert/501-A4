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
	$(CC) $(PFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(PFLAGS) $(LFLAGS) $(IFLAGS) $(FLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)
