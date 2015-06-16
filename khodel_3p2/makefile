CC=g++
CFLAGS=-c -Wall -Wextra -pedantic -std=c++0x -fopenmp -I ../lib/armadillo_5_1/include -I ../lib  -DARMA_DONT_USE_WRAPPER
LIB=-lgsl -lgslcblas -lm -llapack -lblas -larmadillo -fopenmp 
SOURCES:=$(wildcard *.cpp) $(wildcard ../lib/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
OUT=a.out 

all:$(SOURCES) $(OUT)


$(OUT):$(OBJECTS)
	$(CC)  $(LIB) $(OBJECTS) -o $@

%.o:%.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:$(OUT)
	rm -rf *.o ../lib/*.o

.PHONY:test

