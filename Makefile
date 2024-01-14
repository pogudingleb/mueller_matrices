EIGEN_PATH = ./eigen
CFLAGS = -std=c++11 -O2 -I $(EIGEN_PATH)
OBJ_FILES = parsing.o checkers.o runners.o
BINARIES = matrices test

all: $(BINARIES)

parsing.o:
	g++ parsing.cpp $(CFLAGS) -c -o parsing.o

checkers.o:
	g++ checkers.cpp $(CFLAGS) -c -o checkers.o

runners.o:
	g++ runners.cpp $(CFLAGS) -c -o runners.o

matrices: $(OBJ_FILES)
	g++ $(OBJ_FILES) main.cpp $(CFLAGS) -o matrices

test: $(OBJ_FILES)
	g++ $(OBJ_FILES) test.cpp $(CFLAGS) -o test

clean:
	rm -f $(BINARIES)
	rm -f *.o
