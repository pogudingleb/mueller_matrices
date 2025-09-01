EIGEN_PATH = ./eigen
CXX = g++
CXXFLAGS = -std=c++11 -pthread -O2 -I $(EIGEN_PATH)
SRC_DIR = src
SOURCES = $(SRC_DIR)/parsing.cpp $(SRC_DIR)/checkers.cpp $(SRC_DIR)/runners.cpp $(SRC_DIR)/ferrari.cpp
OBJ_FILES = $(SOURCES:.cpp=.o)
BINARIES = matrices test test_ferrari eigenvalues

all: $(BINARIES)

matrices: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/main.cpp $(CXXFLAGS) -o matrices

test: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/test.cpp $(CXXFLAGS) -o test

test_ferrari: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/test_ferrari.cpp $(CXXFLAGS) -o test_ferrari

eigenvalues: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/eigenvalues.cpp $(CXXFLAGS) -o eigenvalues

clean:
	rm -f $(BINARIES)
	rm -f $(SRC_DIR)/*.o
