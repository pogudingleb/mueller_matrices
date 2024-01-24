EIGEN_PATH = ./eigen
CXX = g++
CXXFLAGS = -std=c++11 -pthread -O2 -I $(EIGEN_PATH)
SRC_DIR = src
SOURCES = $(SRC_DIR)/parsing.cpp $(SRC_DIR)/checkers.cpp $(SRC_DIR)/runners.cpp
OBJ_FILES = $(SOURCES:.cpp=.o)
BINARIES = matrices test

all: $(BINARIES)

matrices: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/main.cpp $(CXXFLAGS) -o matrices

test: $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(SRC_DIR)/test.cpp $(CXXFLAGS) -o test

clean:
	rm -f $(BINARIES)
	rm -f $(SRC_DIR)/*.o
