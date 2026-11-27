CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wno-sign-compare

all: bin/vcsolve bin/gen bin/test_bruteforce

bin/vcsolve: src/main.cpp src/*.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

bin/gen: src/gen.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/gen.cpp

bin/test_bruteforce: tests/test_bruteforce.cpp src/*.hpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ tests/test_bruteforce.cpp

test: bin/test_bruteforce
	./bin/test_bruteforce

clean:
	rm -rf bin

.PHONY: all test clean
