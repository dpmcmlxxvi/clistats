# Simple makefile for clistats
CPP=g++
CPPFLAGS=-O2 -std=c++17 -Wall -Wextra -Wconversion
RM=/bin/rm
RMFLAGS=-f

# Compile and build clistats
all: clistats clistatstest

# Compile and build clistats
clistats: src/clistats.cpp
	$(CPP) $(CPPFLAGS) src/clistats.cpp -o $@

# Compile and build clistatstest
clistatstest: test/clistatstest.cpp
	$(CPP) $(CPPFLAGS) test/clistatstest.cpp -o $@

# Removed clistats executable
clean:
	$(RM) $(RMFLAGS) clistats clistatstest

# Run static analysis
analysis:
	cppcheck --enable=warning,performance,portability --std=c++17 src/clistats.cpp test/clistatstest.cpp

.PHONY: analysis
