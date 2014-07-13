
# Simple makefile for clistats
CPP=g++
CPPFLAGS=-O2
RM=/bin/rm
RMFLAGS=-rf

# Compile and build clistats
all: clistats

# Compile and build clistats
clistats: src/clistats.cpp
	$(CPP) $(CPPFLAGS) src/clistats.cpp -o clistats.exe

# Removed clistats executable
clean:
	$(RM) $(RMFLAGS) clistats.exe
