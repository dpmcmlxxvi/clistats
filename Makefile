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
	$(RM) $(RMFLAGS) *.gcda *.gcno coverage.info
	$(RM) $(RMFLAGS) -r coverage_html

# Run static analysis
analysis:
	cppcheck --enable=warning,performance,portability --std=c++17 src/clistats.cpp test/clistatstest.cpp

# Build and run test suite with coverage instrumentation, then generate an lcov report
coverage:
	$(CPP) -O0 --coverage -std=c++17 -D_CLISTATS_TESTING test/clistatstest.cpp -o clistatstest
	./clistatstest
	lcov --capture --directory . --output-file coverage.info
	lcov --remove coverage.info '/usr/*' --output-file coverage.info
	genhtml coverage.info --output-directory coverage_html
	xdg-open coverage_html/index.html

.PHONY: analysis coverage
