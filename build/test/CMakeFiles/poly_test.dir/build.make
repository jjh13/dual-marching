# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/joshua/Projects/sisl

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/joshua/Projects/sisl/build

# Include any dependencies generated for this target.
include test/CMakeFiles/poly_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/poly_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/poly_test.dir/flags.make

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o: test/CMakeFiles/poly_test.dir/flags.make
test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o: ../test/polynomial_test.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/joshua/Projects/sisl/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o"
	cd /Users/joshua/Projects/sisl/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/poly_test.dir/polynomial_test.cpp.o -c /Users/joshua/Projects/sisl/test/polynomial_test.cpp

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poly_test.dir/polynomial_test.cpp.i"
	cd /Users/joshua/Projects/sisl/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/joshua/Projects/sisl/test/polynomial_test.cpp > CMakeFiles/poly_test.dir/polynomial_test.cpp.i

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poly_test.dir/polynomial_test.cpp.s"
	cd /Users/joshua/Projects/sisl/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/joshua/Projects/sisl/test/polynomial_test.cpp -o CMakeFiles/poly_test.dir/polynomial_test.cpp.s

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.requires:
.PHONY : test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.requires

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.provides: test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/poly_test.dir/build.make test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.provides.build
.PHONY : test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.provides

test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.provides.build: test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o

# Object files for target poly_test
poly_test_OBJECTS = \
"CMakeFiles/poly_test.dir/polynomial_test.cpp.o"

# External object files for target poly_test
poly_test_EXTERNAL_OBJECTS =

test/poly_test: test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o
test/poly_test: test/CMakeFiles/poly_test.dir/build.make
test/poly_test: /usr/local/lib/libfftw3.dylib
test/poly_test: test/CMakeFiles/poly_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable poly_test"
	cd /Users/joshua/Projects/sisl/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poly_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/poly_test.dir/build: test/poly_test
.PHONY : test/CMakeFiles/poly_test.dir/build

test/CMakeFiles/poly_test.dir/requires: test/CMakeFiles/poly_test.dir/polynomial_test.cpp.o.requires
.PHONY : test/CMakeFiles/poly_test.dir/requires

test/CMakeFiles/poly_test.dir/clean:
	cd /Users/joshua/Projects/sisl/build/test && $(CMAKE_COMMAND) -P CMakeFiles/poly_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/poly_test.dir/clean

test/CMakeFiles/poly_test.dir/depend:
	cd /Users/joshua/Projects/sisl/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joshua/Projects/sisl /Users/joshua/Projects/sisl/test /Users/joshua/Projects/sisl/build /Users/joshua/Projects/sisl/build/test /Users/joshua/Projects/sisl/build/test/CMakeFiles/poly_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/poly_test.dir/depend

