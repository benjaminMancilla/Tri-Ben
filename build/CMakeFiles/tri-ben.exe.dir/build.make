# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/benjamin/Desktop/Tri-Ben

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/benjamin/Desktop/Tri-Ben/build

# Include any dependencies generated for this target.
include CMakeFiles/tri-ben.exe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tri-ben.exe.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tri-ben.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tri-ben.exe.dir/flags.make

CMakeFiles/tri-ben.exe.dir/src/main.cpp.o: CMakeFiles/tri-ben.exe.dir/flags.make
CMakeFiles/tri-ben.exe.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/tri-ben.exe.dir/src/main.cpp.o: CMakeFiles/tri-ben.exe.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/benjamin/Desktop/Tri-Ben/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tri-ben.exe.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tri-ben.exe.dir/src/main.cpp.o -MF CMakeFiles/tri-ben.exe.dir/src/main.cpp.o.d -o CMakeFiles/tri-ben.exe.dir/src/main.cpp.o -c /home/benjamin/Desktop/Tri-Ben/src/main.cpp

CMakeFiles/tri-ben.exe.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tri-ben.exe.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/benjamin/Desktop/Tri-Ben/src/main.cpp > CMakeFiles/tri-ben.exe.dir/src/main.cpp.i

CMakeFiles/tri-ben.exe.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tri-ben.exe.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/benjamin/Desktop/Tri-Ben/src/main.cpp -o CMakeFiles/tri-ben.exe.dir/src/main.cpp.s

# Object files for target tri-ben.exe
tri__ben_exe_OBJECTS = \
"CMakeFiles/tri-ben.exe.dir/src/main.cpp.o"

# External object files for target tri-ben.exe
tri__ben_exe_EXTERNAL_OBJECTS =

tri-ben.exe: CMakeFiles/tri-ben.exe.dir/src/main.cpp.o
tri-ben.exe: CMakeFiles/tri-ben.exe.dir/build.make
tri-ben.exe: CMakeFiles/tri-ben.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/benjamin/Desktop/Tri-Ben/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tri-ben.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tri-ben.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tri-ben.exe.dir/build: tri-ben.exe
.PHONY : CMakeFiles/tri-ben.exe.dir/build

CMakeFiles/tri-ben.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tri-ben.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tri-ben.exe.dir/clean

CMakeFiles/tri-ben.exe.dir/depend:
	cd /home/benjamin/Desktop/Tri-Ben/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/benjamin/Desktop/Tri-Ben /home/benjamin/Desktop/Tri-Ben /home/benjamin/Desktop/Tri-Ben/build /home/benjamin/Desktop/Tri-Ben/build /home/benjamin/Desktop/Tri-Ben/build/CMakeFiles/tri-ben.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tri-ben.exe.dir/depend

