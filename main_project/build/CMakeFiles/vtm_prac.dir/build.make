# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build

# Include any dependencies generated for this target.
include CMakeFiles/vtm_prac.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/vtm_prac.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/vtm_prac.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vtm_prac.dir/flags.make

CMakeFiles/vtm_prac.dir/src/main.cpp.o: CMakeFiles/vtm_prac.dir/flags.make
CMakeFiles/vtm_prac.dir/src/main.cpp.o: /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/main.cpp
CMakeFiles/vtm_prac.dir/src/main.cpp.o: CMakeFiles/vtm_prac.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vtm_prac.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/vtm_prac.dir/src/main.cpp.o -MF CMakeFiles/vtm_prac.dir/src/main.cpp.o.d -o CMakeFiles/vtm_prac.dir/src/main.cpp.o -c /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/main.cpp

CMakeFiles/vtm_prac.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/vtm_prac.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/main.cpp > CMakeFiles/vtm_prac.dir/src/main.cpp.i

CMakeFiles/vtm_prac.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/vtm_prac.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/main.cpp -o CMakeFiles/vtm_prac.dir/src/main.cpp.s

CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o: CMakeFiles/vtm_prac.dir/flags.make
CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o: /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/householder_for_sym_matrix.cpp
CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o: CMakeFiles/vtm_prac.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o -MF CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o.d -o CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o -c /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/householder_for_sym_matrix.cpp

CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/householder_for_sym_matrix.cpp > CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.i

CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/src/householder_for_sym_matrix.cpp -o CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.s

# Object files for target vtm_prac
vtm_prac_OBJECTS = \
"CMakeFiles/vtm_prac.dir/src/main.cpp.o" \
"CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o"

# External object files for target vtm_prac
vtm_prac_EXTERNAL_OBJECTS =

vtm_prac: CMakeFiles/vtm_prac.dir/src/main.cpp.o
vtm_prac: CMakeFiles/vtm_prac.dir/src/householder_for_sym_matrix.cpp.o
vtm_prac: CMakeFiles/vtm_prac.dir/build.make
vtm_prac: CMakeFiles/vtm_prac.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable vtm_prac"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vtm_prac.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vtm_prac.dir/build: vtm_prac
.PHONY : CMakeFiles/vtm_prac.dir/build

CMakeFiles/vtm_prac.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vtm_prac.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vtm_prac.dir/clean

CMakeFiles/vtm_prac.dir/depend:
	cd /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build /Users/mihailprotopopov/PROJECTS/vtm_prac/main_project/build/CMakeFiles/vtm_prac.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/vtm_prac.dir/depend

