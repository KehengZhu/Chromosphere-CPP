# Building Instructions for `chromo_standalone`

This document provides step-by-step instructions to compile and run the `chromo_standalone` project.

## Source Files

util/chromo_interface.cpp

util/chromo_util.cpp

chromosphere.cpp

chromosphere.hpp

ModChromosphereCPP.f90

chromo_main.f90

## Prerequisites

Ensure that the following software is installed on your system:

- **CMake** (version 3.12 or higher)
- **g++** (a C++ compiler that supports C++14)
- **gfortran or NAGFor** (a Fortran compiler)

## Compilation Instructions

Follow these steps to compile and run the `chromo_standalone` project:

1. **Create a Build Directory:**

   First, create a new directory for building the project. This keeps the build files separate from the source files.

   ```bash
   mkdir chromo_build
   ```

2. **Navigate to the Build Directory:**

   Change your current directory to the newly created `chromo_build` directory.

   ```bash
   cd chromo_build
   ```

3. **Run CMake:**

   In the `chromo_build` directory, run CMake to configure the project. This command will generate the necessary build files.

   ```bash
   cmake ..
   ```

   - CMake will automatically detect and use `gfortran` as the Fortran compiler.
   - It will also set up the project to use the C++14 standard for C++ files.

4. **Build the Project:**

   Use the `make` command to compile the project. This will generate the `chromo_main` executable.

   ```bash
   make
   ```

   - The `make` command will compile the source files and link them together. The `chromo_main` executable will be created in the `chromo_build` directory.

5. **Run the Executable:**

   Once the build process is complete, you can run the `chromo_main` executable.

   ```bash
   ./chromo_main
   ```

   This will execute the `chromo_standalone` program.

## Project Structure Overview

Here’s a brief overview of the project structure based on the `CMakeLists.txt`:

- **Include Directories:**
  - The project includes header files from the `include` directory.

- **Libraries:**
  - A library named `chromosphere` is created from the source files in the `util` directory and `chromosphere.cpp`.

- **Executable:**
  - The main program executable `chromo_main` is built from the Fortran source files `chromo_main.f90` and `ModChromophereCPP.f90`.
  - This executable is linked with the `chromosphere` library.

## Debug Mode

The project is configured to build in **Debug** mode by default, which enables debugging symbols and disables optimizations. If you wish to change this to a different build type (e.g., Release), modify the `CMakeLists.txt` by changing the line:

```cmake
set(CMAKE_BUILD_TYPE DEBUG)
```

to

```cmake
set(CMAKE_BUILD_TYPE RELEASE)
```

Then, re-run the CMake and build steps.

## CMakeLists.txt Breakdown

```cmake
cmake_minimum_required(VERSION 3.12)
project(chromo_standalone)
```

- **`cmake_minimum_required(VERSION 3.12)`:**
  - This command specifies the minimum version of CMake required to build the project. In this case, CMake version 3.12 or higher is required. This ensures that the project takes advantage of features and fixes available in CMake 3.12 and later versions.

- **`project(chromo_standalone)`:**
  - This command names the project `chromo_standalone`. It sets up the environment for the project, such as defining project-specific variables and properties.

```cmake
# debug mode
set(CMAKE_BUILD_TYPE DEBUG)
```

- **`set(CMAKE_BUILD_TYPE DEBUG)`:**
  - This line sets the build type to `DEBUG`. When building in debug mode, the compiler includes debugging symbols in the binaries and disables optimizations. This makes it easier to debug the program but may result in slower performance. You can change this to `RELEASE` for an optimized build.

```cmake
# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
```

- **`set(CMAKE_CXX_STANDARD 14)`:**
  - This command sets the C++ standard to C++14. It ensures that the C++ compiler will use the C++14 standard to compile the project's C++ files.

- **`set(CMAKE_CXX_STANDARD_REQUIRED True)`:**
  - This command enforces the use of the specified C++ standard (C++14). If the compiler does not support C++14, CMake will generate an error instead of falling back to an earlier standard.

```cmake
# set fortran compiler as gfortran
set(CMAKE_Fortran_COMPILER gfortran)
```

- **`set(CMAKE_Fortran_COMPILER gfortran)`:**
  - This line explicitly sets the Fortran compiler to `gfortran`. This is important because your project includes Fortran source files, so you need to specify the appropriate compiler.

```cmake
# export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
```

- **`set(CMAKE_EXPORT_COMPILE_COMMANDS ON)`:**
  - This command tells CMake to generate a `compile_commands.json` file. This file contains detailed information about the compilation commands used for each source file in the project. It is particularly useful for integration with development tools like code editors and static analysis tools.

```cmake
# set languages: cxx and fortran
enable_language(CXX)
enable_language(Fortran)
```

- **`enable_language(CXX)`:**
  - This command enables C++ as one of the languages for the project. It tells CMake to use a C++ compiler to compile C++ source files.

- **`enable_language(Fortran)`:**
  - Similarly, this command enables Fortran as a language for the project, instructing CMake to use a Fortran compiler for Fortran source files.

```cmake
# Add include dir
include_directories(include)
```

- **`include_directories(include)`:**
  - This command adds the `include` directory to the list of directories that the compiler will search for header files. This is necessary if your source files include header files located in this directory.

```cmake
# Add chromosphere library
add_library(chromosphere util/chromo_util.cpp util/chromo_interface.cpp chromosphere.cpp)
```

- **`add_library(chromosphere util/chromo_util.cpp util/chromo_interface.cpp chromosphere.cpp)`:**
  - This command creates a static library named `chromosphere` from the specified source files (`chromo_util.cpp`, `chromo_interface.cpp`, `chromosphere.cpp`). The library will be compiled and archived, and can later be linked to an executable or another library.

```cmake
# Add the executable for chromo_main.cpp
add_executable(chromo_main chromo_main.f90 ModChromophereCPP.f90)
```

- **`add_executable(chromo_main chromo_main.f90 ModChromophereCPP.f90)`:**
  - This command creates an executable named `chromo_main` from the specified Fortran source files (`chromo_main.f90`, `ModChromophereCPP.f90`). The executable is the main program that will be run after compilation.

```cmake
# Link the executable to the chromosphere library
target_link_libraries(chromo_main chromosphere)
```

- **`target_link_libraries(chromo_main chromosphere)`:**
  - This command links the `chromo_main` executable with the `chromosphere` library. This means that the `chromo_main` executable can use the functions and classes defined in the `chromosphere` library.

```cmake
# Set linker language for the executable
set_target_properties(chromo_main PROPERTIES LINKER_LANGUAGE Fortran)
```

- **`set_target_properties(chromo_main PROPERTIES LINKER_LANGUAGE Fortran)`:**
  - This command explicitly sets the linker language for the `chromo_main` executable to Fortran. Since the main entry point is a Fortran source file, this ensures that the correct linker settings are used.
