#!/bin/bash

# Set the directory for the compiled modules and program
BIN_DIR="bin"

# Ensure the bin directory exists
mkdir -p $BIN_DIR

# Compile the utilities module
echo "Compiling utilities module..."
gfortran -c -J$BIN_DIR src/fortran/helper/utilities.f90 -o $BIN_DIR/utilities.o

# Check if utilities compiled successfully
if [ $? -ne 0 ]; then
    echo "Failed to compile utilities module."
    exit 1
fi



# Compile the finite differences module
echo "Compiling finite differences module..."
gfortran -c -J$BIN_DIR src/fortran/helper/finitdifferences.f90 -o $BIN_DIR/finitdifferences.o

# Check if finite differences compiled successfully
if [ $? -ne 0 ]; then
    echo "Failed to compile finite differences module."
    exit 1
fi

# Compile the main program
echo "Compiling the main program..."
gfortran -c -I$BIN_DIR src/fortran/advection2D.f90 -o $BIN_DIR/advection2D.o

# Check if main program compiled successfully
if [ $? -ne 0 ]; then
    echo "Failed to compile the main program."
    exit 1
fi

# Link all object files to create the executable
echo "Linking object files to create executable..."
gfortran $BIN_DIR/utilities.o $BIN_DIR/finitdifferences.o $BIN_DIR/advection2D.o -o $BIN_DIR/advection2D

# Check if linking was successful
if [ $? -ne 0 ]; then
    echo "Failed to link object files."
    exit 1
fi

echo "Compilation and linking completed successfully :0"

# Execute the program
echo "Running the program..."
./$BIN_DIR/advection2D

# Check if the program ran successfully
if [ $? -ne 0 ]; then
    echo "Failed to run the program."
    exit 1
fi

echo "Program executed successfully :~)"
