#!/bin/bash

# Set the directory for the compiled modules and program
BIN_DIR="/home/fillies/Documents/Uni_Projects/metro_model_development/bin"

# Compiler and flags
FC="gfortran"
FFLAGS="-O2"
NETCDF_LIB="/usr/lib/x86_64-linux-gnu"  # Adjust as necessary
NETCDF_INC="/usr/include"               # Adjust as necessary

# Source and helper directories
SRC_DIR="/home/fillies/Documents/Uni_Projects/metro_model_development/src/fortran"
HELPER_DIR="$SRC_DIR/helper"

# Ensure the bin directory exists
mkdir -p $BIN_DIR

# Compile the utilities module
echo "Compiling utilities module..."
$FC $FFLAGS -I$NETCDF_INC -J$BIN_DIR -c $HELPER_DIR/utilities.f90 -o $BIN_DIR/utilities.o

# Compile the finite differences module
echo "Compiling finite differences module..."
$FC $FFLAGS -I$NETCDF_INC -J$BIN_DIR -c $HELPER_DIR/finitdifferences.f90 -o $BIN_DIR/finitdifferences.o

# Compile the Poisson solver utilities module
echo "Compiling Poisson solver utilities module..."
$FC $FFLAGS -I$NETCDF_INC -J$BIN_DIR -c $HELPER_DIR/poisson_solver_utilities.f90 -o $BIN_DIR/poisson_solver_utilities.o

# Compile the main program
echo "Compiling the main program..."
$FC $FFLAGS -I$NETCDF_INC -J$BIN_DIR -c $SRC_DIR/moist_convection.f90 -o $BIN_DIR/moist_convection.o

# Link all object files to create the executable
echo "Linking object files to create executable..."
$FC -o $BIN_DIR/moist_convection $BIN_DIR/utilities.o $BIN_DIR/finitdifferences.o $BIN_DIR/poisson_solver_utilities.o $BIN_DIR/moist_convection.o -L$NETCDF_LIB -lnetcdf -lnetcdff

# Check if the executable was created successfully
if [ -f "$BIN_DIR/moist_convection" ]; then
    echo "Compilation and linking completed successfully :)"
    # Execute the program
    echo "Running the program..."
    $BIN_DIR/moist_convection
else
    echo "Failed to create the executable."
fi
