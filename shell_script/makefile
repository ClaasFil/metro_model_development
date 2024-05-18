# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2

# Source files
SRC_DIR = ../src/fortran
HELPER_DIR = $(SRC_DIR)/helper

# Source files
SOURCES = $(SRC_DIR)/dry_convection.f90 \
          $(HELPER_DIR)/utilities.f90 \
          $(HELPER_DIR)/finitdifferences.f90 \
          $(HELPER_DIR)/poisson_solver_utilities.f90

# Object files
OBJECTS = $(SOURCES:.f90=.o)

# Executable name
EXEC = dry_convection

# Default target
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

# Rule to compile Fortran files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJECTS) $(EXEC)

# Phony targets
.PHONY: all clean