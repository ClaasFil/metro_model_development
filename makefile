# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2

# Source files
SRC_DIR = /home/fillies/Documents/Uni_Projects/metro_model_development/src/fortran
HELPER_DIR = $(SRC_DIR)/helper
MOD_DIR = /home/fillies/Documents/Uni_Projects/metro_model_development/bin

# Source files
SOURCES = $(HELPER_DIR)/utilities.f90 \
		  $(HELPER_DIR)/finitdifferences.f90 \
		  $(HELPER_DIR)/poisson_solver_utilities.f90\
		  $(SRC_DIR)/dry_convection.f90 

# Object files
OBJECTS = $(SOURCES:.f90=.o)

# Executable name
EXEC = dry_convection

# Default target
all: $(EXEC)

# Dependencies for modules
$(HELPER_DIR)/utilities.o: $(HELPER_DIR)/utilities.f90
$(HELPER_DIR)/finitdifferences.o: $(HELPER_DIR)/finitdifferences.f90 $(HELPER_DIR)/utilities.o
$(HELPER_DIR)/poisson_solver_utilities.o: $(HELPER_DIR)/poisson_solver_utilities.f90 $(HELPER_DIR)/utilities.o
$(SRC_DIR)/dry_convection.o: $(SRC_DIR)/dry_convection.f90 $(HELPER_DIR)/utilities.o $(HELPER_DIR)/finitdifferences.o $(HELPER_DIR)/poisson_solver_utilities.o

# Rule to create the executable
$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) -I$(MOD_DIR) -o $@ $(OBJECTS)

# Rule to compile Fortran files
%.o: %.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR)  -J$(MOD_DIR) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJECTS) $(EXEC)

# Phony targets
.PHONY: all clean