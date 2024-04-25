@echo off
rem Set the directory for the compiled modules and program
set BIN_DIR=bin

rem Ensure the bin directory exists
if not exist "%BIN_DIR%" mkdir "%BIN_DIR%"

rem Compile the utilities module
echo Compiling utilities module...
gfortran -c -J%BIN_DIR% src\fortran\helper\utilities.f90 -o %BIN_DIR%\utilities.o
if %ERRORLEVEL% neq 0 (
    echo Failed to compile utilities module.
    exit /b 1
)

rem Compile the finite differences module
echo Compiling finite differences module...
gfortran -c -J%BIN_DIR% src\fortran\helper\finitdifferences.f90 -o %BIN_DIR%\finitdifferences.o
if %ERRORLEVEL% neq 0 (
    echo Failed to compile finite differences module.
    exit /b 1
)

rem Compile the main program
echo Compiling the main program...
gfortran -c -I%BIN_DIR% src\fortran\poisson.f90 -o %BIN_DIR%\poisson.o
if %ERRORLEVEL% neq 0 (
    echo Failed to compile the main program.
    exit /b 1
)

rem Link all object files to create the executable
echo Linking object files to create executable...
gfortran %BIN_DIR%\utilities.o %BIN_DIR%\finitdifferences.o %BIN_DIR%\poisson.o -o %BIN_DIR%\poisson
if %ERRORLEVEL% neq 0 (
    echo Failed to link object files.
    exit /b 1
)

echo Compilation and linking completed successfully :0

rem Execute the program
echo Running the program...
%BIN_DIR%\poisson
if %ERRORLEVEL% neq 0 (
    echo Failed to run the program.
    exit /b 1
)

echo Program executed successfully :~)
