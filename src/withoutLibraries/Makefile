 ## FC = the compiler to use
FC=mpiifort

PROGRAM = run.exe

#### Compiler options

#### Optimization
##FFLAGS = -O3 -ipo

# Suppressing warning about too long routine names
FFLAGS += -diag-disable 5462 -diag-disable 6536

# Turning on preprocessor
FFLAGS += -fpp

#FFLAGS += -g -traceback

#### Debugging options
# Checking everything
FFLAGS += -O0 -traceback -g -implicitnone -check uninit -check bounds -ftrapuv #-check arg_temp_created

# Robust debug result
FFLAGS += -init=snan -traceback -warn interfaces -fp-speculation safe -check uninit -check bounds -ftrapuv -g

## -- Disables optimizations. Default is -O2
#FFLAGS += -O0

## -- Produces symbolic debug information in object file
#FFLAGS += -g

## -- Enables debug information and control output of enhanced debug information. To use this option, you must also specify the -g option.
#FFLAGS += -debug all

## -- Tells the compiler to generate extra information in the object file to provide source file traceback information when a severe error occurs at runtime.
## Specifying -traceback will increase the size of the executable program, but has no impact on runtime execution speeds.
#FFLAGS += -traceback

## -- Checks for all runtime failures
#FFLAGS += -check all

## Trapping unitialized variables
#FFLAGS += -ftrapuv

## -- Generates code to perform runtime checks on array subscript and character substringexpressions
## Once the program is debugged, omit this option to reduce executable program size and slightly improve runtime performance.
#FFLAGS += -check bounds
## -- Checks for uninitialized scalar variables without the SAVE attribute.
#FFLAGS += -check uninit

## -- Puts temporary variables on the heap. This should be prevented due to performance
## Consider defining a local array with contiguous memory
#FFLAGS += -fp-stack-check -heap-arrays
## -- Checks if temporary arguments are created
#FFLAGS += -check arg_temp_created

#### Suffix-rules:
## Begin by throwing away all old suffix-
## rules, and then create new ones for compiling 
## *.f90-files.
.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FFLAGS) $<

#### Include the dependency-list created by makedepf90 below
include .depend

#### target 'clean' for deleting object- *.mod- and other unwanted files
clean:
	rm -f *.o *.mod *.exe *~

#### Create a dependency list using makedepf90.
depend .depend:
	makedepf90 -o $(PROGRAM) *.f90 > .depend

#install:
#	$(FC) $(FFLAGS) $(LIBS) "$(MKLROOT)/include/mkl_cdft.f90" -c -o mkl_cdft.o
