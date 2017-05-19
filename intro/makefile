# Makefile for Fortran Tutorial

# Type "make clean" (and hit enter) to remove executable, object files
# Type "make" into the terminal to compile and create executable
# Type "./ex_prog" to run executable


.SUFFIXES: .f90 .o .mod
#File endings for Fortran 90 and mod files
OBJECTS = main.o \
          modf.o \
          cels2kelv.o \
          my_avg.o
#list of object files
FLAGS = -g -fdefault-real-8 -fdefault-integer-8 \
        -Wall -Wextra -fbounds-check
#FLAGS can be modified for each compiler, depending on use
#These are some debugging flags that will help catch errors
EXEC = ex_prog
#Name of executable
FC = gfortran
# GNU Fortran compiler

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) $(FLAGS) -o $(EXEC) $(OBJECTS) 

modf.o: modf.f90
	$(FC) $(FLAGS) -c $<;

%.o:	%.f90 modf.o
	$(FC) $(FLAGS) -c $<;

clean:
	rm -f *.o
	rm -f *.mod
	rm -f ex_prog
