# Makefile to build MPI Example 
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2016 
#

F90 = mpifort  
CFLAGS = 
COPTFLAGS = -O3 

# MPI Independent Targets

f90: f90.o
	$(F90) -o mpifortran mpifortran.o -fcray-pointer
  
f90.o:
	$(F90) -c mpifortran.f90 -fcray-pointer 

# Default Targets for Cleaning up the Environment
clean :
	rm *.o

pristine :
	rm *.o
	touch *.c  

ctags :
	ctags *.c

