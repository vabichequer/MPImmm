# Makefile to build MPI Example 
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2016 
#

CC = mpicc 
CFLAGS = 
COPTFLAGS = -O3 

# MPI Independent Targets

mpimmm : mpimmm.o  
	$(CC) mpimmm.o -o mpimmm -lm $(LIBFLAGS) 	

mpimmm.o : mpimmm.c
	$(CC) $(CFLAGS) $(COPTFLAGS) -c mpimmm.c

# Default Targets for Cleaning up the Environment
clean :
	rm *.o

pristine :
	rm *.o
	touch *.c  

ctags :
	ctags *.c

