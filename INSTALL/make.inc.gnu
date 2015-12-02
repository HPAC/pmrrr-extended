# Compiler for C and Fortran
CC = mpicc

# Compiler flags
CFLAGS = -pthread -O3

# Archiver and flags used when building the archive
AR = /usr/bin/ar 
ARFLAGS = rcs

# On some systems 'spinlocks' are not supported, therefore 
# here the flag to use 'mutexes' instead; default value is 1
SPINLOCK_SUPPORT = 1

