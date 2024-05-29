# .c makefile

# Home path
SU_HOME =/home/william/cwp
SU_HOME1 = /home/william/cwp/src/Complex
JOB_HOME = ..

# Compile path
GCC = /usr/bin/gcc
MPICC = /usr/bin/mpicc

# Include path
CFLAGS = -I$(SU_HOME)/include
CFLAGS1 = -I$(JOB_HOME)/include
CFLAGS2 = -I$(SU_HOME1)/include

# Library path
LDFLAGS = -L$(SU_HOME)/lib
LDFLAGS1 = -L$(JOB_HOME)/lib
LDFLAGS2 = -L$(SU_HOME1)/lib

LIB = -lsu -lpar -lcwp -lm
LIB1 = -lmyfun

CFILES0=main
CFILES1=DemRun
CFILES2=Element
CFILES3=ElementKernel
CFILES4=Contact
CFILES5=ContactKernel

EXECNAME=$(JOB_HOME)/bin/DEMEXE

all:
	-$(GCC) -fopenmp -o $(EXECNAME) $(CFILES0).c $(CFILES1).c $(CFILES2).c $(CFILES3).c $(CFILES4).c $(CFILES5).c \
	$(CFLAGS) $(CFLAGS1) \
	$(LDFLAGS) $(LDFLAGS1) \
	$(LIB) $(LIB1)
	rm -f *.o
clean::
	rm -f $(EXECNAME)

