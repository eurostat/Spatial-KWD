#
# @fileoverview Copyright (c) 2019-2020, by Stefano Gualandi, UniPv,
#               via Ferrata, 1, Pavia, Italy, 27100
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)
#
#

COMPILER    = g++ ${OPTFLAG}
LINKER      = g++ ${LDFLAGS}

# Directory for my files
MYHOME          = ${PWD}
BIN             = ${MYHOME}/bin
INCLUDE         = ${MYHOME}/include
LIB             = ${MYHOME}/lib
SRC             = ${MYHOME}/src

OPTFLAG = -O3 -ffast-math -march=native -DNDEBUG -Wall -std=c++11 -fopenmp -DLINUX -Wall
LDFLAGS = -O3 -DNDEBUG -lm -pthread -std=c++11 -fopenmp

# Directory for output files
OUT_DIR=bin lib

MKLROOT = /opt/intel/mkl
MKL_INC = -DMKL_ILP64 -m64 -I${MKLROOT}/include #-DUSE_MKL
MKL_LIB = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl

laptop: ${OUT_DIR} ${SRC}/SolverCLI.cpp
	${COMPILER} -c -g -pg ${SRC}/SolverCLI.cpp -o ${LIB}/SolverCLI.o -I${INCLUDE} ${MKL_INC} -I./externs
	${LINKER} -o ${BIN}/solver ${LIB}/SolverCLI.o ${MKL_LIB}

# Create subdirectory for output files (bin,lib)
MKDIR_P = mkdir -p

lib:
	${MKDIR_P} lib

bin:
	${MKDIR_P} bin

# Be careful!
clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~