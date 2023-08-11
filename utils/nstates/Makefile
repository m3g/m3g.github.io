# configure generated Makefile
#
# Makefile for Packmol: Read the comments if you have some
#                       problem while compiling.
#
# You may use the ./configure script to search automatically for
# some fortran compiler.
#
# This make file will try to compile packmol with the default
# fortran compiler, defined by the FC directive. For doing this,
# just type
#
#          make 
#
# If you want to compile with some specific fortran compiler, you must 
# change the line below to the path of your fortran compiler. 
#
FORTRAN = /usr/bin/gfortran
#
# Change "AUTO" to the fortran command you want. 
#
# Change the flags of the compilation if you want:
#
FLAGS= -O3 -ffast-math 
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################
#
# Get the default fortran compiler
#
ifeq ($(FORTRAN),AUTO)
FORTRAN = $(FC)
endif 
#
# Flags for compiling development version
#
GENCANFLAGS := $(FLAGS)
ifeq ($(MAKECMDGOALS),devel)
FLAGS = -Wall -fcheck=bounds -g -fbacktrace -ffpe-trap=zero,overflow,underflow
GENCANFLAGS = -fcheck=bounds -g -fbacktrace -ffpe-trap=zero,overflow,underflow 
endif

all : numberofstates.f90
	$(FORTRAN) -o numberofstates numberofstates.f90 $(FLAGS)

devel : numberofstates.f90
	$(FORTRAN) -o numberofstates numberofstates.f90 $(FLAGS)
