# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
#
#  Makefile for Libraries
#  Some systems only need the stuff in dc_lapack.f (Divide and Conquer)
#
default: linalg.a
#
include ../arch.mk
RANLIB=ranlib
#
ALL_OBJS= lapack.o dc_lapack.o blas.o 
DC_OBJS=  dc_lapack.o   
#
default: linalg.a
#
libblas.a: blas.o
	ar $(ARFLAGS_EXTRA) cru libblas.a blas.o
	$(RANLIB) libblas.a
	mv libblas.a ..
liblapack.a: lapack.o machine.o
	ar $(ARFLAGS_EXTRA) cru liblapack.a lapack.o machine.o
	$(RANLIB) liblapack.a
	mv liblapack.a ..
linalg.a: $(ALL_OBJS)
	ar $(ARFLAGS_EXTRA) cru linalg.a $(ALL_OBJS)
	$(RANLIB) linalg.a
	@(chmod +rx linalg.a)
	mv linalg.a ..
dc_lapack.a: $(DC_OBJS)
	ar $(ARFLAGS_EXTRA) cru dc_lapack.a $(DC_OBJS)
	$(RANLIB) dc_lapack.a
	mv dc_lapack.a ..
#
clean:
	rm -f $(ALL_OBJS) *.a




