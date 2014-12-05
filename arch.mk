# Modify these for your Fortran compiler:
.SUFFIXES:
.SUFFIXES:.f .F .o .a .f90 .F90

# GFortran
FC = ifort 
#---------origin--------------------
#FFLAGS = -Wall -std=f95 -Wextra -Wimplicit-interface -fPIC
#FFLAGS += -g -fbounds-check
#---------siesta xmqin config-------
#FFLAGS=-I/share/apps/intel/mkl/include -I/share/apps/intel/mkl/include/intel64/lp64/ -g -O2
#LDFLAGS= -static-intel -L/share/apps/intel/mkl/lib/intel64 
#----------siesta auto confit---------
FFLAGS=-g -O1
FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT
#----------------------------------------
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)
# Release flags:
#F90FLAGS += -O3 -march=native -ffast-math -funroll-loops
#---------------------------LIBS-----------------------------------
BLAS_LIBS=
LAPACK_LIBS=
linalg= linalg.a 
libdft= libdftatom.a
libfdf= libfdf.a
COMP_LIB=$(libfdf) #$(linalg)  $(libdft)

LIBS= -Wl ,--start-group $(libdft) $(BLAS_LIBS) $(LAPACK_LIBS) $(linalg) -Wl, --end group

FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT #-DMPI
FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

#Dependency rules are created by autoconf according to whether
#discrete preprocessing is necessary or not
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

