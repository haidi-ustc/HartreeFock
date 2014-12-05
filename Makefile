.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

include arch.mk
#
default: what  $(COMP_LIB) HF final 

what:
	@echo
	@echo "Compilation is about to running !"
	@echo "If something is wrong ,pls modify the right Makefile"
	@echo "arch.mk file should be suitable for your OS"
	@echo
	@echo  "Hit ^C to abort..."
	@sleep 1
#--------------------------------------------------------
# all libdft obj files should be placed here
libdft_OBJS =  constants.o dft.o dft_data.o dftatom.o drivers.o \
        energies.o mesh.o mixings.o ode1d.o rdirac.o \
        reigen.o rpoisson.o rschroed.o \
        states.o types.o utils.o
$(libdft):$(libdft_OBJS)
	ar cru $(libdft) $(libdft_OBJS)
	@(chmod +rx libdftatom.a)
#---------------------------------------------------------
# new files should be placed here 
OBJS=  precision.o constant.o cc.o matrix.o init_mat.o debug.o readdata.o \
       io.o sys.o time.o print.o ao2mo.o  scf_diis.o  rhf.o   \
       m_fdf_global.o  calc.o  mp2.o ccpt.o HF.o  ccsd.o  \
       $(libfdf) 

 # cholesky.o io.o  $(libfdf) $(libdft) 

HF  : $(OBJS)
	$(FC) -o hf.x $(OBJS) -mkl 
#========================================================
#other libraries that might need to be compiled
$(linalg):
	@echo "==> Compiling $(linalg) in Math_lib..."
	@(cd Math_lib; make )
#---------------------------------------------------------
$(libfdf):
	@echo "==> Compiling $(libfdf) in fdf..."
	@(cd fdf ; make)
#---------------------------------------------------------
move:
	@(rm -rf objdir; mkdir objdir)
	@(mv *.mod *.a *.o objdir)
	@(cp hf.x objdir)
final:
	@echo ""
	@echo "make success!"
clean:
	@(rm -rf hf.x *.a *.o *.mod  atom *Standard )
#	@(rm -rf $(prog))
#	@(cd Math_lib ; make clean; cd ..)
	@(cd fdf ; make clean; cd ..)
	@(rm -rf ./bin ERI_mol*)
install:
	rm -rf ./bin
	mkdir ./bin
	cp -r data ./bin
	cp *.dat input hf.x  ./bin
pref=`date +"%B-%d-%Y-%X"`
tar:
	@(rm -rf  HartreeFock;mkdir HartreeFock;)
	@(cp -r README *.f90 *.F *.f input data fdf Math_lib arch.mk Makefile  HartreeFock)
	@(tar zcvf  HartreeFock.tar.gz HartreeFock/ ;rm -rf HartreeFock/)
	@(mv HartreeFock.tar.gz ../${pref}-HartreeFock.tar.gz)

#---------------------------------------------------------
#main program atom file denpendency rules
HF.o       : precision.o constant.o calc.o debug.o matrix.o ao2mo.o ccsd.o
#cholesky.o : precision.o
mp2.o      : precision.o calc.o
ccsd.o     : precision.o 
cc.o     : precision.o 
rhf.o      : precision.o constant.o calc.o matrix.o
#uhf.o      : precision.o constant.o calc.o
constant.o : precision.o
io.o       : sys.o
calc.o     : precision.o matrix.o constant.o
matrix.o   : precision.o
m_fdf_global.o: precision.o
#init_mat.o : precision.o matrix.o  
readdata.o :  precision.o matrix.o debug.o calc.o
ao2mo.o    :  precision.o matrix.o calc.o
timer.o    :  parallel.o sys.o
init_t2.o  : precision.o matrix.o
ccspt.o    : precision.o
