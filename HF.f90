program Hartree_Fock_Roothan
!--------------------------------------------------------------------------------
! This program is used to solve hartree-Fock-Roothan Equation
! All of Matrix, including S (overlap matrix),V(electron-neucleus interaction),
! T(Kenetic Matrix),ERI( double electron integrals ),H_core(core hamiltonian) are
! obtained from files in the disk.
! 
!---------------------------------------------------------------------------------
!integer       :: i,j,m,n    general index number for atomic orbitals
!integer       :: ip,iq,ir,is    general index number for molecular orbitals
!integer       :: istat      used to check file read status
!integer       :: IOP        used to define the output verbosity 
!                             IOP=1  just output the simple information
!                             IOP=2  just output the detail information
!integer       :: Nelec      nuber of total electron
!integer       :: vir        number of virtual orbitals
!integer       :: occp       number of occupied orbitals
!integer       :: Nbasis         number of  the basis
!integer       :: max_iter   max iteration number 
!real(dp)      :: S(Nbasis,Nbasis)  overlap matrix
!real(dp)      :: V(Nbasis,Nbasis)  electron-neucleus interaction matrix
!real(dp)      :: T(Nbasis,Nbasis)  Kenetic matrix
!real(dp)      :: P(Nbasis,Nbasis)  convergeced density matrix
!real(dp)      :: H_core(Nbasis,Nbasis) core Hamitonian matrix
!real(dp)      :: ERI(Nbasis,Nbasis,Nbasis,Nbasis) electron repulsion integrals based on atomic orbitals
!real(dp)      :: ERI_mol(Nbasis,Nbasis,Nbasis,Nbasis) electron repulsion integrals based on molecular orbitals
!real(dp)      :: F_atom(Nbasis,Nbasis) Fock matrix based on atomic orbitals
!real(dp)      :: F_mol(Nbasis,Nbasis)  Fock matrix based on Molecular orbitals
!real(dp)      :: L(Nbasis,Nbasis)      Cholesky decomposition of S matrix S=trans(L)*L
!                               here L is Lower triangle matrix
!real(dp)      :: L_inv(Nbasis,Nbasis)  inverse of L matrix
!real(dp)      :: L_inv_T(Nbasis,Nbasis) inverse and transpose of L matrix
!real(dp)      :: S_inv(Nbasis,Nbasis) inverse of overlap matrix S
!real(dp)      :: Eig(Nbasis)  Eigenvalues of Fock matrix
!real(dp)      :: eps_E    Energy convergence precision
!real(dp)      :: eps_rho  Density matrix convergence precision
!real(dp)      :: E_hf     Hartree Fock Energy  
!real(dp)      :: E_mp2    MP2 total Energy    
!real(dp)      :: E_ccsd   CCSD total Energy 
!real(dp)      :: Ecorr_mp2 MP2 correlation Energy 
!real(dp)      :: Ecorr_mp2 CCSD correlation Energy 
!real(dp)      :: time_* timer  Variables
!real(dp)      :: 
!---------------------------------------------------------------------------------
use m_precision,only  : dp
use constant, only    : occp,Nelec,Lwmax,eps_e,eps_rho,max_iter,charlen,CC_max_iter,ccsd_eps
use calc
use sys 
use FCS
use fdf,only          : fdf_string,fdf_integer,fdf_double
use m_debug           
use m_ao2mo 
use ccsd_solver       
use coupled_cluster
use m_CCPT
use m_fdf_global,only : fdf_global_get
use alloc_dealloc,only: alloc,dealloc
implicit none
integer               ::i,j,m,n,ij,mn,ijmn,istat,fileid,ip,iq,ir,is,Nbasis,&
                        nmo,nvirt,vlen,nopen,nclosed
real(kind=dp)         ::val,E_hf,E_mp2,E_ccsd,Nelec_check,work(Lwmax),&
                        delta_Ne,nalpha,nbeta,Ecorr_ccsd,Ecorr_mp2
real(kind=dp)         ::temp,temp1,temp2,temp3,time_1,time_2
real(dp),allocatable  :: py_ERI(:,:,:,:)
integer               ::Lwork,Iwork(Lwmax),Liwork,info,atom_Z,isp,lun,&
                        a,b
character(len=charlen)::filein,fileout,sname_default,sname,&
                        slabel_default, slabel,scftype_default,&
                        scftype
character             :: fmts*15
logical               :: hf_convge,spin,ccsd_converg,NonIterTriple,test
integer,intrinsic     :: int,min
logical,external      :: LEQI
call timestamp("Beginning  >>>")

!----------------------for fdf test ---------------------------------
!for fdf init 
      filein ="input"
      fileout ='fdf.log'
      call fdf_init(filein,fileout)
!find the system name
      sname_default = 'Hartree Fock'
      sname = fdf_string('SystemName',sname_default)
!find the system label
      slabel_default  = 'HF'
      slabel = fdf_string('SystemLabel',slabel_default)
      call fdf_global_get(spin,'SpinPolarized',.false.)
      if(spin)then
        isp=2
      else
        isp=1
      end if
      
      Nbasis = fdf_integer('NumberOfBasis',0)
      if(Nbasis==0)then
          call  die("Number of basis error !")
      end if
      delta_Ne=fdf_double('TotalSpin',0.d0)
       Nalpha=(Nelec+delta_Ne)/2.d0
       Nbeta=(Nelec-delta_Ne)/2.d0 
!      call fdf_global_get(Atom_Z,'AtomicNumber', 0)
!      if(Atom_Z<=0 .or. Atom_Z>=119)then
!         call  die("Atomic Number error !")
!      end if

print*,"SystemName    ",sname
print*,"SystemLabel   ",slabel

call alloc(Nbasis,isp,IOP)

if(IOP>=5) print*, "ok alloc"
call readdata()

if(IOP>=5) print*, "ok readata"
!-------------------output matrix----------
CALL PRINT_MATRIX( 'E_nuc value :', 1,1,E_nuc,1 )
CALL PRINT_MATRIX( 'S Matrix :', Nbasis, Nbasis,S,Nbasis )
CALL PRINT_MATRIX( 'T Matrix :', Nbasis, Nbasis,T,Nbasis )
CALL PRINT_MATRIX( 'V Matrix :', Nbasis, Nbasis,V,Nbasis )
CALL PRINT_MATRIX( 'H_core Matrix :', Nbasis, Nbasis,H_Core,Nbasis )

deallocate(V)
deallocate(T)

scftype_default='RHF'
scftype=fdf_string('CalcType',scftype_default)

      if(leqi(scftype,'RHF')) scftype='rhf'
      if(leqi(scftype,'ROHF'))  scftype='rohf'
      if(leqi(scftype,'UHF')) scftype='uhf'

      if( (trim(scftype).ne.'rhf').and.(trim(scftype).ne.'rohf').and. &
         (trim(scftype).ne.'uhf') ) then

         write(6,'(/,2a,/,a)') &
         'size_name: Incorrect Calculation Type option specified,',&
         ' active options are:',&
         'RHF, ROHF, UHF '
         call die
      endif
!nuber of MO equal to Nbasis

nmo=Nbasis*2

  if(isp==1 .and. mod(int(Nelec),2)==0)then
     
     if(trim(scftype).eq.'rhf')then
        call rhf_scf(H_core,S,E_nuc,Nbasis,max_iter,eps_e,eps_rho,IOP,F_atom(:,:,isp),&
                 Coeff(:,:,isp),P(:,:,isp),Eig(:,isp),E_hf,hf_convge)

               if(hf_convge)then

                 print*,""
                 print*,"--------------Pre-process for MP2 and CCSD-----------"
                  CALL PRINT_MATRIX( 'Coeff Matrix :', Nbasis, Nbasis,Coeff(:,:,isp),Nbasis )
                  CALL PRINT_MATRIX( 'Fock  Matrix based on Atom:', Nbasis, Nbasis,F_atom(:,:,isp),Nbasis )
                  F_mol(:,:,isp)=matmul(matmul(transpose(Coeff(:,:,isp)),transpose(F_atom(:,:,isp))),Coeff(:,:,isp))
                  CALL PRINT_MATRIX( 'Fock  Matrix based on Atom:', Nbasis, Nbasis,F_mol(:,:,isp),Nbasis )

                  print*,"Number of electrons:"
                  Nelec_check=2*trace(matmul(P(:,:,isp),S),Nbasis)  !because of double occupied orbitals
                  print*,"Nelec=",Nelec
                  print*,""

!conver ao ERI to MO eRI 

                  call ao2mo(twoe,Nbasis,Coeff(:,:,isp),IOP,twoe_mo)
                  call mp2(twoe_mo,Nbasis,occp,Eig(:,isp),E_hf,Ecorr_mp2,IOP)
                  call spinmo(twoe_mo,Nbasis*2,IOP)
                  call spinFockEig(H_core,Coeff,Nbasis,Nbasis*2,int(Nelec),Eig(:,isp),IOP)

!                  call ccsd(E_hf,Ecorr_mp2,Nbasis*2,int(Nelec),CC_max_iter,IOP,E_ccsd)
                   allocate(ts(2*Nbasis-int(Nelec),Nbasis*2))
                   allocate(td(2*Nbasis-int(Nelec),2*Nbasis-int(Nelec),Nbasis*2,Nbasis*2))

                   call ccsd(Ecorr_mp2,E_hf,int(Nelec),Nbasis*2,Spin_Fock,ansym_ERI,&
                          CC_max_iter,CCSD_eps,IOP,ts,td,Ecorr_ccsd,ccsd_converg)
                            NonIterTriple=.true.
                            if(ccsd_converg .and. NonIterTriple) then
                            call ccpt(E_hf,Ecorr_ccsd,Spin_Fock,ansym_ERI,ts,td,int(Nelec),Nbasis*2)
                            end if


                   test=.false.
                   if(test)then 
                    allocate(py_ERI(Nbasis*2,Nbasis*2,Nbasis*2,Nbasis*2))
                    call io_assign(lun)
                    open(unit=lun,file='py_ERI.dat',status='old')
                    py_ERI=0.0_dp
                    do while(.true.)
                      if(istat .ne. 0)  exit
                      read(lun,*,iostat=istat)a,b,i,j,py_ERI(a,b,i,j)
                      write(*,"(4i4,f12.8)")i,j,a,b,py_ERI(i,j,a,b)
                    end do
                    call ccsd(Ecorr_mp2,E_hf,int(Nelec),Nbasis*2,Spin_Fock,py_ERI,&
                            CC_max_iter,ccsd_eps,IOP,ts,td,Ecorr_ccsd,ccsd_converg)
                   end if


               else
                   continue
               end if
      ! closed shell RHF 
     elseif(trim(scftype).eq.'uhf')then
!       call uhf_scf(H_core,S,ERI_atom,E_nuc,Nbasis,max_iter,eps_e,eps_rho,IOP,F_atom(:,:,isp),&
!                Coeff(:,:,isp),P(:,:,isp),Eig(:,isp),E_hf,hf_convge)      
     else
       call die("Unable run ROHF for Unpolarized Calc!")
     end if 
  end if
  if(isp==2 .and. mod(int(Nelec),2)/=0)then
     if(trim(scftype).eq.'rhf')then
       write(6,"(/,a,/)")'Unable to NonPolarized RHF !'
       call die
      ! closed shell RHF from system
     elseif((trim(scftype).eq.'uhf'))then
      ! closed shell URHF 
     else
      ! ROHF
     end if
  end if




 call dealloc()
call timestamp("End   >>")

end program 


