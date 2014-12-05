subroutine rhf_scf(H,S,E_nuc,Nb,max_iter,eps_e,eps_rho,IOP,F,C,P,E,E_total,isconve)
use m_precision,only:  dp
use constant, only  :  occp , Lwmax
use calc
use FCS,only        :  twoe
implicit none
integer :: i,j,m,n,Nb,IOP
real(dp):: H(Nb,Nb),C_prim(Nb,Nb),E(Nb),work(Lwmax),&
           F(Nb,Nb),F_prim(Nb,Nb),P(Nb,Nb),OLDP_Cp(Nb,Nb),OLDP(Nb,Nb),&
           G(Nb,Nb),C(Nb,Nb),S(Nb,Nb),X(Nb,Nb)
real(dp)::delta_rho,delta_E,E_elec,E_nuc,E_total,eps_e,eps_rho,&
          NewE,E_elec_old,OldE
integer :: max_iter,iter,Lwork,Iwork(Lwmax),Liwork,info,N_diis,lun
logical :: isconve,do_diis
isconve=.false.
do_diis=.false.
!--------------debug-----------------------------

call calc_S_half(S,Nb,X)

if(do_diis)then
 write(*,"(a)")"DIIS is ON !"
else
 write(*,"(a)")"DIIS is OFF!"
end if 
!init P and G matrix 
P=0.d0
OldE=0.d0
G=0.d0

if(IOP>=5)then
  CALL PRINT_MATRIX( 'Init X Matrix :', Nb, Nb,X,Nb )
end if
       write(6,'(/2a)') 'haidi :                 ',&
     &                    '=============================='
       write(6,'(30(" "),a)') 'Begin '//"SCF"//' Module '                    
       write(6,'(2a)') '                        ',&
     &                    '==============================' 
write(*,"(5(a12,5x))")"Iter","E(elec)","E(tot)","Delta(E)","RMS(D)"

!Begin SCF LOOP
do iter=1,max_iter
   call calc_G(P,IOP,twoe,Nb,G)
   F=G+H
   if(do_diis)then
     ! for DIIS module
       call scf_diis()
   end if

   if(IOP >= 5)then
      CALL PRINT_MATRIX( 'G Matrix :', Nb, Nb,G,Nb )
      CALL PRINT_MATRIX( 'F Matrix :', Nb, Nb,F,Nb )
   end if
!------solve the eigen equation F'C'=C'E -----------------------------
F_prim=matmul(transpose(X),matmul(F,X))
call Eigfunc(F_prim,Nb,C_prim,E)
C=matmul(X,C_prim)
!---------------------------------------------------------------------

OLDP_Cp=OLDP
call calc_Density(C,P,Nb,occp,OLDP)

NewE=calc_E_elec(P,H,F,Nb)
E_total   =NewE+E_nuc
delta_rho=rms_rho(OLDP_Cp,OLDP,Nb)
delta_e=abs(NewE-OldE)
OldE=NewE

write(*,fmt=1000)iter,NewE,E_total,Delta_E,Delta_rho
1000 format(5x,i3,7x,f15.10,3x,f15.10,3x,f14.10,3x,f14.10,3x)
!--------------OUTPUT F' C' E C P------------------
if(IOP>5)then
  CALL PRINT_MATRIX( 'P  Matrix :', Nb, Nb,P,Nb )
  CALL PRINT_MATRIX( 'F  Matrix :', Nb, Nb,F,Nb )
  CALL PRINT_MATRIX( 'Elec  Energy/Hartree:', 1, 1,NewE,1 )
end if

if(delta_E < eps_e .and. delta_rho< eps_rho)then
   write(*,"(/,a)")"Density Matrix and Energy convergenced !"
   isconve=.true.
   CALL PRINT_MATRIX( 'Elec  Energy/Hartree:', 1, 1,NewE,1 )
   CALL PRINT_MATRIX( 'Total Energy/Hartree:', 1, 1,E_total,1 )
   CALL PRINT_MATRIX( 'G Matrix :', Nb, Nb,G,Nb )
   CALL PRINT_MATRIX( 'F Matrix :', Nb, Nb,F,Nb )
   CALL PRINT_MATRIX( 'C Matrix :', Nb, Nb,C,Nb )
   CALL PRINT_MATRIX( 'Eig Matrix :',1, Nb,E,1 )

     if(.true.)then
       call io_assign(lun)
       open(lun,file='myC.dat')
       do i=1,Nb
          do j=1,Nb
             write(lun,"(2i3,f15.8)")i,j,C(i,j)
          end do
       end do
       call io_close(lun)
       call io_assign(lun)
       open(lun,file='myP.dat')
       do i=1,Nb
          do j=1,Nb
             write(lun,"(2i3,f15.8)")i,j,P(i,j)
          end do
       end do
       call io_close(lun)
       call io_assign(lun)
       open(lun,file='myF.dat')
       do i=1,Nb
          do j=1,Nb
             write(lun,"(2i3,f15.8)")i,j,F(i,j)
          end do
       end do
       call io_close(lun)
     end if
   return
end if

 if(iter == max_iter-1) then
     write(*,"(/,a)")"Density Matrix and Energy disconvergenced !"
     return
 end if

end do

end subroutine rhf_scf 

