module  ccsd_solver
use m_precision,only : dp
implicit none

public  ccsd
private

contains 
subroutine ccsd(E_Corr_mp2,E_hf,Nocc,Nmo,Fock,ERI,iter_max,CCSD_eps,IOP,&
                ts_out,td_out,E_Corr_CCSD,isconverge)
implicit none
!-------------------------------------------------------------------------------
!integer  Nocc         Number of the occupied orbitals
!integer  Nmo          Number of the molecular orbitals
!integer  Nvir         Number of the virtual orbitals
!integer  IOP          I/O print level
!integer  N_DIIS       dimension of the DIIS subspace
!integer  iter_max     max iteration times for CCSD
!real(dp) Fock         MO basis Spin Fock matrix
!real(dp) ERI          MO basis ERI matrix
!
!
!-------------------------------------------------------------------------------
integer       :: a,b,i,j,n,IOP,iter,iter_max,N_DIIS
integer       :: Nocc,Nvir,Nmo  
logical       :: isconverge,DO_DIIS,flag
real(dp)      :: E_Corr_CCSD ,DE_Corr_CCSD,CCSD_EPS,E_Corr_mp2,E_hf,Old_E_Corr_CCSD
real(dp)      :: Fock(Nmo,Nmo),ERI(Nmo,Nmo,Nmo,Nmo)
real(dp),intent(out):: ts_out(Nmo-Nocc,Nocc),td_out(Nmo-Nocc,Nmo-Nocc,Nocc,Nocc)
real(dp),allocatable:: ts(:,:),ts_new(:,:),td(:,:,:,:),td_new(:,:,:,:),Fae(:,:),&
                       Fmi(:,:),Fme(:,:),Dai(:,:),Dabij(:,:,:,:),Wmnij(:,:,:,:),&
                       Wabef(:,:,:,:),Wmbej(:,:,:,:)

       write(6,'(/2a)') 'haidi :                 ',&
     &                    '=============================='
       write(6,'(30(" "),a)') 'Begin '//"CCSD"//' Module '
       write(6,'(2a)') '                        ',&
     &                    '=============================='


Nvir=Nmo-Nocc
allocate(ts(Nvir,Nocc))
allocate(ts_new(Nvir,Nocc))
allocate(Fae(Nvir,Nvir))
allocate(Fmi(Nocc,Nocc))
allocate(Fme(Nocc,Nvir))
allocate(Dai(Nvir,Nocc))
allocate(Dabij(Nvir,Nvir,Nocc,Nocc))
allocate(td(Nvir,Nvir,Nocc,Nocc))
allocate(td_new(Nvir,Nvir,Nocc,Nocc))
allocate(Wmnij(Nocc,Nocc,Nocc,Nocc))
allocate(Wmbej(Nocc,Nvir,Nvir,Nocc))
allocate(Wabef(Nvir,Nvir,Nvir,Nvir))

isconverge=.false.
if(IOP>=5)then
  print*,"Nmo=",Nmo
  print*,"Nocc=",Nocc
  print*,"Nvir=",Nvir
end if

! Make denominator arrays
Dai=0.0_dp
do a=1,nvir
   do i=1,nocc
      Dai(a,i)=Fock(i,i)-Fock(a+Nocc,a+Nocc)
   end do
end do

Dabij=0.0_dp
do i=1,Nocc
   do j=1,Nocc
      do a=1,Nvir
         do b=1,Nvir
            Dabij(a,b,i,j)=Fock(i,i)+Fock(j,j)-Fock(a+Nocc,a+Nocc)-Fock(b+Nocc,b+Nocc)
         end do
      end do
   end do
end do
!Init guess T1
Do i=1,Nocc
   Do a=1,Nvir
      ts(a,i)=Fock(a+Nocc,i)/Dai(a,i)
   end do
end do
!Init guess T2
do i=1,Nocc
   do j=1,Nocc
      do a=1,Nvir
         do b=1,Nvir
            td(a,b,i,j)=td(a,b,i,j)+ERI(i,j,a+Nocc,b+Nocc)/Dabij(a,b,i,j)
         end do
      end do
   end do
end do
!  E_ccsd=0.0_dp
   N_DIIS=4
!  Error1Set=[]  
!  Error2Set-[]
!  T1Set=[]
!  T2Set=[]
CALL PRINT_MATRIX( 'Init MP2 CorrE /Hartree :', 1, 1,E_Corr_mp2,1 )
CALL PRINT_MATRIX( 'Init MP2 Energy  /Hartree :', 1, 1,E_Corr_mp2+E_hf,1 )
print*,""
   E_Corr_ccsd=E_Corr_mp2
write(*,"(a)")"Iter      E_corr      Delata_E"

! Begin CCSD Loop
do iter=1,iter_max
   Old_E_Corr_CCSD=E_Corr_CCSD
   call updateintermediates(Nvir,Nocc,.true.,ts,td,Fock,ERI,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
   call makeT1(Nvir,Nocc,.true.,ts,td,Fock,ERI,Fae,Fmi,Fme,Dai,ts_new)
        ts=ts_new
   call makeT2(Nvir,Nocc,.true.,ts,td,Fock,ERI,Fae,Fmi,Fme,ts_new,Wmnij,Wabef,Wmbej,Dabij,td_new)
        td=td_new
   if(Do_DIIS)then
     !for diis
   end if
   E_Corr_CCSD=ccsd_E(Nocc,Nvir,ts,td,Fock,ERI)
   DE_Corr_CCSD=abs(E_Corr_CCSD-Old_E_Corr_CCSD)
   write(*,"(i3,f15.10,f15.10)")iter,E_Corr_CCSD,DE_Corr_CCSD

   if (DE_Corr_CCSD < CCSD_EPS)then
     print*,"Total Iteration:",iter
     isconverge=.true.

     CALL PRINT_MATRIX( 'CCSD CorrE /Hartree :', 1, 1,E_Corr_CCSD,1 )
     CALL PRINT_MATRIX( 'CCSD Energy /Hartree :', 1, 1,E_Corr_CCSD+E_hf,1 )

     ts_out=ts
     td_out=td
     deallocate(ts)
     deallocate(ts_new)
     deallocate(Fae)
     deallocate(Fmi)
     deallocate(Fme)
     deallocate(Dai)
     deallocate(Dabij)
     deallocate(td)
     deallocate(td_new)
     deallocate(Wmnij)
     deallocate(Wmbej)
     deallocate(Wabef)

     return
   end if

end do
! End Loop

write(*,"(a)") "CCSD failed to isconvergee to disired precision !"

deallocate(ts)
deallocate(ts_new)
deallocate(Fae)
deallocate(Fmi)
deallocate(Fme)
deallocate(Dai)
deallocate(Dabij)
deallocate(td)
deallocate(td_new)
deallocate(Wmnij)
deallocate(Wmbej)
deallocate(Wabef)

return
end subroutine CCSD

 
real(dp) function  taus(ts,td,Nocc,Nvir,a,b,i,j)
implicit none
integer :: a,b,i,j,Nvir,Nocc
real(dp) :: td(Nvir,Nvir,Nocc,Nocc),ts(Nvir,Nocc)
taus=td(a,b,i,j)+0.5_dp*(ts(a,i)*ts(b,j)-ts(b,i)*ts(a,j))
return
end function taus

real(dp) function  tau(ts,td,Nocc,Nvir,a,b,i,j)
implicit none
integer  :: a,b,i,j,Nvir,Nocc
real(dp) :: td(Nvir,Nvir,Nocc,Nocc),ts(Nvir,Nocc)
tau=td(a,b,i,j)+ts(a,i)*ts(b,j)-ts(b,i)*ts(a,j)
return 
end function tau

subroutine updateintermediates(Nvir,Nocc,flag,ts,td,Fock,ERI,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
implicit none
integer             :: i,j,a,b,e,f,m,n
integer,intent(in)  :: Nvir,Nocc
real(dp),intent(out):: Fae(Nvir,Nvir),Fmi(Nocc,Nocc),Fme(Nocc,Nvir),&
                       Wmnij(Nocc,Nocc,Nocc,Nocc),Wabef(Nvir,Nvir,Nvir,Nvir),&
                       Wmbej(Nocc,Nvir,Nvir,Nocc)
real(dp),intent(in) :: ts(Nvir,Nocc),td(Nvir,Nvir,Nocc,Nocc),Fock(Nvir+Nocc,Nvir+Nocc),&
                       ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc)
logical             :: flag
if(flag)then
  Fae=0.0_dp
  do a=1,Nvir
     do e=1,Nvir
        Fae(a,e)=(1.0_dp-delta(a,e))*Fock(a+Nocc,e+Nocc)
        do m=1,Nocc
           Fae(a,e)=Fae(a,e)-0.50_dp*Fock(m,e+Nocc)*ts(a,m)
           do f=1,Nvir
              Fae(a,e)=Fae(a,e)+ts(f,m)*ERI(m,a+Nocc,f+Nocc,e+Nocc)
              do n=1,Nocc
                 Fae(a,e)=Fae(a,e)-0.5_dp*taus(ts,td,Nocc,Nvir,a,f,m,n)*ERI(m,n,e+Nocc,f+Nocc)
              end do
           end do 
        end do
      end do
  end do
  
  Fmi=0.0_dp
  do m=1,Nocc
     do i=1,Nocc
        Fmi(m,i)=(1.0_dp-delta(m,i))*Fock(m,i)
        do e=1,Nvir
           Fmi(m,i)=Fmi(m,i)+0.50_dp*Fock(m,e+Nocc)*ts(e,i)
           do n=1,Nocc
              Fmi(m,i)=Fmi(m,i)+ts(e,n)*ERI(m,n,i,e+Nocc)
              do f=1,Nvir
                 Fmi(m,i)=Fmi(m,i)+0.5_dp*taus(ts,td,Nocc,Nvir,e,f,i,n)*ERI(m,n,e+Nocc,f+Nocc)
              end do
           end do
        end do
      end do
  end do
  
  Fme=0.0_dp
 do m=1,Nocc
    do e=1,Nvir
       Fme(m,e)=Fock(m,e+Nocc)
       do n=1,Nocc
          do f=1,Nvir
             Fme(m,e)=Fme(m,e)+ts(f,n)*ERI(m,n,e+Nocc,f+Nocc)
          end do
       end do
    end do
 end do

 Wmnij=0.0_dp
 do m=1,Nocc
    do n=1,Nocc
       do i=1,Nocc
          do j=1,Nocc
             Wmnij(m,n,i,j)=ERI(m,n,i,j)
             do e=1,Nvir
                Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+ts(e,j)*ERI(m,n,i,e+Nocc)-ts(e,i)*ERI(m,n,j,e+Nocc)
                do f=1,Nvir
                   Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+0.25_dp*tau(ts,td,Nocc,Nvir,e,f,i,j)&
                                  *ERI(m,n,e+Nocc,f+Nocc)
                end do
             end do
          end do
       end do
     end do
  end do
  
  Wabef=0.0_dp
  do a=1,Nvir
     do b=1,Nvir
        do e=1,Nvir
           do f=1,Nvir
              Wabef(a,b,e,f)=ERI(a+Nocc,b+Nocc,e+Nocc,f+Nocc)
              do m=1,Nocc
                 Wabef(a,b,e,f)=Wabef(a,b,e,f)-ts(b,m)*ERI(a+Nocc,m,e+Nocc,f+Nocc)&
                                 +ts(a,m)*ERI(b+Nocc,m,e+Nocc,f+Nocc)
                 do n=1,Nocc                            ! tau(ts,td,Nocc,Nvir,a,b,i,j)
                    Wabef(a,b,e,f)=Wabef(a,b,e,f)+0.25_dp*tau(ts,td,Nocc,Nvir,a,b,m,n)&
                                   *ERI(m,n,e+Nocc,f+Nocc)
                 end do
              end do
            end do
         end do
      end do
  end do
  
 Wmbej=0.0_dp
  do m=1,Nocc
     do b=1,Nvir
        do e= 1,Nvir
           do j=1,Nocc
              Wmbej(m,b,e,j)=ERI(m,b+Nocc,e+Nocc,j)
              do f=1,Nvir
                 Wmbej(m,b,e,j)= Wmbej(m,b,e,j)+ts(f,j)*ERI(m,b+Nocc,e+Nocc,f+Nocc)
              end do
              do n=1,Nocc
                 Wmbej(m,b,e,j)= Wmbej(m,b,e,j)-ts(b,n)*ERI(m,n,e+Nocc,j)
                 do f=1,Nvir
                    Wmbej(m,b,e,j)= Wmbej(m,b,e,j)-(0.5_dp*td(f,b,j,n)+ts(f,j)*ts(b,n))&
                                   *ERI(m,n,e+Nocc,f+Nocc)
                 end do
              end do
              
            end do
        end do
     end do
  end do
end if
return
end subroutine updateintermediates

real(dp) function delta(i,j) 
implicit none
integer  ::i,j
if(i==j)then
  delta=1.0_dp
else
  delta=0.0_dp
end if
return
end function delta

subroutine makeT1(Nvir,Nocc,flag,ts,td,Fock,ERI,Fae,Fmi,Fme,Dai,ts_new)
implicit none
integer             :: i,j,a,b,e,f,m,n
integer,intent(in)  :: Nvir,Nocc
real(dp),intent(in) :: ts(Nvir,Nocc),td(Nvir,Nvir,Nocc,Nocc),Fae(Nvir,Nvir),Fmi(Nocc,Nocc),&
                       Fme(Nocc,Nvir),Fock(Nvir+Nocc,Nvir+Nocc),Dai(Nvir,Nocc),&
                       ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc)
real(dp),intent(out):: ts_new(Nvir,Nocc)
logical             :: flag

if(flag)then
 ts_new=0.0_dp
 do a=1,Nvir
    do i=1,Nocc
       ts_new(a,i)=Fock(i,a+Nocc)
       do e=1,Nvir
          ts_new(a,i)=ts_new(a,i)+ts(e,i)*Fae(a,e)
       end do
       do m=1,Nocc
          ts_new(a,i)=ts_new(a,i)-ts(a,m)*Fmi(m,i)
          do e=1,Nvir
             ts_new(a,i)=ts_new(a,i)+td(a,e,i,m)*Fme(m,e)
             do f=1,Nvir
                ts_new(a,i)=ts_new(a,i)-0.5_dp*td(e,f,i,m)*ERI(m,a+Nocc,e+Nocc,f+Nocc)
             end do
             do n=1,Nocc
                ts_new(a,i)=ts_new(a,i)-0.5_dp*td(a,e,m,n)*ERI(n,m,e+Nocc,i)
             end do
          end do
       end do
       do n=1,Nocc
          do f=1,Nvir
             ts_new(a,i)=ts_new(a,i)-ts(f,n)*ERI(n,a+Nocc,i,f+Nocc)
          end do
       end do
       ts_new(a,i)=ts_new(a,i)/Dai(a,i)
    end do     
 end do     
end if
return
end subroutine makeT1


subroutine makeT2(Nvir,Nocc,flag,ts,td,Fock,ERI,Fae,Fmi,Fme,ts_new,Wmnij,Wabef,Wmbej,Dabij,td_new)
implicit none
integer              :: i,j,a,b,e,f,m,n
integer,intent(in)   ::Nvir,Nocc
real(dp),intent(in)  :: ts(Nvir,Nocc),td(Nvir,Nvir,Nocc,Nocc),Fae(Nvir,Nvir),Fmi(Nocc,Nocc),&
                        Fme(Nocc,Nvir),ts_new(Nvir,Nocc),Fock(Nvir+Nocc,Nvir+Nocc),&
                        ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc),&
                        Wmnij(Nocc,Nocc,Nocc,Nocc),Wabef(Nvir,Nvir,Nvir,Nvir),&
                        Wmbej(Nocc,Nvir,Nvir,Nocc),Dabij(Nvir,Nvir,Nocc,Nocc)
real(dp),intent(out) :: td_new(Nvir,Nvir,Nocc,Nocc)
logical  :: flag

if(flag)then
 td_new=0.0_dp
 do a=1,Nvir
    do b=1,Nvir
       do i=1,Nocc
          do j=1,Nocc
             td_new(a,b,i,j)=td_new(a,b,i,j)+ERI(i,j,a+Nocc,b+Nocc)
             do e=1,Nvir
                 td_new(a,b,i,j)=td_new(a,b,i,j)+td(a,e,i,j)*Fae(b,e)-td(b,e,i,j)*Fae(a,e)
                 do m=1,Nocc
                    td_new(a,b,i,j)=td_new(a,b,i,j)-0.5_dp*td(a,e,i,j)*ts(b,m)*Fme(m,e)+&
                                                    0.5_dp*td(a,e,i,j)*ts(a,m)*Fme(m,e)
                 end do
             end do
             do m=1,Nocc
                td_new(a,b,i,j)=td_new(a,b,i,j)-td(a,b,i,m)*Fmi(m,j)+td(a,b,j,m)*Fmi(m,i)
                do e=1,Nvir
                   td_new(a,b,i,j)=td_new(a,b,i,j)-0.5_dp*td(a,b,i,m)*ts(e,j)*Fme(m,e)+&
                                                   0.5_dp*td(a,b,i,m)*ts(e,i)*Fme(m,e)
                end do
             end do
             do e=1,Nvir
                td_new(a,b,i,j)=td_new(a,b,i,j)+ts(e,i)*ERI(a+Nocc,b+Nocc,e+Nocc,j)-ts(e,j)&
                                *ERI(a+Nocc,b+Nocc,e+Nocc,i)
                do f=1,Nvir
                   td_new(a,b,i,j)=td_new(a,b,i,j)+0.5_dp*tau(ts,td,Nocc,Nvir,e,f,i,j)*Wabef(a,b,e,f)
                end do
             end do
             do m=1,Nocc
                td_new(a,b,i,j)=td_new(a,b,i,j)-ts(a,m)*ERI(m,b+Nocc,i,j)+ts(b,m)*ERI(m,a+Nocc,i,j)
                do e=1,Nvir
                   td_new(a,b,i,j)=td_new(a,b,i,j)+td(a,e,i,m)*Wmbej(m,b,e,j)-ts(e,i)*ts(a,m)&
                                   *ERI(m,b+Nocc,e+Nocc,j)
                   td_new(a,b,i,j)=td_new(a,b,i,j)-td(a,e,j,m)*Wmbej(m,b,e,i)+ts(e,j)*ts(a,m)&
                                   *ERI(m,b+Nocc,e+Nocc,i)
                   td_new(a,b,i,j)=td_new(a,b,i,j)-td(b,e,i,m)*Wmbej(m,a,e,j)-ts(e,i)*ts(b,m)&
                                   *ERI(m,a+Nocc,e+Nocc,j)
                   td_new(a,b,i,j)=td_new(a,b,i,j)+td(b,e,j,m)*Wmbej(m,a,e,i)-ts(e,j)*ts(b,m)&
                                   *ERI(m,a+Nocc,e+Nocc,i)
                end do
                do n=1,Nocc
                   td_new(a,b,i,j)=td_new(a,b,i,j)+0.5_dp*tau(ts,td,Nocc,Nvir,a,b,m,n)*Wmnij(m,n,i,j)
                end do
             end do
             td_new(a,b,i,j)=td_new(a,b,i,j)/Dabij(a,b,i,j)
          end do
       end do
    end do
 end do
end if
return
end subroutine makeT2

real(dp) function ccsd_E(Nocc,Nvir,ts,td,Fock,ERI)
implicit none
integer               :: a,b,i,j
integer,intent(in)    :: Nvir,Nocc
real(dp),intent(in)   :: ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc),&
                         ts(Nvir,Nocc),td(Nvir,Nvir,Nocc,Nocc),Fock(Nvir+Nocc,Nvir+Nocc)
real(dp)              :: tmp
 tmp=0.0_dp
 do i=1,Nocc
     do a=1,Nvir
         tmp=tmp+Fock(i,a+Nocc)*ts(a,i)
            do j = 1,Nocc
               do b=1,Nvir
                  tmp=tmp+0.25_dp*ERI(i,j,a+Nocc,b+Nocc)*td(a,b,i,j)&
                         +0.5_dp*ERI(i,j,a+Nocc,b+Nocc)*ts(a,i)*ts(b,j)
               end do
            end do
     end do
 end do
 ccsd_E=tmp
 return
end function ccsd_E

end module
