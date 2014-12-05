module m_CCPT
use m_precision,only : dp
public ccpt
private

contains
subroutine ccpt(E_hf,Ecorr_ccsd,Fock,ERI,ts,td,Nocc,Nmo)
implicit none
integer,intent(in) :: Nmo,Nocc
real(dp),intent(in) :: E_hf,Ecorr_ccsd,ts(Nmo-Nocc,Nocc),td(Nmo-Nocc,Nmo-Nocc,Nocc,Nocc),&
                       Fock(Nmo,Nmo),ERI(Nmo,Nmo,Nmo,Nmo)
integer             :: i,j,k,a,b,c,Nvir
real(dp),allocatable:: Dabcijk(:,:,:,:,:,:),Tabcijk_C(:,:,:,:,:,:),Tabcijk_D(:,:,:,:,:,:)
real(dp)            :: E_CCPT
       write(6,'(/2a)') 'haidi :                 ',&
     &                    '=============================='
       write(6,'(30(" "),a)') 'Begin '//"CCSD(T)"//' Module '
       write(6,'(2a)') '                        ',&
     &                    '=============================='

Nvir=Nmo-Nocc
allocate(Dabcijk(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc))
allocate(Tabcijk_C(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc))
allocate(Tabcijk_D(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc))

Dabcijk=0.0_dp
do i=1,Nocc
   do j=1,Nocc
      do k=1,Nocc
         do a=1,Nvir
            do b=1,Nvir
               do c=1,Nvir
                  Dabcijk(a,b,c,i,j,k)=Fock(i,i)+Fock(j,j)+Fock(k,k)&
                           -Fock(a+Nocc,a+Nocc)-Fock(b+Nocc,b+Nocc)-Fock(c+Nocc,c+Nocc)
               end do
            end do
        end do
      end do
   end do
end do
call Calc_Tabcijk_C(ERI,Dabcijk,td,Nvir,Nocc,Tabcijk_C)
call Calc_Tabcijk_D(ERI,Dabcijk,ts,Nvir,Nocc,Tabcijk_D)
E_CCPT=Calc_CCPT_E(Nocc,Nvir,Dabcijk,Tabcijk_C,Tabcijk_D)

CALL PRINT_MATRIX( 'CCSD(T) CorrE /Hartree :', 1, 1,E_CCPT,1 )
CALL PRINT_MATRIX( 'CCSD(T) Energy /Hartree :',1, 1,E_CCPT+E_hf+Ecorr_ccsd,1 )

deallocate(Dabcijk,Tabcijk_C,Tabcijk_D)
end subroutine


subroutine Calc_Tabcijk_C(ERI,Dabcijk,td,Nvir,Nocc,Tabcijk_C)
implicit none
integer            :: Nocc,Nvir
real(dp),intent(in):: Dabcijk(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(out):: Tabcijk_C(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(in):: ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc)
real(dp),intent(in):: td(Nvir,Nvir,Nocc,Nocc)
integer            :: i,j,k,a,b,c,e,m
Tabcijk_C=0.d0
do i=1,Nocc
    do j=1,Nocc
       do k=1,Nocc
          do a=1,Nvir
             do b=1,Nvir
                do c=1,Nvir
                   do e=1,Nvir
                     Tabcijk_C(a,b,c,i,j,k)= Tabcijk_C(a,b,c,i,j,k)&
                     +td(a,e,j,k)*ERI(e+Nocc,i,b+Nocc,c+Nocc)-td(a,e,i,k)*ERI(e+Nocc,j,b+Nocc,c+Nocc)-td(a,e,j,i)*ERI(e+Nocc,k,b+Nocc,c+Nocc)&
                     -td(b,e,j,k)*ERI(e+Nocc,i,a+Nocc,c+Nocc)+td(b,e,i,k)*ERI(e+Nocc,j,a+Nocc,c+Nocc)+td(b,e,j,i)*ERI(e+Nocc,k,a+Nocc,c+Nocc)&
                     -td(c,e,j,k)*ERI(e+Nocc,i,b+Nocc,a+Nocc)+td(c,e,i,k)*ERI(e+Nocc,j,b+Nocc,a+Nocc)+td(c,e,j,i)*ERI(e+Nocc,k,b+Nocc,a+Nocc)
                      
                   end do         
                   do m=1,Nocc
                     Tabcijk_C(a,b,c,i,j,k)= Tabcijk_C(a,b,c,i,j,k)&
                     -td(b,c,i,m)*ERI(m,a+Nocc,j,k)+td(b,c,j,m)*ERI(m,a+Nocc,i,k)+td(b,c,k,m)*ERI(m,a+Nocc,j,i)&
                     +td(a,c,i,m)*ERI(m,b+Nocc,j,k)-td(a,c,j,m)*ERI(m,b+Nocc,i,k)-td(a,c,k,m)*ERI(m,b+Nocc,j,i)&
                     +td(b,a,i,m)*ERI(m,c+Nocc,j,k)-td(b,a,j,m)*ERI(m,c+Nocc,i,k)-td(b,a,k,m)*ERI(m,c+Nocc,j,i)
                   end do
                     Tabcijk_C(a,b,c,i,j,k)=Tabcijk_C(a,b,c,i,j,k)/Dabcijk(a,b,c,i,j,k)
                end do
             end do
          end do
       end do
    end do
end do

return
end subroutine


subroutine Calc_Tabcijk_D(ERI,Dabcijk,ts,Nvir,Nocc,Tabcijk_D)
implicit none
integer            :: Nocc,Nvir
real(dp),intent(in):: Dabcijk(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(out):: Tabcijk_D(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(in):: ERI(Nvir+Nocc,Nvir+Nocc,Nvir+Nocc,Nvir+Nocc)
real(dp),intent(in):: ts(Nvir,Nocc)
integer            :: i,j,k,a,b,c
Tabcijk_D=0.0_dp
do i=1,Nocc
    do j=1,Nocc
       do k=1,Nocc
          do a=1,Nvir
             do b=1,Nvir
                do c=1,Nvir
                   Tabcijk_D(a,b,c,i,j,k)= ( &
                   ts(a,i)*ERI(j,k,b+Nocc,c+Nocc)-ts(a,j)*ERI(i,k,b+Nocc,c+Nocc)-ts(a,k)*ERI(j,i,b+Nocc,c+Nocc)&
                  -ts(b,i)*ERI(j,k,a+Nocc,c+Nocc)+ts(b,j)*ERI(i,k,a+Nocc,c+Nocc)+ts(b,k)*ERI(j,i,a+Nocc,c+Nocc)&
                  -ts(c,i)*ERI(j,k,b+Nocc,a+Nocc)+ts(c,j)*ERI(i,k,b+Nocc,a+Nocc)+ts(c,k)*ERI(j,i,b+Nocc,a+Nocc)&
                   )/Dabcijk(a,b,c,i,j,k)
                 end do
              end do
           end do
       end do
    end do
end do

return
end subroutine Calc_Tabcijk_D


real(dp) function Calc_CCPT_E(Nocc,Nvir,Dabcijk,Tabcijk_C,Tabcijk_D)
implicit none
integer            :: Nocc,Nvir
real(dp),intent(in):: Dabcijk(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(in):: Tabcijk_D(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
real(dp),intent(in):: Tabcijk_C(Nvir,Nvir,Nvir,Nocc,Nocc,Nocc)
integer            :: i,j,k,a,b,c
Calc_CCPT_E=0.0_dp
do i=1,Nocc
    do j=1,Nocc
       do k=1,Nocc
          do a=1,Nvir
             do b=1,Nvir
                do c=1,Nvir
                   Calc_CCPT_E=Calc_CCPT_E+1.0_dp/36.0_dp*Tabcijk_C(a,b,c,i,j,k)*Dabcijk(a,b,c,i,j,k)*&
                               (Tabcijk_C(a,b,c,i,j,k)+Tabcijk_D(a,b,c,i,j,k))
                end do
             end do
          end do
       end do
    end do
end do
return
end function Calc_CCPT_E


end module m_CCPT
