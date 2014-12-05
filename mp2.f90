subroutine mp2(twoe_mo,Nb,occp,Eig,E_hf,Ecorr_mp2,IOP)
use m_precision ,only:dp
!use FCS,         only:ERI_MO
use calc,        only:indexs
implicit none
integer  ::i,j,a,b,Nb,vir,occp,ia,ja,ib,jb,iajb,ibja,flag
real(dp) ::twoe_mo(*),Eig(Nb)
real(dp) ::tmp,fm,fz,E_mp2,Ecorr_mp2,E_hf
integer  ::IOP


       write(6,'(/2a)') 'haidi :                 ',&
     &                    '=============================='
       write(6,'(30(" "),a)') 'Begin '//"MP2"//' Module '
       write(6,'(2a)') '                        ',&
     &                    '=============================='


do i=1,occp
     do j=1,occp
        do a=occp+1,Nb
             ia=indexs(i,a)
             ja=indexs(j,a)
            do b=occp+1,Nb
                 ib=indexs(i,b)
                 jb=indexs(j,b)
                 iajb=indexs(ia,jb)
                 ibja=indexs(ib,ja)
                fm=Eig(i)+Eig(j)-Eig(a)-Eig(b)
!               fz=ERI_MO(i,a,j,b)*(2*ERI_MO(i,a,j,b)-ERI_MO(i,b,j,a))
                fz=twoe_mo(iajb)*(2*twoe_mo(iajb)-twoe_mo(ibja))
                tmp=tmp+fz/fm
            end do
        end do
    end do
end do
Ecorr_mp2= tmp

E_mp2=E_hf+Ecorr_mp2
  CALL PRINT_MATRIX( 'HF /Hartree :', 1, 1,E_hf,1 )
  CALL PRINT_MATRIX( 'MP2 correlation Energy /Hartree :', 1, 1,Ecorr_mp2,1 )
  CALL PRINT_MATRIX( 'E-MP2 /Hartree :', 1, 1,E_mp2,1 )

end subroutine

