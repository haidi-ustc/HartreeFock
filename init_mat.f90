module alloc_dealloc 
use m_precision,only : dp
use FCS
public alloc,dealloc
contains
subroutine alloc(Nbasis,isp,IOP)
implicit none

integer:: Nbasis,isp,I,numbers,IOP
allocate(S(Nbasis,Nbasis))
allocate(V(Nbasis,Nbasis))
allocate(T(Nbasis,Nbasis))
allocate(P(Nbasis,Nbasis,isp))
allocate(Coeff(Nbasis,Nbasis,isp))
allocate(Eig(Nbasis,isp))
allocate(H_core(Nbasis,Nbasis))
!allocate(H_mo(Nbasis,Nbasis))
allocate(F_atom(Nbasis,Nbasis,isp))
allocate(F_mol(Nbasis,Nbasis,isp))
allocate(ERI_ao(Nbasis,Nbasis,Nbasis,Nbasis))
allocate(ERI_mo(Nbasis,Nbasis,Nbasis,Nbasis))
allocate(ansym_ERI(2*Nbasis,2*Nbasis,2*Nbasis,2*Nbasis))
allocate(Spin_Fock(2*Nbasis,2*Nbasis))
allocate(H_mo(Nbasis,Nbasis))

Numbers=Nbasis*(Nbasis+1)/2*((Nbasis*(Nbasis+1)/2)+1)/2
allocate(twoe(Numbers))
allocate(twoe_mo(Numbers))
allocate(Ioff(Numbers+2))

Ioff(1)=0
do I=1,Numbers
 Ioff(I+1)=Ioff(I)+I
end do

if(IOP>=5)then
 print*,"Ioff"
 print*,Ioff
 print*,"size(twoe)=",size(twoe) ! 406 for sto-3g h2o
end if
S=0.0_dp
V=0.0_dp
T=0.0_dp
ERI_ao=0.0_dp
ERI_mo=0.0_dp
ansym_ERI=0.0_dp
H_core=0.0_dp
twoe=0.0_dp
Spin_Fock=0.0_dp
end subroutine

subroutine dealloc()

deallocate(Ioff)
deallocate(S)
deallocate(P)


end subroutine

end module
