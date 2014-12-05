module FCS
use m_precision, only: dp
implicit none

real(kind=dp),allocatable ::S(:,:),H_core(:,:),V(:,:),T(:,:),H_mo(:,:),&  
                            F_atom(:,:,:),Coeff(:,:,:),Eig(:,:), F_mol(:,:,:),& 
                            P(:,:,:),twoe(:),twoe_mo(:),Spin_Fock(:,:),&
                            ERI_ao(:,:,:,:),ERI_mo(:,:,:,:),ansym_ERI(:,:,:,:)
integer,allocatable       ::Ioff(:)
!                          L(:,:),L_inv(:,:),L_inv_T(:,:),S_inv(:,:),&
real(dp)                  :: E_nuc
end module FCS
