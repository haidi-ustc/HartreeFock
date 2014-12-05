module constant
use m_precision
implicit none
integer,parameter:: occp=5,Lwmax=5000
real(dp)         :: eps=1.0d-8
integer,parameter::Nbasis=7
integer,parameter::max_iter=40
integer,parameter::CC_max_iter=60
real(dp)         ::eps_rho=1.0d-8,eps_e=1.0d-8
real(dp)         ::ccsd_eps=1.0D-8
integer,parameter::charlen=30
real(dp),parameter::Nelec=10.d0
end 
