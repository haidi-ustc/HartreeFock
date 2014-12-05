module calc
use m_precision,only:dp
public calc_G,scal_mat,rms_rho,calc_E_elec, calc_Density, calc_Fock,trace,Eigfunc,&
       indexs,calc_S_half
private
contains

subroutine calc_S_half(S,Nb,X)
use constant,   only:Lwmax
implicit none
integer :: i,j,Nb
integer :: Lwork,Iwork(Lwmax),Liwork,info
real(dp):: S(Nb,Nb),LS(Nb,Nb),Lmda(Nb),work(Lwmax),&
           Lmda_hf(Nb),L_diag(Nb,Nb),L_diag_inv(Nb,Nb),S_hf_inv(Nb,Nb),&
           X(Nb,Nb)
!call random_seed()
!firstly: obtain the EigenValues and eigenVectors of S Matrix
LS=0.d0
do i=1,Nb
  do j=1,Nb
    if(i<=j)then
           LS(i,j)=S(i,j)
     end if
  end do
end do

!CALL PRINT_MATRIX( 'Upper Fock_prim Matrix :', Nb, Nb,C_prim,Nb )

Lwork=-1
iwork=-1
call DSYEVD('Vectors','Upper',Nb,LS,Nb,Lmda,work,Lwork,Iwork,Liwork,info)
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK(1) )
!solve eigenproblem
call DSYEVD('Vectors','Upper',Nb,LS,Nb,Lmda,work,Lwork,Iwork,Liwork,info)

if( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
END IF


if( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
END IF

CALL PRINT_MATRIX( 'LS Matrix :', Nb, Nb,LS,Nb )
CALL PRINT_MATRIX( 'Lmda Matrix :', 1, Nb,Lmda,1 )

L_diag=0.0d0
L_diag_inv=0.0d0
do i=1,Nb
   Lmda_hf(i)=dsqrt(Lmda(i))
   L_diag(i,i)=Lmda_hf(i)
   L_diag_inv(i,i)=1.0d0/Lmda_hf(i)
end do

CALL PRINT_MATRIX( 'Lmda_diag Matrix :', Nb, Nb,L_diag,Nb )
CALL PRINT_MATRIX( 'Lmda_diag_inv Matrix :', Nb, Nb,L_diag_inv,Nb )
X=matmul(LS,matmul((L_diag_inv),transpose(LS)))
!X=matmul(matmul(LS,L_diag_inv),transpose(LS))

CALL PRINT_MATRIX( 'X matrix :', Nb, Nb,X,Nb )
return
end subroutine 


real(dp) function trace(mat,dims)
implicit none
integer::i,dims
real(dp):: mat(dims,dims),sums
sums=0.d0
do i=1,dims
  sums=sums+mat(i,i)
end do
trace=sums
return
end function trace

subroutine  calc_G(P,IOP,twoe,Nb,G)
use FCS,        only:Ioff
implicit none 
integer::i,j,m,n,Nb,ij,mn,ijmn,im,jn,imjn,IOP
real(dp),intent(in) ::P(Nb,Nb),twoe(*)
real(dp),intent(out)::G(Nb,Nb)
real(dp)            ::sums
if(IOP>=5) print*,"  I   J   M   N   IJ  MN  IJMN  twoe"
do i=1,Nb
  do j=1,Nb
    G(i,j)=0.d0
    do m=1,Nb
      do n=1,Nb
         ij=indexs(i,j)
         mn=indexs(m,n)
         ijmn=indexs(ij,mn)
         im=indexs(i,m)
         jn=indexs(j,n)
         imjn=indexs(im,jn)
         G(i,j)=G(i,j)+(2.0_dp*twoe(ijmn)-twoe(imjn))*P(m,n)
         if(IOP>=5)then
            write(*,"(7i4,f12.8)")i,j,m,n,ij,mn,ijmn,twoe(ijmn)
         end if
!       we can also use the following method to read
!       ERI data
!       sums=sums+(2*ERI(i,j,m,n)-ERI(i,m,j,n))*P(m,n)
      end do
    end do
  end do
end do
return
end subroutine

subroutine scal_mat(alpha,Mat,Nb,aMat)
implicit none 
integer::i,j,Nb
real(dp)::alpha,Mat(Nb,Nb),aMat(Nb,Nb)
do i=1,Nb
  do j=1,Nb
     aMat(i,j)=alpha*Mat(i,j)
  end do
end do
return
end subroutine

real(dp) function rms_rho(D2,D1,Nb)
use  m_precision,only:dp
implicit none
integer::i,j,Nb
real(dp):: D1(Nb,Nb),D2(Nb,Nb),sums
real(dp),intrinsic:: dsqrt
sums=0.d0
do i=1,Nb
  do j=1,Nb
   sums=sums+(D2(i,j)-D1(i,j))**2.0_dp
  end do
end do
 rms_rho=dsqrt(sums)
return
end function


real(dp) function calc_E_elec(D,H,F,Nb)
implicit none
integer::i,j,Nb
real(dp):: sums,D(Nb,Nb),H(Nb,Nb),F(Nb,Nb)
sums=0.0_dp
do i=1,Nb
  do j=1,Nb
   sums=sums+D(i,j)*(H(i,j)+F(i,j))
  end do
end do
calc_E_elec=sums
return
end function


subroutine calc_Fock(H,Nb,D,F,ERI)
implicit none
integer::Nb,i,j,m,n
real(dp)::H(Nb,Nb),D(Nb,Nb),F(Nb,Nb),ERI(Nb,Nb,Nb,Nb),sums

do i=1,Nb
  do j=1,Nb
      sums=0.d0
      do m=1,Nb
         do n=1,Nb
           sums=sums+D(m,n)*(2*ERI(i,j,m,n)-ERI(i,m,j,n))
         end do
      end do
      F(i,j)=H(i,j)+sums
  end do
end do
return
end subroutine

subroutine calc_Density(C,P,Nb,occp,OLDP) !calc_Density(C,Nb,occp,D)
implicit none
integer  ::Nb,occp,i,j,m
real(dp) ::C(Nb,Nb),P(Nb,Nb),OLDP(Nb,Nb)

do i=1,Nb
  do j=1,Nb
      OLDP=P(i,j)
      P(i,j)=0.d0
      do m=1,occp
          P(i,j)= P(i,j)+C(i,m)*C(j,m)
      end do
  end do
end do
return
end subroutine

subroutine Eigfunc(F_prim,Nb,C_prim,E)
use constant, only  :  Lwmax
implicit none
integer :: i,j,m,n,Nb,IOP
real(dp):: C_prim(Nb,Nb),E(Nb),work(Lwmax),F_prim(Nb,Nb)
integer :: Lwork,Iwork(Lwmax),Liwork,info
C_prim=0.d0
do i=1,Nb
  do j=1,Nb
    if(i<=j)then
           C_prim(i,j)=F_prim(i,j)
     end if
  end do
end do

Lwork=-1
iwork=-1
call DSYEVD('Vectors','Upper',Nb,C_prim,Nb,E,work,Lwork,Iwork,Liwork,info)
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK(1) )
!solve eigenproblem
call DSYEVD('Vectors','Upper',Nb,C_prim,Nb,E,work,Lwork,Iwork,Liwork,info)

if( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
END IF
if(IOP==3)then
  CALL PRINT_MATRIX( 'Eigenvalues E :', 1, Nb,E,1 )
  CALL PRINT_MATRIX( 'Guess Coefficient C_prim  Matrix :', Nb, Nb,C_prim,Nb )
end if
!------end solve the eigen equation F'C'=C'E ---------------------------
return
end subroutine Eigfunc

integer function indexs(i,j)
use FCS,only : Ioff
implicit none
integer::i,j

if(i>j)then
  indexs=ioff(i)+j
else
  indexs=ioff(j)+i
end if
return
end function
end module calc
