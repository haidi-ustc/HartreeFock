module M_ao2mo
use m_precision,only  : dp
use calc       ,only  : indexs
public ao2mo,SpinFockEig,SpinMo
private
integer               :: i,j,k,l,ij,kl,klij,ijkl,p,q,r,s,pq,rs,pqrs,&
                          pr,qs,prqs,ps,qr,psqr
contains
subroutine ao2mo(twoe,Nb,C,IOP,twoe_mo)
use FCS        ,only  : ERI_ao,ERI_MO
implicit none
real(dp)              :: temp(Nb*(Nb+1)/2,Nb*(Nb+1)/2),X(Nb,Nb),Y(Nb,Nb),&
                         C(Nb,Nb),twoe_mo(*),twoe(*)
real(dp)              :: time_1,time_2
integer               :: flag,Nb,lun,IOP
!-------------2E transform ERI_AO to ERI_MOl------------------
!flag 1 2 3 are all right
! Algorithm    Time
!     1        1.27E-4
!     2        1.29E-2
!     3        6.67E-2
!------------------------------------------------------------- 
flag=1
if(flag==1)then
call cpu_time(time_1)
i=1
ij=1
do while(i<=Nb)
   j=1
   do while(j<=i)
      k=1
      kl=1
      do while(k<=Nb)
         l=1
         do while(l<=k)
            ijkl=indexs(ij,kl)
            X(k,l)=twoe(ijkl)
            X(l,k)=twoe(ijkl)
            kl=kl+1
            l=l+1
         end do
         k=k+1
      end do
         Y=0.d0
         Y=matmul(transpose(C),X)
         X=0.d0
         X=matmul(Y,C)
         kl=1
         do k=1,Nb
            do l=1,k
              temp(kl,ij)=X(k,l)
              kl=kl+1
            end do
         end do
      j=j+1
      ij=ij+1
   end do
   i=i+1
end do

k=1
kl=1
do while(k<=Nb)
   l=1
   do while(l<=k)
      X=0.d0
      y=0.d0
      i=1
      ij=1
      do while(i<=Nb)
         j=1
         do while(j<=i)
            X(i,j)=temp(kl,ij)
            X(j,i)=temp(kl,ij)
            ij=ij+1
            j=j+1
         end do
         i=i+1
      end do
         Y=0.d0
         Y=matmul(transpose(C),X)
         X=0.d0
         X=matmul(Y,C)
         ij=1
         do i=1,Nb
            do j=1,i
               klij=indexs(kl,ij)
               twoe_mo(klij)=X(i,j)
               ij=ij+1
            end do
         end do
      l=l+1
      kl=kl+1
   end do
   k=k+1
end do
call cpu_time(time_2)
print*,"flag1=",time_2-time_1
return
end if

if(flag==2)then
call cpu_time(time_1)
 ijkl=1
 do i=1,Nb
     do j=1,i
       do k=1,i
          if(i==k)then
            do l=1,j
              do p=1,Nb
                 do q=1,Nb
                   pq=indexs(p,q)
                      do r=1,Nb
                        do s=1,Nb
                           rs=indexs(r,s)
                           pqrs=indexs(pq,rs)
                           twoe_mo(ijkl)=twoe_mo(ijkl)+C(p,i)*C(q,j)*C(r,k)*C(s,l)*twoe(pqrs)
                        end do
                      end do
                 end do
              end do
              ijkl=ijkl+1
            end do  
          else
            do l=1,k
              do p=1,Nb
                 do q=1,Nb
                   pq=indexs(p,q)
                      do r=1,Nb
                        do s=1,Nb
                           rs=indexs(r,s)
                           pqrs=indexs(pq,rs)
                           twoe_mo(ijkl)=twoe_mo(ijkl)+C(p,i)*C(q,j)*C(r,k)*C(s,l)*twoe(pqrs)
                        end do
                      end do
                 end do
              end do
             ijkl=ijkl+1
            end do
          end if
         end do    
      end do
  end do
call cpu_time(time_2)
print*,"flag2=",time_2-time_1
return
end if

if(flag==3)then   !right for ERI_MO format 
call io_assign(lun)
!open(unit=lun,file='MO.dat')
call cpu_time(time_1)
 do p=1,Nb
   do q=1,Nb
       pq=indexs(p,q)
       do r=1,Nb
          do s=1,Nb
               rs=indexs(r,s)
               pqrs=indexs(pq,rs)
               ERI_MO(p,q,r,s)=0.d0
                do i=1,Nb
                    do j=1,Nb
                        ij=indexs(i,j)
                       do k=1,Nb
                          do l=1,Nb
                              kl=indexs(k,l)
                              ijkl=indexs(ij,kl)
!                             ERI_MO(p,q,r,s)= ERI_MO(p,q,r,s)+C(i,p)*C(j,q)*ERI_AO(i,j,k,l)&
                              ERI_MO(p,q,r,s)= ERI_MO(p,q,r,s)+C(i,p)*C(j,q)*twoe(ijkl)&
                                                    *C(k,r)*C(l,s)
                          end do
                       end do
                    end do
                end do
               twoe_mo(pqrs)=ERI_MO(p,q,r,s)
!               write(lun,fmt="(4i3,f15.8)")p,q,r,s,ERI_MO(p,q,r,s)
          end do
        end do
    end do
 end do
!call io_close(lun)
call cpu_time(time_2)
print*,"flag3=",time_2-time_1
  if(IOP>=5)then
  call io_assign(lun)
  open(unit=lun,file='ERI_MO.dat')
        do p=1,Nb
            do q=1,p
                pq=indexs(p,q)
                do r=1,Nb
                   do s=1,r
                      rs=indexs(r,s)
                      pqrs=indexs(pq,rs)
                      write(unit=lun,fmt="(7i4,f12.8)") p,q,r,s,pq,rs,pqrs,ERI_mo(p,q,r,s)
                   end do
                end do
            end do
        end do
   call io_close(lun)
   end if
return  
end if

end subroutine


subroutine spinmo(twoe_mo,Nmo,IOP)
use FCS        ,only  : ansym_ERI
! Nmo  Number of spin molecular obitals
! No   Number of occupied spin obritals
implicit none
integer               :: Nmo,No
real(dp)              :: time_1,time_2,value1,value2,twoe_mo(*)
integer               :: flag,Nb,lun,IOP
integer,intrinsic     :: ceiling
call io_assign(lun)
!open(unit=lun,file='Ansym.dat')
!open(unit=lun,file='./cmp/dat',status='old')
do p=1,Nmo
  do q=1,Nmo
     do r=1,Nmo
           pr=indexs(ceiling(p/2.0),ceiling(r/2.0))
           qr=indexs(ceiling(q/2.0),ceiling(r/2.0))
        do s=1,Nmo
           !pr=indexs(ceiling(p/2.0),ceiling(r/2.0))
           qs=indexs(ceiling(q/2.0),ceiling(s/2.0))
           prqs=indexs(pr,qs)
           value1=twoe_mo(prqs)*ones(p,r)*ones(q,s)
           ps=indexs(ceiling(p/2.0),ceiling(s/2.0))
           !qr=indexs(ceiling(q/2.0),ceiling(r/2.0))
           psqr=indexs(ps,qr)
           value2=twoe_mo(psqr)*ones(p,s)*ones(q,r)
           ansym_ERI(p,q,r,s)=value1-value2
!           write(lun,fmt="(4i3,f20.15)")p,q,r,s,ansym_ERI(p,q,r,s)
!           read(lun,fmt="(4i3,f20.15)")i,j,k,l,ansym_ERI(p,q,r,s)
        end do
     end do
  end do
end do
end subroutine  spinmo

integer function ones(a,b)
implicit none
integer::a,b
if(mod(a,2)==mod(b,2))then
   ones=1
else
   ones=0
end if
return
end function

subroutine spinFockEig(H,C,Nb,nmo,no,Eig,IOP)
! form spin-orbital Fock Matrix !
use FCS        ,only  : Spin_Fock,ansym_ERI,H_mo
implicit none
real(dp)              :: time_1,time_2
integer               :: flag,Nb,lun,nmo,no,IOP,m
real(dp)              :: H(Nb,Nb),C(Nb,Nb),Eig(Nb)
real(dp),intrinsic    :: matmul,transpose
integer,intrinsic     :: ceiling
character(len=13)     :: fmts
H_mo=0.d0
H_mo=matmul(transpose(C),H)
H_mo=matmul(H_mo,C)
  !equivalent to above
  !H=matmul(matmul(transpose(C),H),C)
  !H=matmul(matmul(C,H),transpose(C))
 if(IOP>=5)then
   print*,"Nmo=",Nmo
   print*,"No=",No
 end if
do p=1,nmo
  do q=1,nmo
     Spin_Fock(p,q)=H_mo(ceiling(p/2.0),ceiling(q/2.0))*ones(p,q)
     do m=1,no
        Spin_Fock(p,q)=Spin_Fock(p,q)+ansym_ERI(p,m,q,m) 
     end do
  end do
end do
!Spin_Fock=0.d0
!do p=1,nmo
!   do q=1,nmo
!       Spin_Fock(p,p)=Eig(ceiling(p/2.0))
!end do
      
!__________________output Spin_Fock______________
if(IOP>=5)then
  print*,"Fock matrix :"
  write(fmts,"(a,i2,7a)")"(",Nmo,'f7.3)'
  do i = 1,Nmo
    write(*,fmt=fmts)(Spin_Fock(i,j),j=1,Nmo)
  end do
end if
return
end subroutine spinFockEig

end module 
