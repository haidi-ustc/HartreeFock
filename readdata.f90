subroutine readdata()
use m_precision,only : dp
use m_debug
use FCS        
use calc       ,only : indexs
implicit none

integer         :: i,j,m,n,ij,mn,ijmn,istat,fileid
real(dp)        :: val,Temp(2000)
!for debug
integer         :: counts,flag,lun
flag=1

!------------read overlap Matrix-----------------
call io_close(fileid)
open(unit=fileid,file="./data/s.dat",status='old')
do while(.true.)
read(fileid,*,iostat=istat)i,j,val
 if(istat .ne. 0) exit
 s(i,j)=val
 s(j,i)=val
end do
call io_close(fileid)

!------------read V_core Matrix------------------
open(unit=fileid,file="./data/v.dat",status='old')
do while(.true.)
read(fileid,*,iostat=istat)i,j,val
 if(istat .ne. 0) exit
 V(i,j)=val
 V(j,i)=V(i,j)
end do
call io_close(fileid)

!------------read T_core Matrix------------------
open(unit=fileid,file="./data/t.dat",status='old')
do while(.true.)
read(fileid,*,iostat=istat)i,j,val
 if(istat .ne. 0) exit
 T(i,j)=val
 T(j,i)=T(i,j)
end do
call io_close(fileid)

!------------read ERI Matrix------------------
open(unit=fileid,file="./data/eri.dat",status='old')
if(IOP>=5)then
   print*,"  I   J    M    N    IJ   MN   IJMN"
end if

call io_assign(lun)

open(unit=lun,file="twoe.dat")
do while(.true.)
read(fileid,*,iostat=istat)i,j,m,n,val
             if(istat .ne. 0) exit
             IJ=indexs(i,j)
             mn=indexs(m,n)
             ijmn=indexs(ij,mn)
             if(IOP>=5)then
               write(*,"(6i5,i7)") i,j,m,n,IJ,mn,ijmn
             end if
             twoe(ijmn)=val

             write(lun,fmt="(i3,f15.8)")ijmn,val

if(flag==3)then
 ERI_AO(i,j,m,n)=val
 ERI_AO(j,i,m,n)=val
 ERI_AO(i,j,n,m)=val
 ERI_AO(j,i,n,m)=val
 ERI_AO(m,n,i,j)=val
 ERI_AO(m,n,j,i)=val
 ERI_AO(n,m,i,j)=val
 ERI_AO(n,m,j,i)=val
end if
end do
call io_close(lun)
call io_close(fileid)

!------- read enuc------------------------
open(unit=fileid,file="./data/enuc.dat",status='old')
read(fileid,*,iostat=istat)val
E_nuc=val
call io_close(fileid)
H_core=T+V

end subroutine
