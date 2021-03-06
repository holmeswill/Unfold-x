
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING in Elk distribution for license details.

!BOP
! !ROUTINE: vecfbz
! !INTERFACE:
subroutine vecfbz(eps,bvec,vpl,iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps  : zero component tolerance (in,real)
!   bvec : reciprocal lattice vectors (in,real(3,3))
!   vpl  : input vector in lattice coordinates (inout,real(3))
!   iv   : integer parts of vpl (out,integer(3))
! !DESCRIPTION:
!   Maps a vector in lattice coordinates to the first Brillouin zone. This is
!   done by first removing its integer components and then adding primitive
!   reciprocal lattice vectors until the shortest vector is found.
!
! !REVISION HISTORY:
!   Created September 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(in) :: bvec(3,3)
real(8), intent(inout) :: vpl(3)
integer, intent(out) :: iv(3)
! local variables
integer i1,i2,i3,j1,j2,j3
real(8) v0(3),v1(3),v2(3),v3(3),t1,t2
! map vector to [0,1) interval
call r3frac(eps,vpl,iv)
v0(:)=bvec(:,1)*vpl(1)+bvec(:,2)*vpl(2)+bvec(:,3)*vpl(3)
t1=v0(1)**2+v0(2)**2+v0(3)**2
j1=0;j2=0;j3=0
do i1=-1,0
  v1(:)=v0(:)+dble(i1)*bvec(:,1)
  do i2=-1,0
    v2(:)=v1(:)+dble(i2)*bvec(:,2)
    do i3=-1,0
      v3(:)=v2(:)+dble(i3)*bvec(:,3)
      t2=v3(1)**2+v3(2)**2+v3(3)**2
      if (t2.lt.t1+eps) then
        j1=i1;j2=i2;j3=i3
        t1=t2
      end if
    end do
  end do
end do
vpl(1)=vpl(1)+dble(j1)
vpl(2)=vpl(2)+dble(j2)
vpl(3)=vpl(3)+dble(j3)
iv(1)=iv(1)-j1
iv(2)=iv(2)-j2
iv(3)=iv(3)-j3
return
end subroutine
!EOC

! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3frac
! !INTERFACE:
subroutine r3frac(eps,v,iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
!   iv  : integer parts of v (out,integer(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!   The integer components of {\tt v} are returned in the variable {\tt iv}.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
integer, intent(out) :: iv(3)
! local variables
integer i
do i=1,3
  iv(i)=floor(v(i))
  v(i)=v(i)-iv(i)
  if ((1.d0-v(i)).lt.eps) then
    v(i)=0.d0
    iv(i)=iv(i)+1
  end if
  if (v(i).lt.eps) v(i)=0.d0
end do
return
end subroutine
!EOC
