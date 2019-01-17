



subroutine set_nlpp_parameters( nfile, myid, nodes,  &
& lrela, ntype, lchk, lvandi, lvflag, lvand, mxl, nylmmx_, mxref, lmax, nrefe )
!-----------------------------------------------------------------------
!     set parameters in nlpp_parameters
!-----------------------------------------------------------------------
use outfile
use psvariables
use nlpp_parameters
implicit none
integer :: nfile(*), myid, nodes
logical :: lrela
integer :: ntype, mxl, nylmmx_, mxref, lmax(ntype), nrefe(0:mxl,ntype)
logical :: lchk(0:mxl,ntype), lvandi(ntype), lvflag(ntype), lvand

!-----declare local variables
integer :: it, l, lmxx, lmpxx, lmptxx, laxx, lfmxx, lftmxx, i, j
integer :: lnlnl, ni, nj, ii, l1, l2, ll, j2


if( lrela ) then
    nylmmx = 2*(mxl+1)*(mxl+1) - 1
else
    nylmmx = mxl*(mxl+2) + 1
end if
nylmmx_ = nylmmx
nqlmmx = (mxl + 1)*(mxl + 2)/2
nrefxx = mxref
ncmxx  = nrefxx*(nrefxx-1)/2 + nrefxx
ncmtxx = nrefxx*nrefxx
!lmmax  = nrefxx*nylmmx
!lpmax  = ncmxx *nylmmx
!lptmax = ncmtxx*nylmmx
!lfmax  = lmmax*(lmmax + 1)/2
!lftmax = lmmax*lmmax
!lnref  = nrefxx*(mxl + 1)

lmmax  = 0
lpmax  = 0
lptmax = 0
lfmax  = 0
lftmax = 0
lnref  = 0

typedo: do it = 1, ntype
!typeif: if( lvandi(it) ) then
!--- pairs for  l.eq.l' .and.  m.eq.m' ---
   lmxx   = 0
   lmpxx  = 0
   lmptxx = 0
   laxx   = 0
   do l = 0, lmax(it)
      if( lchk(l,it) ) then
          do j = MJTAB1(l,it), MJTAB2(l,it)
             lmxx   = lmxx   + ( 2*l+1 ) * NREFE(l,it)
             lmpxx  = lmpxx  + ( 2*l+1 )  &
&                   * ( NREFE(l,it)*( NREFE(l,it) - 1 )/2 + NREFE(l,it) )
             lmptxx = lmptxx + ( 2*l+1 ) * NREFE(l,it) * NREFE(l,it)
             laxx = laxx + NREFE(l,it)
          end do
      end if
   end do
   lmmax  = max( lmmax,  lmxx )
   lpmax  = max( lpmax,  lmpxx )
   lptmax = max( lptmax, lmptxx )
   lnref  = max( lnref,  laxx )

!   lmx(it)   = lmxx
!   lpmx(it)  = lmpxx
!   lptmx(it) = lmptxx
!   lax(it)   = laxx

!--- pairs for  l.eq.l' .and.  m.ne.m' .or. l.ne.l' ---
   lfmxx  = lmxx*( lmxx + 1 )/2
   lftmxx = lmxx*lmxx
   lfmax  = max( lfmax,  lfmxx )
   lftmax = max( lftmax, lftmxx )

!   lfmx(it)  = lfmxx
!   lftmx(it) = lftmxx

!end if typeif
end do typedo

!do i = 1, 2
!   if( loutfile(i) ) then
!       write(nfile(i),*) 'lmmax  =', lmmax,  ' cf. ', nrefxx*nylmmx
!       write(nfile(i),*) 'lpmax  =', lpmax,  ' cf. ', ncmxx *nylmmx
!       write(nfile(i),*) 'lptmax =', lptmax, ' cf. ', ncmtxx*nylmmx
!       write(nfile(i),*) 'lfmax  =', lfmax,  ' cf. ', nrefxx*nylmmx*(nrefxx*nylmmx + 1)/2
!       write(nfile(i),*) 'lftmax =', lftmax, ' cf. ', nrefxx*nylmmx*nrefxx*nylmmx
!       write(nfile(i),*) 'lnref  =', lnref,  ' cf. ', nrefxx*(mxl + 1)
!   end if
!end do


if( lvand ) then

MXL2  = 2*mxl
MXL21 = MXL2 + 1
nqylmx = (MXL21 + 1)*(MXL21 + 2)/2
nqabm1 = ncmtxx*MXL*(MXL + 1)*(MXL + 2)/6
nqabm2 = ncmxx*(MXL + 1)*(MXL + 2)/2
!nqabmx = nqabm1 + nqabm2

nqabmx = 0

typedo2: do it = 1, ntype
typeif2: if( lvandi(it) ) then

   lnlnl = 0
   do l = 0, lmax(it)
   if( lchk(l,it) ) then
   do j = MJTAB1(l,it), MJTAB2(l,it)
       do NI =  1, NREFE(l,it)
       do NJ = NI, NREFE(l,it)
       do ll = 0, 2*l
          if( mod(l+l+ll,2).eq.0 ) then
              lnlnl = lnlnl + 1
          end if
       end do
       end do
       end do
    end do
    end if
    end do

    !   for full calculation
    if( lvflag(it) ) then
        !--- loop of l1 ---
        do l1 = 0, lmax(it)
        if( lchk(l1,it) ) then
        do j = MJTAB1(l1,it), MJTAB2(l1,it)
            DO NI = 1, NREFE(l1,it)
            !--- loop of l2 ---
!            do l2 = l1 + 1, lmax(it)
            do l2 = l1 , lmax(it)
            if( lchk(l2,it) ) then
            do j2 = MJTAB1(l2,it), MJTAB2(l2,it)
               if( l1 == l2 .and. j >= j2 ) cycle
               do NJ = 1, NREFE(l2,it)
                  do ll = abs(l1 - l2), l1 + l2
                     if( mod(l1+l2+ll,2).eq.0 ) then
                         lnlnl = lnlnl + 1
                     end if
                  end do
               end do
            end do
            end if
            end do
            !--- end of loop of l2 ---
            end do
        end do
        end if
        end do
        !--- end of loop of l1 ---
    end if
!   end of loop for full calculation

    nqabmx = max( nqabmx, lnlnl )

end if typeif2
end do typedo2

do i = 1, 2
   if( loutfile(i) ) then
       write(nfile(i),*) 'nqabmx =', nqabmx,  ' cf. ', nqabm1 + nqabm2
   end if
end do

else

!-----Q functions are not needed for kbpp
MXL2  = 0
MXL21 = 0
nqylmx = 0
nqabm1 = 0
nqabm2 = 0
nqabmx = 0

end if


return
end subroutine




subroutine sphrb2( l2, c1, s1, c2, s2, ylmr, ylmi, mxl )
!-----------------------------------------------------------------------
!     Spherical harmonics and its derivatives
!-----------------------------------------------------------------------
implicit real*8(a-h,o-z)
dimension ylmr(3,-mxl:mxl,0:mxl), ylmi(3,-mxl:mxl,0:mxl)
data pi / 3.141592653589793d0   /
data pi00 / 2.820947917738781d-01 /
data pi10 / 4.886025119029199d-01 /
data pi11 / 3.454941494713355d-01 /
data pi20 / 3.153915652525200d-01 /
data pi21 / 7.725484040463791d-01 /
data pi22 / 3.862742020231896d-01 /
data pi30 / 3.731763325901154d-01 /
data pi31 / 3.231801841141507d-01 /
data pi32 / 1.021985476433282d0   /
data pi33 / 4.172238236327841d-01 /

!      ylmr(1,0,0) = sqrt( 1.d0/(4.d0*pi) )
ylmr(1,0,0) = pi00
ylmi(1,0,0) = 0.d0
ylmr(2,0,0) = 0.d0
ylmi(2,0,0) = 0.d0
ylmr(3,0,0) = 0.d0
ylmi(3,0,0) = 0.d0
if( l2.le.0 ) return

!      ylmr(1,1, 0) = sqrt( 3.d0/(4.d0*pi) )*c1
ylmr(1, 0, 1) = pi10*c1
ylmi(1, 0, 1) = 0.d0
ylmr(2, 0, 1) = -pi10*s1
ylmi(2, 0, 1) = 0.d0
ylmr(3, 0, 1) = 0.d0
ylmi(3, 0, 1) = 0.d0
!      a1 = sqrt( 3.d0/(8.d0*pi) )*s1
a1  = pi11*s1
ac1 = pi11*c1
ylmr(1, 1, 1) = a1*c2
ylmi(1, 1, 1) = a1*s2
ylmr(2, 1, 1) = ac1*c2
ylmi(2, 1, 1) = ac1*s2
ylmr(3, 1, 1) = -pi11*s2
ylmi(3, 1, 1) =  pi11*c2
ylmr(1,-1, 1) = -ylmr(1, 1, 1)
ylmi(1,-1, 1) =  ylmi(1, 1, 1)
ylmr(2,-1, 1) = -ylmr(2, 1, 1)
ylmi(2,-1, 1) =  ylmi(2, 1, 1)
ylmr(3,-1, 1) = -ylmr(3, 1, 1)
ylmi(3,-1, 1) =  ylmi(3, 1, 1)
if( l2.le.1 ) return

c12 = c1*c1
s12 = s1*s1
c22p = c2*c2 - s2*s2
s22p = s2*c2 + c2*s2
!      ylmr(1, 0, 2) = sqrt( 5.d0/(16.d0*pi) )*(3.d0*c12 - 1.d0)
ylmr(1, 0, 2) = pi20*(3.d0*c12 - 1.d0)
ylmi(1, 0, 2) = 0.d0
ylmr(2, 0, 2) = -pi20*6.d0*c1*s1
ylmi(2, 0, 2) = 0.d0
ylmr(3, 0, 2) = 0.d0
ylmi(3, 0, 2) = 0.d0
!      a1 = sqrt( 15.d0/(8.d0*pi) )*s1*c1
a01 = pi21*c1
a1  =  a01*s1
ac1 = pi21*( c12 - s12 )
ylmr(1, 1, 2) = a1*c2
ylmi(1, 1, 2) = a1*s2
ylmr(2, 1, 2) = ac1*c2
ylmi(2, 1, 2) = ac1*s2
ylmr(3, 1, 2) = -a01*s2
ylmi(3, 1, 2) =  a01*c2
ylmr(1,-1, 2) = -ylmr(1, 1, 2)
ylmi(1,-1, 2) =  ylmi(1, 1, 2)
ylmr(2,-1, 2) = -ylmr(2, 1, 2)
ylmi(2,-1, 2) =  ylmi(2, 1, 2)
ylmr(3,-1, 2) = -ylmr(3, 1, 2)
ylmi(3,-1, 2) =  ylmi(3, 1, 2)
!      a1 = sqrt( 15.d0/(32.d0*pi) )*s12
a1  = pi22*s12
a01 = 2.d0*pi22*s1
ac1 =  a01*c1
ylmr(1, 2, 2) = a1*c22p
ylmi(1, 2, 2) = a1*s22p
ylmr(2, 2, 2) = ac1*c22p
ylmi(2, 2, 2) = ac1*s22p
ylmr(3, 2, 2) = -a01*s22p
ylmi(3, 2, 2) =  a01*c22p
ylmr(1,-2, 2) =  ylmr(1, 2, 2)
ylmi(1,-2, 2) = -ylmi(1, 2, 2)
ylmr(2,-2, 2) =  ylmr(2, 2, 2)
ylmi(2,-2, 2) = -ylmi(2, 2, 2)
ylmr(3,-2, 2) =  ylmr(3, 2, 2)
ylmi(3,-2, 2) = -ylmi(3, 2, 2)
if( l2.le.2 ) return


!      write(*,*) ' ***  l=3 is not supported yet .'
c13 = c12*c1
s13 = s12*s1
c23p = c22p*c2 - s22p*s2
s23p = s22p*c2 + c22p*s2
!      ylmr(1,0, 3) = sqrt( 7.d0/(16.d0*pi) )*c1*(5.d0*c12 - 3.d0)
ylmr(1, 0, 3) = pi30*c1*(5.d0*c12 - 3.d0)
ylmi(1, 0, 3) = 0.d0
ylmr(2, 0, 3) = -pi30*s1*( 15.d0*c12 - 3.d0 )
ylmi(2, 0, 3) = 0.d0
ylmr(3, 0, 3) = 0.d0
ylmi(3, 0, 3) = 0.d0
!      a1 = sqrt( 21.d0/(64.d0*pi) )*s1*(5.d0*c12 - 1.d0)
a01 = pi31*(5.d0*c12 - 1.d0)
a1  = a01*s1
ylmr(1, 1, 3) = a1*c2
ylmi(1, 1, 3) = a1*s2
ac1 = pi31*c1*(4.d0 - 15.d0*s12)
ylmr(2, 1, 3) = ac1*c2
ylmi(2, 1, 3) = ac1*s2
ylmr(3, 1, 3) = -a01*s2
ylmi(3, 1, 3) =  a01*c2
ylmr(1,-1, 3) = -ylmr(1, 1, 3)
ylmi(1,-1, 3) =  ylmi(1, 1, 3)
ylmr(2,-1, 3) = -ylmr(2, 1, 3)
ylmi(2,-1, 3) =  ylmi(2, 1, 3)
ylmr(3,-1, 3) = -ylmr(3, 1, 3)
ylmi(3,-1, 3) =  ylmi(3, 1, 3)
!      a1 = sqrt( 105.d0/(32.d0*pi) )*s12*c1
a1 = pi32*s12*c1
ylmr(1, 2, 3) = a1*c22p
ylmi(1, 2, 3) = a1*s22p
ac1 = pi32*s1*(2.d0*c12 - s12)
ylmr(2, 2, 3) = ac1*c22p
ylmi(2, 2, 3) = ac1*s22p
a01 = 2.d0*pi32*s1*c1
ylmr(3, 2, 3) = -a01*s22p
ylmi(3, 2, 3) =  a01*c22p
ylmr(1,-2, 3) =  ylmr(1, 2, 3)
ylmi(1,-2, 3) = -ylmi(1, 2, 3)
ylmr(2,-2, 3) =  ylmr(2, 2, 3)
ylmi(2,-2, 3) = -ylmi(2, 2, 3)
ylmr(3,-2, 3) =  ylmr(3, 2, 3)
ylmi(3,-2, 3) = -ylmi(3, 2, 3)
!      a1 = sqrt( 35.d0/(64.d0*pi) )*s13
a1 = pi33*s13
ylmr(1, 3, 3) = a1*c23p
ylmi(1, 3, 3) = a1*s23p
a01 = 3.d0*pi33*s12
ac1 = a01*c1
ylmr(2, 3, 3) = ac1*c23p
ylmi(2, 3, 3) = ac1*s23p
ylmr(3, 3, 3) = -a01*s23p
ylmi(3, 3, 3) =  a01*c23p
ylmr(1,-3, 3) = -ylmr(1, 3, 3)
ylmi(1,-3, 3) =  ylmi(1, 3, 3)
ylmr(2,-3, 3) = -ylmr(2, 3, 3)
ylmi(2,-3, 3) =  ylmi(2, 3, 3)
ylmr(3,-3, 3) = -ylmr(3, 3, 3)
ylmi(3,-3, 3) =  ylmi(3, 3, 3)
!cc      if( l2.le.3 ) return

return
end




