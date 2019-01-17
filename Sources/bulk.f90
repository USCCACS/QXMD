



module para_ewald
!-----------------------------------------------------------------------
! type declaration of variables for Ewald method
!-----------------------------------------------------------------------
implicit none


real*8 :: alloc_mem

real*8  :: gamma     ! Ewald parameter
real*8  :: dkmax2    ! maximum reciprocal vector

real*8  ::  pi

real*8  ::  volume
integer, allocatable, dimension(:) :: nga, ngb, ngc
real*8,  allocatable, dimension(:) :: gx, gy, gz
real*8,  allocatable, dimension(:) :: recnrm
real*8,  allocatable, dimension(:) :: recnrmex ! for the double-grid method
integer :: numk2
integer :: KCSMA1, KCSMA2, KCSMA3

real*8,  allocatable, dimension(:,:) :: DCQ1, DCQ2, DCQ3
real*8,  allocatable, dimension(:,:) :: DSQ1, DSQ2, DSQ3
real*8,  allocatable, dimension(:,:) :: DCm1, DCm2, DCm3
real*8,  allocatable, dimension(:,:) :: DSm1, DSm2, DSm3

real*8,  allocatable, dimension(:) :: sumtmp, buff

integer :: numk21, numk22, nnumk2

integer :: numk2nod1, numk2nod2, numk2nod
integer, allocatable, dimension(:) :: numk2ncnt, numk2ndsp

integer :: mxatom = 0
integer :: ifsetm = 0


save


end module




subroutine ewaldini( nfile, myid, nodes,  &
& ntype, zv, ecutdens, talloc_mem, mx1, llclpp_r, llclpp_g,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& rctflc, llking, rlking, glkgmax, glkgexct, hcell, h_MD, gamma_tr,  &
& rccc2, lhfull, vlocli ,xitgrd, rr, nvlcl, nion, mshx, mshy, mshz,  &
& ldouble_grid_recip, nd1v, nd1vks, dgalpha, pwscale, lvshape )
!-----------------------------------------------------------------------
!     initialize Ewald method for bulk calculations
!-----------------------------------------------------------------------
use outfile
use para_ewald
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
real*8,  dimension(ntype) :: zv
real*8 :: ecutdens
real*8 :: talloc_mem
real*8,  dimension(0:mx1) :: tberf, tberfa
real*8 :: drcut
real*8,  dimension(0:mx1) :: tbeff, tbeffa
real*8 :: drcutf
real*8,  dimension(0:mx1,ntype) :: tablc, tablca
real*8,  dimension(ntype)       :: dltlc, rmxlc
real*8,  dimension(0:mx1,ntype) :: tbflc, tbflca
real*8,  dimension(ntype)       :: dltflc, rmxflc
real*8 :: rctflc
logical :: llclpp_r, llclpp_g
logical, dimension(ntype) :: llking
real*8,  dimension(ntype) :: rlking, glkgmax, glkgexct
real*8  :: hcell(3,3), h_MD(3,3)
logical :: lhfull
real*8,  dimension(0:nvlcl) :: vlocli ,xitgrd, rr
integer :: nion
integer :: mshx, mshy, mshz
logical :: ldouble_grid_recip
integer, dimension(3) :: nd1v, nd1vks
real*8  :: dgalpha
real*8  :: pwscale
logical :: lvshape

real*8,  dimension(3,3) :: b
real*8,  dimension(3,3) :: rba

!-----variables for Ewald parameter
integer, parameter :: NAPMAX = 50
real*8,  dimension(NAPMAX) :: ALPX2, ALPX3, ALPXN
real*8  :: XBASE, XKANR
integer :: ifcall = 0
save ifcall, ALPX2, ALPX3, ALPXN, XBASE, XKANR


alloc_mem = talloc_mem

    ACCURR = 0.1000d-09
    ACCURY = 0.1000d-09

pi = acos(-1.d0)

!--- volume and reciprocal lattice vectors -----------------------------
    if( .not.ldouble_grid_recip ) then
        do i = 1, 3
           do ix = 1, 3
              rba(ix,i) = hcell(ix,i)
           end do
        end do
      else
        do i = 1, 3
           do ix = 1, 3
              rba(ix,i) = hcell(ix,i)*dble(nd1vks(i))/dble(nd1v(i))
           end do
        end do
    end if

call RCIPRL( rba, B, volume )


!--- cutoff length in real space 
    rccc1 = 1.d+10
    do i = 1, 3
       RK1Di = B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i)
       rk1di = 0.5d0/sqrt( rk1di )
       rccc1  = min( rccc1, rk1di )
    end do
    rccc = rccc1*0.95d0

!--- modify cutoff length for filtering local pseudopotentials
rlkgmax = -1.d0
do it = 1, ntype
   if( llking(it) ) rlkgmax = max( rlkgmax, rlking(it) )
end do
lhfull = .false.
if( rlkgmax.gt.1.d0 ) then
        if( rccc/rlkgmax.ge.rctflc ) rccc = rccc/rlkgmax
        lhfull = rccc*rlkgmax .gt. rccc1
end if


do i = 1, 2
if( loutfile(i) ) then
    write(nfile(i),*) '   Ewald:   rccc', rccc
    write(nfile(i),*) '   Ewald:   lhfull', lhfull
end if
end do
rccc2 = rccc*rccc


!---  determination of Ewald potential parameter ----------------------
if( ifcall .eq. 0 ) then
    ifcall = 1
    call SETALP( ALPX2, ALPX3, ALPXN, XBASE, XKANR, NAPMAX )
end if

ACRYLN = ACCURR*rccc
ACRYLN = - LOG(ACRYLN) - XBASE

    Dewa  = ACRYLN*XKANR
    Lewa  = Dewa
    DLewa = Lewa
    Dewa  = Dewa - DLewa

    Lewa  = Lewa + 1
    ALP = Dewa*((Dewa-1.0)*ALPX2(Lewa) + ALPX3(Lewa)) +ALPXN(Lewa)


gamma    = ALP/rccc
gamma_tr = gamma

do i = 1, 2
if( loutfile(i) ) then
    write(nfile(i),*) '   Ewald:   ALP  ', alp
    write(nfile(i),*) '   Ewald:   gamma', gamma
end if
end do
!-----------------------------------------------------------------------

!--- truncation length in reciprocal space -----------------------------
gamma2 = gamma*gamma
CNST1  = pi*pi/gamma2
CNST1R = 1.0D0/CNST1
VACURY = VOLUME*ACCURY

    DKMAX2 = - CNST1R*LOG(VACURY)

    ACLNPI = VACURY*DKMAX2
    DKMAX2 = - CNST1R*LOG(ACLNPI)

    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '   Ewald:   DKMAX2', DKMAX2
    end if
    end do


!--- set reciprocal vectors --------------------------------------------
call clmrcp( nfile, myid, nodes, b, nion, mshx, mshy, mshz,  &
& ldouble_grid_recip, dgalpha, llclpp_r, pwscale, lvshape )
!cc      call stgexp


!--- set error function table ------------------------------------------
call seterf(  &
& alp, tberf, tberfa, drcut, tbeff, tbeffa, drcutf, mx1 )


!--- filtering local pseudopotentials ----------------------------------
if( llclpp_r )  &
& call locking( nfile, myid, nodes,  &
& ntype, zv, ecutdens, gamma, rccc2, mx1,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& llking, rlking, glkgmax, glkgexct, vlocli ,xitgrd, rr, nvlcl )


talloc_mem = alloc_mem


return
end




subroutine locking( nfile, myid, nodes,  &
& ntype, zv, ecutdens, gamma, rccc2, mx1,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& llking, rlking, glkgmax, glkgexct, vlocli ,xitgrd, rr, nvlcl )
!-----------------------------------------------------------------------
!    filtering local pseudopotentials
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
real*8,  dimension(ntype) :: zv
real*8, dimension(0:mx1) :: tberf, tberfa
real*8 :: drcut
real*8, dimension(0:mx1) :: tbeff, tbeffa
real*8 :: drcutf
real*8,  dimension(0:mx1,ntype) :: tablc, tablca
real*8,  dimension(ntype)       :: dltlc, rmxlc
real*8,  dimension(0:mx1,ntype) :: tbflc, tbflca
real*8,  dimension(ntype)       :: dltflc, rmxflc
logical, dimension(ntype) :: llking
real*8,  dimension(ntype) :: rlking, glkgmax, glkgexct

real*8, dimension(0:nvlcl)  :: vlocli ,xitgrd, rr
!-----declare local variables
integer, parameter :: noofdt = 200
real*8, allocatable, dimension(:) :: betaq
real*8 :: gdel
real*8, allocatable, dimension(:) :: betar
real*8 :: rspdel
real*8,  allocatable, dimension(:,:) :: amatrx
real*8,  allocatable, dimension(:,:) :: vmatrx
integer, allocatable, dimension(:)   :: ivmatrx
integer :: status

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 1
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip



!------allocate memory
allocate( betaq(0:noofdt), betar(0:noofdt),  &
& amatrx(0:noofdt,0:noofdt), vmatrx(noofdt,noofdt+1),  &
& ivmatrx(noofdt), stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in locking' )

alocmem = 8.d0 * ( noofdt + 1 ) * ( 2 + noofdt + 1 + noofdt )  &
&       + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, allocated (locking)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, allocated (locking)'


rccc  = sqrt(rccc2)
deltr = rccc/dble(mx1)

gamma2 = gamma*gamma
!cc      rlkgmax = -1.d0
typedo: do it = 1, ntype
if( mod(it-1,nodes) == myid ) then
   zvit2  = 2.d0*zv(it)
   zvit2g = zvit2*gamma
   typeif: if( llking(it) ) then
!--- filtering ---------------------------------------------------------
!cc             rlkgmax = max( rlkgmax, rlking(it) )
       lltk = 0
       iitk = it

       xrr = 0.d0
       do ir = 0, mx1
          rr(ir) = xrr
          r = xrr
          r2 = r*r

          r2g = r2*gamma2
          m = r2g/(2.d0*drcut)
          m = 2*m
          d = 0.5d0*( r2g/drcut - dble(m) )
          cerf = d*( (d-1.d0)*tberfa(m) + tberf(m+2)  &
&                        - tberf(m)   ) + tberf(m)
          vlocef = zvit2g*cerf

          if( r.lt.rmxlc(it) ) then
              m = r/(2.d0*dltlc(it))
              m = 2*m
              d = 0.5d0*( r/dltlc(it) - dble(m) )
              vloc = d*( (d-1.d0)*tablca(m,it) + tablc(m+2,it)  &
&                              - tablc(m,it) ) + tablc(m,it)
          else
              vloc = - zvit2/r
          end if
          vlocli(ir) = (vlocef+vloc)*xrr*xrr
!      write(95,'(3e14.6)') r, vlocef+vloc
          xrr = xrr + deltr
       end do
       dkmax4 = sqrt(ecutdens) * glkgmax(it)
       dkmax  = sqrt(ecutdens) * glkgexct(it)
       eg4   = dkmax4*dkmax4
       egmax = dkmax*dkmax
       rmxnlc     = rccc * rlking(it)
       rmxlc(it)  = rmxnlc
       rmxflc(it) = rmxnlc
!            --- Gaussian filtering ------------------------------
       call kingsm( nfile, myid, nodes,  &
& lltk, iitk, deltr, egmax, eg4, .false., rmxnlc,  &
& vlocli ,xitgrd, rr, nvlcl, betaq, gdel, betar, rspdel, noofdt,  &
& amatrx, vmatrx, ivmatrx )
!            -----------------------------------------------------
       dltlc(it)  = rmxlc(it)*rmxlc(it)/dble(mx1-2)
       dltflc(it) = rmxflc(it)*rmxflc(it)/dble(mx1-2)

!        --- interpolation by  spline ---
       LL   = 5
       N0   = 30
       NOVERH = LL*2
       NOVER  = NOVERH*2
       INI  = 0
       ILSM = 0
       xrr2 = 0.d0
       do ir = 1, mx1 - 2
          xrr2 = xrr2 + dltlc(it)
          xrr = sqrt(xrr2)
          if( xrr > rspdel*dble(ILSM) ) then
              do
                 ILS = MIN( INI + N0 - 1, noofdt )
                 IF( ILS-INI+1 < LL+2 ) LL = ILS - INI - 1
                 IF( ILS == noofdt ) THEN
                     ILSM = noofdt + 1
                   ELSE
                     ILSM = ILS - NOVERH
                 end if
                 IF( rspdel*dble(ILSM) >= xrr ) exit
                 INI = INI + NOVERH
              end do

              Nm = 0
              do LP = INI, ILS
                 Nm = Nm + 1
                 xs(Nm) = rspdel*dble(LP)
                 ys(Nm) = betar(LP)
              end do
              CALL SPLINS( LL, nm, XS, YS, QS, BS, AS, WS, IP,  &
&                          KQM, MM, MM1, NDV )
              INI  = ILS - NOVER
          end if

          CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS,  &
&                      KQM, MM, NDV, 0 )
          CALL SPLINV( LL, Nm, xrr, ppfsv, QS, AS, WS,  &
&                      KQM, MM, NDV, 1 )
          tablc(ir,it) = ppsv
          tbflc(ir,it) = ppfsv/xrr
!      write(96,'(3e14.6)') xrr, ppsv
       end do



       tablc(0,it)     = betar(0)
       tablc(mx1-2,it) = betar(noofdt)
       tablc(mx1-1,it) = 0.d0
       tablc(mx1  ,it) = 0.d0
       tbflc(0,it)     = 2.d0*tbflc(1,it) - tbflc(2,it)
       tbflc(mx1-2,it) = tbflc(mx1-3,it)
       tbflc(mx1-1,it) = 0.d0
       tbflc(mx1,it)   = 0.d0

   end if typeif
end if
end do typedo


do it = 1, ntype
   iroot = mod(it-1,nodes)
   if( llking(it) ) then

   !--- unify variables by broadcast
       call dbcast( rmxlc(it),1,iroot)
       call dbcast(rmxflc(it),1,iroot)
       call dbcast( dltlc(it),1,iroot)
       call dbcast(dltflc(it),1,iroot)
       call dbcast(tablc(0,it),mx1+1,iroot)
       call dbcast(tbflc(0,it),mx1+1,iroot)

       do ir = 0, mx1 - 2
          tablca(ir,it) = 2.D0*( tablc(ir,it) + tablc(ir+2,it)  &
&                                        - 2.d0*tablc(ir+1,it) )
          tbflca(ir,it) = 2.D0*( tbflc(ir,it) + tbflc(ir+2,it)  &
&                                        - 2.d0*tbflc(ir+1,it) )
       end do

   end if
end do


!------deallocate memory
deallocate( betaq, betar, amatrx, vmatrx, ivmatrx, stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in locking' )

alocmem = 8.d0 * ( noofdt + 1 ) * ( 2 + noofdt + 1 + noofdt )  &
&       + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, deallocated (locking)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, deallocated (locking)'


return
end




SUBROUTINE SETALP( ALPX2, ALPX3, ALPXN, XBASE, XKANR, NAPMAX )
!-----------------------------------------------------------------------
!             set parameters for Ewald potential
!-----------------------------------------------------------------------
implicit none

integer :: NAPMAX
real*8, dimension(NAPMAX) :: ALPX2, ALPX3, ALPXN
real*8 :: XBASE, XKANR
real*8 :: XKAN

integer, parameter :: NALMAX = 101
real*8,  dimension(NALMAX) :: ALPX =  &
&(/0.271476099D+01,0.274928930D+01,0.278343255D+01,0.281720299D+01, &
& 0.285061226D+01,0.288367138D+01,0.291639088D+01,0.294878072D+01, &
& 0.298085044D+01,0.301260911D+01,0.304406537D+01,0.307522752D+01, &
& 0.310610344D+01,0.313670070D+01,0.316702654D+01,0.319708789D+01, &
& 0.322689139D+01,0.325644342D+01,0.328575011D+01,0.331481734D+01, &
& 0.334365074D+01,0.337225576D+01,0.340063763D+01,0.342880138D+01, &
& 0.345675187D+01,0.348449376D+01,0.351203157D+01,0.353936965D+01, &
& 0.356651220D+01,0.359346327D+01,0.362022678D+01,0.364680652D+01, &
& 0.367320616D+01,0.369942924D+01,0.372547919D+01,0.375135934D+01, &
& 0.377707290D+01,0.380262300D+01,0.382801264D+01,0.385324479D+01, &
& 0.387832226D+01,0.390324784D+01,0.392802419D+01,0.395265393D+01, &
& 0.397713958D+01,0.400148361D+01,0.402568838D+01,0.404975625D+01, &
& 0.407368948D+01,0.409749024D+01,0.412116070D+01,0.414470290D+01, &
& 0.416811893D+01,0.419141075D+01,0.421458022D+01,0.423762931D+01, &
& 0.426055980D+01,0.428337352D+01,0.430607216D+01,0.432865741D+01, &
& 0.435113099D+01,0.437349446D+01,0.439574942D+01,0.441789757D+01, &
& 0.443994010D+01,0.446187885D+01,0.448371437D+01,0.450544960D+01, &
& 0.452708473D+01,0.454862134D+01,0.457006072D+01,0.459140502D+01, &
& 0.461265431D+01,0.463380859D+01,0.465487444D+01,0.467584658D+01, &
& 0.469673095D+01,0.471752296D+01,0.473823020D+01,0.475885131D+01, &
& 0.477938713D+01,0.479983802D+01,0.482020507D+01,0.484048788D+01, &
& 0.486070312D+01,0.488081054D+01,0.490085693D+01,0.492081970D+01, &
& 0.494072265D+01,0.496052734D+01,0.498026260D+01,0.499992187D+01, &
& 0.501950805D+01,0.503903320D+01,0.505843749D+01,0.507789062D+01, &
& 0.509716796D+01,0.511642578D+01,0.513552734D+01,0.515472656D+01, &
& 0.517364746D+01/)
save ALPX

integer :: i, j


 XBASE = 9.0D0
 XKAN  = 0.4D0
 XKANR = 1.0D0/XKAN

 J = 1
 DO I = 1, NAPMAX
    ALPXN(I) = ALPX(J)
    ALPX2(I) = 2.0D0*(ALPX(J+2) + ALPX(J)) - 4.0D0*ALPX(J+1)
    ALPX3(I) = ALPX(J+2) - ALPX(J)

    J = J + 2
 end do


RETURN
END subroutine




module for_clmrcp
!-----------------------------------------------------------------------
! type declaration of allocatable variables in clmrcp
!-----------------------------------------------------------------------
implicit none

integer :: ifftk0 = 0
integer :: lplnumk2 = 0
integer :: lplKCSMA1 = 0
integer :: lplKCSMA2 = 0
integer :: lplKCSMA3 = 0
integer :: mplKCSMA1 = 0
integer :: mplKCSMA2 = 0
integer :: mplKCSMA3 = 0
integer :: mplmshx = 0
integer :: mplmshy = 0
integer :: mplmshz = 0
integer, allocatable, dimension(:) :: ktmp1, ktmp2, ktmp3
real*8,  allocatable, dimension(:) :: gss

!--- for sorting ------------------------------
integer, allocatable, dimension(:) :: key, ix, dstkey, x1, x2
real*8 :: dkeymx = 1000000000.d0
save ifftk0, lplnumk2, lplKCSMA1, lplKCSMA2, lplKCSMA3, dkeymx
save mplKCSMA1, mplKCSMA2, mplKCSMA3, mplmshx, mplmshy, mplmshz

end module




subroutine clmrcp( nfile, myid, nodes, b, nion_nod,  &
& mshx, mshy, mshz, ldouble_grid_recip, dgalpha, llclpp_r, pwscale,  &
& lvshape )
!-----------------------------------------------------------------------
!       determination of range of reciprocal lattice vectors
!                                         for Ewald potential
!-----------------------------------------------------------------------
use outfile
use para_ewald
use for_clmrcp
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
real*8, dimension(3,3) :: b
integer :: nion_nod
integer :: mshx, mshy, mshz
logical :: ldouble_grid_recip
real*8  :: dgalpha
logical :: llclpp_r
real*8  :: pwscale
logical :: lvshape

!------declare local variables
real*8, dimension(3)   :: rk1d
real*8, dimension(3,3) :: cross
integer :: mshxx, mshyx, mshzx
integer :: status
integer :: nion
integer, dimension(3) :: IIP1 = (/ 2, 3, 1 /)
integer, dimension(3) :: IIP2 = (/ 3, 1, 2 /)
integer, dimension(3) :: IIP3 = (/ 1, 2, 3 /)
integer, dimension(3) :: kmaxm1
integer, dimension(24) :: tmp1, tmp2, tmp3
integer :: NUMg
save IIP1, IIP2, IIP3


!--- determination of maximum integer : KMAXM --------------------------
do i = 1, 3
   RK1D(i) = B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i)
end do

KMAXM = 0
do i = 1, 3
   IP1 = IIP1(i)
   IP2 = IIP2(i)
   IP3 = IIP3(i)
   A1X = B(2,IP1)*B(3,IP2) - B(3,IP1)*B(2,IP2)
   A1Y = B(3,IP1)*B(1,IP2) - B(1,IP1)*B(3,IP2)
   A1Z = B(1,IP1)*B(2,IP2) - B(2,IP1)*B(1,IP2)
   AAA = A1X*A1X + A1Y*A1Y + A1Z*A1Z
   AAA = 1.0D0/SQRT(AAA)
   A1X = A1X*AAA
   A1Y = A1Y*AAA
   A1Z = A1Z*AAA
   PK3 = B(1,IP3)*A1X + B(2,IP3)*A1Y + B(3,IP3)*A1Z

       DMMZ  = sqrt(dkmax2)/PK3
       kmaxm1(i) = int(ABS(DMMZ))
   KMAXM = max( KMAXM, kmaxm1(i) )
end do

if(loutfile(1)) write(nfile(1),*) '   Ewald:   KMAXM', KMAXM
if(loutfile(2)) write(nfile(2),*) '   Ewald:   KMAXM', KMAXM
!-----------------------------------------------------------------------


!      ifftk = 2*(KMAXM + 1)
!      ifftk = ifftk*ifftk*ifftk/2

ifftk1 = 2*(kmaxm1(1) + 2)
ifftk2 = 2*(kmaxm1(2) + 2)
ifftk3 = 2*(kmaxm1(3) + 2)
ifftk = ifftk1*ifftk2*ifftk3/2

!-----check array sizes
if( ifftk > ifftk0 ) then

    !-----if already allocated, deallocate arrays
    if( allocated(gss) ) then
        deallocate( gss, ktmp1, ktmp2, ktmp3, key,  &
& ix, dstkey, x1, x2, stat=status )
        dealocmem = 8.d0 * ifftk0 * 1.d0 + 4.d0 * ifftk0 * 8.d0
    if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d6),' MB, deallocated (clmrcp)'
    if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d6),' MB, deallocated (clmrcp)'
    end if

    !------allocate memory
    ifftk0 = ifftk * pwscale
    allocate( gss(ifftk0),  &
& ktmp1(ifftk0), ktmp2(ifftk0), ktmp3(ifftk0), key(ifftk0),  &
& ix(ifftk0), dstkey(ifftk0), x1(ifftk0), x2(ifftk0),  &
& stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in clmrcp' )

    alocmem = 8.d0 * ifftk0 * 1.d0 + 4.d0 * ifftk0 * 8.d0
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d6),' MB, allocated (clmrcp)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d6),' MB, allocated (clmrcp)'

end if


numg = 0
do i = 2, KMAXM + 1
   ih  = i - 1
   do j = 1, i
      jh  = j - 1
      do k = 1, j
         kh  = k - 1

         num = 0
         kji = 1
         if( k == 1 ) kji = 2
         if( k == j ) kji = kji + 2
         if( i == j ) kji = kji + 4

         num = num + 1
         tmp1(num) =  ih
         tmp2(num) =  jh
         tmp3(num) =  kh
         if( kji /= 7 ) then
             num = num + 1
             tmp1(num) =  kh
             tmp2(num) =  ih
             tmp3(num) =  jh
             num = num + 1
             tmp1(num) =  jh
             tmp2(num) =  kh
             tmp3(num) =  ih
         end if
         if( kji /= 4 ) then
             num = num + 1
             tmp1(num) = -ih 
             tmp2(num) =  jh 
             tmp3(num) =  kh
             num = num + 1
             tmp1(num) =  kh 
             tmp2(num) = -ih 
             tmp3(num) =  jh
             num = num + 1
             tmp1(num) =  jh 
             tmp2(num) =  kh 
             tmp3(num) = -ih
             if( kji /= 7 .and. kji /= 6 ) then
             if( kji /= 2 ) then
                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) = -jh 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  ih 
                 tmp3(num) = -jh
                 num = num + 1
                 tmp1(num) = -jh 
                 tmp2(num) =  kh 
                 tmp3(num) =  ih

                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) =  jh 
                 tmp3(num) = -kh
                 num = num + 1
                 tmp1(num) = -kh 
                 tmp2(num) =  ih 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) = -kh 
                 tmp3(num) =  ih
             end if
             if( kji /= 3 .and. kji /= 5 ) then
                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) =  kh 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) =  ih 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  jh 
                 tmp3(num) =  ih

                 num = num + 1
                 tmp1(num) = -ih 
                 tmp2(num) =  kh 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) = -ih 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  jh 
                 tmp3(num) = -ih
                 if( kji /= 2 ) then
                     num = num + 1
                     tmp1(num) =  ih 
                     tmp2(num) =  kh 
                     tmp3(num) = -jh
                     num = num + 1
                     tmp1(num) = -jh 
                     tmp2(num) =  ih 
                     tmp3(num) =  kh
                     num = num + 1
                     tmp1(num) =  kh 
                     tmp2(num) = -jh 
                     tmp3(num) =  ih

                     num = num + 1
                     tmp1(num) =  ih 
                     tmp2(num) = -kh 
                     tmp3(num) =  jh
                     num = num + 1
                     tmp1(num) =  jh 
                     tmp2(num) =  ih 
                     tmp3(num) = -kh
                     num = num + 1
                     tmp1(num) = -kh 
                     tmp2(num) =  jh 
                     tmp3(num) =  ih
                 end if
             end if
             end if
         end if
         do mkk = 1, num
            if( abs(tmp1(mkk)) .le. kmaxm1(1)+1 .and.  &
&               abs(tmp2(mkk)) .le. kmaxm1(2)+1 .and.  &
&               abs(tmp3(mkk)) .le. kmaxm1(3)+1      ) then
                numg = numg + 1
                ktmp1(numg) = tmp1(mkk) 
                ktmp2(numg) = tmp2(mkk) 
                ktmp3(numg) = tmp3(mkk) 
            end if
         end do
      end do
   end do
end do
num = numg
if(loutfile(1)) write(nfile(1),*) '   Ewald:   num =', num
if(loutfile(2)) write(nfile(2),*) '   Ewald:   num =', num
!-----------------------------------
!      dpi = acos(-1.d0)*2.d0
do i=1, 3
do j=1, 3
!          cross(j,i) = dpi*B(j,i)
    cross(j,i) = B(j,i)
end do
end do
do mkk = 1, num
   d1h = ktmp1(mkk)
   d2h = ktmp2(mkk)
   d3h = ktmp3(mkk)
   gs1mkk = cross(1,1)*d1h + cross(1,2)*d2h + cross(1,3)*d3h
   gs2mkk = cross(2,1)*d1h + cross(2,2)*d2h + cross(2,3)*d3h
   gs3mkk = cross(3,1)*d1h + cross(3,2)*d2h + cross(3,3)*d3h
!*                ( gs1, gs2, gs3 ) is reciprocal lattice vector

   gss(mkk) = gs1mkk*gs1mkk + gs2mkk*gs2mkk + gs3mkk*gs3mkk
end do


!--- sorting by length of vector : gss --------------------------------
!     ix     : data pointer of record

ddmx = 0.d0
do mkk = 1, num
   ddmx = max( ddmx, abs(gss(mkk)) )
end do

fctmax = dkeymx/ddmx
do mkk = 1, num
   key(mkk) = fctmax*gss(mkk)
   ix(mkk)  = mkk
end do

call vsgar( num, key, ix, ifftk, dstkey, x1, x2, x2 )
call vsoe(  num, key, ix, ifftk )
!-----------------------------------------------------------------------


!-----count number of reciprocal vectors
numk2 = 0
do i = 1, num
   if( gss(i).le.dkmax2 ) numk2 = numk2 + 1
end do
if(loutfile(1)) write(nfile(1),*) '   Ewald:   numk2 =', numk2
if(loutfile(2)) write(nfile(2),*) '   Ewald:   numk2 =', numk2


!-----allocate memory for reciprocal vectors--------------------------
!-----check array sizes
if( numk2 > lplnumk2 ) then

    !-----if already allocated, deallocate arrays
    if( allocated(nga) ) then
        deallocate( nga, ngb, ngc, gx, gy, gz, recnrm,  &
& recnrmex, sumtmp, buff, stat=status )

        the_mem = 4.d0 * lplnumk2 * 3 + 8.d0 * lplnumk2 * 13

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp', .true. )
    end if

    !------allocate memory
    lplnumk2 = numk2 * pwscale
    allocate( nga(lplnumk2), ngb(lplnumk2), ngc(lplnumk2),  &
& gx(lplnumk2), gy(lplnumk2), gz(lplnumk2), recnrm(lplnumk2),  &
& recnrmex(0:lplnumk2),  &
& sumtmp(4*lplnumk2), buff(4*lplnumk2),  &
& stat=status )

    the_mem = 4.d0 * lplnumk2 * 3 + 8.d0 * lplnumk2 * 13

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp', .true. )

end if
!-----allocate memory for reciprocal vectors--------------------------


do i = 1, numk2
   j = ix(i)
   nga(i) = ktmp1(j)
   ngb(i) = ktmp2(j)
   ngc(i) = ktmp3(j)
   d1h = nga(i)
   d2h = ngb(i)
   d3h = ngc(i)
   gx(i) = cross(1,1)*d1h + cross(1,2)*d2h + cross(1,3)*d3h
   gy(i) = cross(2,1)*d1h + cross(2,2)*d2h + cross(2,3)*d3h
   gz(i) = cross(3,1)*d1h + cross(3,2)*d2h + cross(3,3)*d3h
   recnrm(i) = gss(j)
end do
if( .not.ldouble_grid_recip ) then
    do i = 0, numk2
       recnrmex(i) = 0.d0
    end do
  else
    const2 = -(pi/dgalpha)**2
    recnrmex(0) = pi/dgalpha**2
    do i = 1, numk2
       recnrmex(i) = exp(const2*recnrm(i))
    end do
end if


KCSMA1 = 0
KCSMA2 = 0
KCSMA3 = 0
do i = 1, numk2
   KCSMA1 = MAX( KCSMA1, ABS(nga(i)) )
   KCSMA2 = MAX( KCSMA2, ABS(ngb(i)) )
   KCSMA3 = MAX( KCSMA3, ABS(ngc(i)) )
end do
!cc      KCSMAX = MAX( KCSMA1, KCSMA2, KCSMA3 )
if(loutfile(1)) write(nfile(1),*) '   Ewald:   KCSMAX1,2,3',  &
&                                  KCSMA1, KCSMA2, KCSMA3
if(loutfile(2)) write(nfile(2),*) '   Ewald:   KCSMAX1,2,3',  &
&                                  KCSMA1, KCSMA2, KCSMA3
!cc          write(*,*) '   Ewald:   KCSMAX     ', KCSMAX


if( .not.lvshape ) then

    !------deallocate arrays
    deallocate( gss, ktmp1, ktmp2, ktmp3, key,  &
& ix, dstkey, x1, x2, stat=status )
    dealocmem = 8.d0 * ifftk0 * 1.d0 + 4.d0 * ifftk0 * 8.d0
    if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d6),' MB, deallocated (clmrcp)'
    if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d6),' MB, deallocated (clmrcp)'

    ifftk0 = 0

end if


!-----allocate memory for reciprocal vectors----------------------------

nion = max( nion_nod, 1 )

    !-----check array sizes
    if( KCSMA1 > lplKCSMA1 .or.  &
&       KCSMA2 > lplKCSMA2 .or.  &
&       KCSMA3 > lplKCSMA3     ) then

        !-----if already allocated, deallocate arrays
        if( allocated(DCQ1) ) then
            deallocate( DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& stat=status )

            the_mem =  &
&  8.d0 * ( 2*lplKCSMA1 +  2*lplKCSMA2 +  2*lplKCSMA3 + 3 )*nion*2

            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp(2)', .true. )
        end if

        !------allocate memory
        lplKCSMA1 = KCSMA1 * pwscale*pwscale
        lplKCSMA2 = KCSMA2 * pwscale*pwscale
        lplKCSMA3 = KCSMA3 * pwscale*pwscale

        allocate(  &
& DCQ1(-lplKCSMA1:lplKCSMA1,nion), DCQ2(-lplKCSMA2:lplKCSMA2,nion),  &
& DCQ3(-lplKCSMA3:lplKCSMA3,nion), DSQ1(-lplKCSMA1:lplKCSMA1,nion),  &
& DSQ2(-lplKCSMA2:lplKCSMA2,nion), DSQ3(-lplKCSMA3:lplKCSMA3,nion),  &
& stat=status )

        the_mem =  &
&  8.d0 * ( 2*lplKCSMA1 +  2*lplKCSMA2 +  2*lplKCSMA3 + 3 )*nion*2

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp(2)', .true. )

    end if




if( llclpp_r ) then

    mshxx = mshx
    mshyx = mshy
    mshzx = mshz
    call gimax(mshxx)
    call gimax(mshyx)
    call gimax(mshzx)

    !-----check array sizes
    if( KCSMA1 > mplKCSMA1 .or.  &
&       KCSMA2 > mplKCSMA2 .or.  &
&       KCSMA3 > mplKCSMA3 .or.  &
&       mshxx  > mplmshx   .or.  &
&       mshyx  > mplmshy   .or.  &
&       mshzx  > mplmshz       ) then

        !-----if already allocated, deallocate arrays
        if( allocated(DCm1) ) then
            deallocate( DCm1, DCm2, DCm3, DSm1, DSm2, DSm3,  &
& stat=status )

            the_mem = 8.d0 * ( ( 2*mplKCSMA1 + 1 )*mplmshx  &
&                            + ( 2*mplKCSMA2 + 1 )*mplmshy  &
&                            + ( 2*mplKCSMA3 + 1 )*mplmshz )*2

            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp(3)', .true. )
        end if

        !------allocate memory
        mplKCSMA1 = KCSMA1 * pwscale
        mplKCSMA2 = KCSMA2 * pwscale
        mplKCSMA3 = KCSMA3 * pwscale
        mplmshx   = mshxx  * pwscale
        mplmshy   = mshyx  * pwscale
        mplmshz   = mshzx  * pwscale

        allocate( DCm1(-mplKCSMA1:mplKCSMA1,mplmshx),  &
&                 DCm2(-mplKCSMA2:mplKCSMA2,mplmshy),  &
&                 DCm3(-mplKCSMA3:mplKCSMA3,mplmshz),  &
&                 DSm1(-mplKCSMA1:mplKCSMA1,mplmshx),  &
&                 DSm2(-mplKCSMA2:mplKCSMA2,mplmshy),  &
&                 DSm3(-mplKCSMA3:mplKCSMA3,mplmshz),  &
& stat=status )

        the_mem = 8.d0 * ( ( 2*mplKCSMA1 + 1 )*mplmshx  &
&                        + ( 2*mplKCSMA2 + 1 )*mplmshy  &
&                        + ( 2*mplKCSMA3 + 1 )*mplmshz )*2

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'clmrcp(3)', .true. )

    end if

end if


RETURN
END




!      subroutine stgexp
!c-----------------------------------------------------------------------
!c     array for (4/volume/pi) * exp((-pi*G/gamma)**2)/G**2
!c-----------------------------------------------------------------------
!      use para_ewald
!      implicit real*8 ( a-h, o-z )
!
!      const  = 4.d0/(volume*pi)
!      const2 = -(pi/gamma)**2
!      do i = 1, numk2
!         recnrm(i) = const*dexp(const2*recnrm(i))/recnrm(i)
!      end do
!
!
!      return
!      end




subroutine seterf(  &
& alp, tberf, tberfa, drcut, tbeff, tbeffa, drcutf, mx1 )
!-----------------------------------------------------------------------
!     set error function table
!
!      table for error function = ( 1 - complementary error function )/r
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
parameter( PAIRUB  = 1.128379167096d0  )
parameter( PAIRUC  = 0.7522527780637d0 )
!  piarub = 2/sqrt(pi)
!  piaruc = 4/(3*sqrt(pi))
real*8, dimension(0:mx1) :: tberf, tberfa
real*8 :: drcut
real*8, dimension(0:mx1) :: tbeff, tbeffa
real*8 :: drcutf

rcut  = alp
rcut2 = rcut*rcut
drcut = rcut2/dble(mx1-2)
drcutf = drcut

xrr = 0.d0
do ir = 1, mx1 - 2
   xrr = xrr + drcut
   r2 = sqrt(xrr)
   tberf(ir) = (1.d0 - DERFNC(r2,1.d-11))/r2
   tbeff(ir) = ( tberf(ir) - PAIRUB*dexp(-xrr) )/xrr
end do
tberf(0)     = PAIRUB
tberf(mx1-1) = 0.d0
tberf(mx1)   = 0.d0
tbeff(0)     = PAIRUC
tbeff(mx1-1) = 0.d0
tbeff(mx1)   = 0.d0
do ir = 0, mx1 - 2
   tberfa(ir) = 2.D0*( tberf(ir) + tberf(ir+2)  &
&                                - 2.d0*tberf(ir+1) )
   tbeffa(ir) = 2.D0*( tbeff(ir) + tbeff(ir+2)  &
&                                - 2.d0*tbeff(ir+1) )
end do

return
end




subroutine revloc(  &
& ntype, nhk1, nhk2, zv, nel, iatoit, iatmpt, natom, nion, llking,  &
& vext, mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
&      mshx1, mshy1, mshz1, mshx, mshy, mshz, mshnod )
!-----------------------------------------------------------------------
!     local pseudopotential in reciprocal space
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )

integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(ntype) :: zv
integer, dimension(natom) :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
logical, dimension(ntype) :: llking

dimension mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension vext(mshnod)


!--- set sine, cosine function for atoms -------------------------------
!      call stscfn

!--- set sine, cosine function for mesh points -------------------------
!      call stscms( mshx1, mshy1, mshz1, mshx, mshy, mshz )


const  = -4.d0/(volume*pi)
const2 = -(pi/gamma)**2
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)

   sumtm1 = 0.d0
   sumtm2 = 0.d0
   do it = 1, ntype
   if( llking(it) ) then

      sum1 = 0.d0
      sum2 = 0.d0
      do k = nhk1(it), nhk2(it)
         DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
         DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)

         DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
         DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
         sum1 = sum1 + DYCK
         sum2 = sum2 + DYSK
      end do
      sumtm1 = sumtm1 + zv(it)*sum1
      sumtm2 = sumtm2 + zv(it)*sum2

   end if
   end do
   sumtmp(2*M-1) = sumtm1
   sumtmp(2*M)   = sumtm2
end do
call gdsum(sumtmp,2*numk2,buff)

#ifndef VECTOR
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)
   dckm = sumtmp(2*M-1)
   dskm = sumtmp(2*M)
!****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )

   hep  = const*dexp(const2*recnrm(m))/recnrm(m)
   hckm = hep*dckm
   hskm = hep*dskm

   do ix = 1, mshx
   do iy = 1, mshy
      DCC12= DCm1(K1M,ix)*DCm2(K2M,iy) - DSm1(K1M,ix)*DSm2(K2M,iy)
      DSC12= DSm1(K1M,ix)*DCm2(K2M,iy) + DCm1(K1M,ix)*DSm2(K2M,iy)
   do iz = 1, mshz
      DYCK = DCC12*DCm3(K3M,iz) - DSC12*DSm3(K3M,iz)
      DYSK = DSC12*DCm3(K3M,iz) + DCC12*DSm3(K3M,iz)
      m1 = mshxyz(ix,iy,iz)
      vext(m1) = vext(m1) + DYCK*hckm + DYSK*hskm
   end do
   end do
   end do
end do
#else
DO M = 1, numk2
   dckm = sumtmp(2*M-1)
   dskm = sumtmp(2*M)
!****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )

   hep  = const*dexp(const2*recnrm(m))/recnrm(m)
   sumtmp(2*M-1) = hep*dckm
   sumtmp(2*M) = hep*dskm
end do

do ix = 1, mshx
do iy = 1, mshy
   DO M = 1, numk2
      K1M = nga(M)
      K2M = ngb(M)
      buff(2*M-1)=  &
&      DCm1(K1M,ix)*DCm2(K2M,iy) - DSm1(K1M,ix)*DSm2(K2M,iy)
      buff(2*M)=  &
&      DSm1(K1M,ix)*DCm2(K2M,iy) + DCm1(K1M,ix)*DSm2(K2M,iy)
   end do
do iz = 1, mshz
   m1 = mshxyz(ix,iy,iz)
   DO M = 1, numk2
      K3M = ngc(M)
      DYCK = buff(2*M-1)*DCm3(K3M,iz) - buff(2*M)*DSm3(K3M,iz)
      DYSK = buff(2*M)*DCm3(K3M,iz) + buff(2*M-1)*DSm3(K3M,iz)
      vext(m1) = vext(m1) + DYCK*sumtmp(2*M-1) + DYSK*sumtmp(2*M)
   end do
end do
end do
end do
#endif

!--- constant term
const3  = 2.d0*pi*dble(nel)/(volume*gamma*gamma)
do m1 = 1, mshnod
   vext( m1 ) = vext( m1 ) + const3
end do

return
end




#ifndef VECTOR
subroutine refloc( nfile, myid, nodes,  &
& floc, strloc, ntype, nhk1, nhk2, zv, nel, nelv, iatoit,  &
& iatmpt, natom, nion, llking, rdelv, lcstress, rho, nmt0,  &
& dyc, dys, mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, mshnod )
!-----------------------------------------------------------------------
!     atomic force by local pseudopotential in reciprocal space
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
real*8,  dimension(3,natom) :: floc
real*8,  dimension(3,3)  :: strloc
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(ntype) :: zv
integer, dimension(natom) :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
logical, dimension(ntype) :: llking
real*8,  dimension(3*natom) :: dyc, dys

dimension      rho(nmt0)
dimension      mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)

dimension bufst(6)
logical lcstress


!--- set sine, cosine function for atoms -------------------------------
!cc      call stscfn

!--- set sine, cosine function for mesh points -------------------------
!cc      call stscms( mshx1, mshy1, mshz1, mshx, mshy, mshz )


do i = 1, 6
   bufst(i) = 0.d0
end do
const  = -8.d0/volume
const2 = -(pi/gamma)**2
strdiag = 0.d0
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)

   sumtmp(1) = 0.d0
   sumtmp(2) = 0.d0
   do it = 1, ntype
   if( llking(it) ) then
   do k = nhk1(it), nhk2(it)
      DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
      DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)

      DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
      DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
      dyck = zv(it)*DYCK
      dysk = zv(it)*DYSK
      sumtmp(1) = sumtmp(1) + dyck
      sumtmp(2) = sumtmp(2) + dysk
      dyc(k) = dyck
      dys(k) = dysk
   end do
   end if
   end do
!cc         call gdsum(sumtmp,2,buff)
!****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )

   sumtmp(3) = 0.d0
   sumtmp(4) = 0.d0
   do ix = 1, mshx
   do iy = 1, mshy
      DCC12= DCm1(K1M,ix)*DCm2(K2M,iy) - DSm1(K1M,ix)*DSm2(K2M,iy)
      DSC12= DSm1(K1M,ix)*DCm2(K2M,iy) + DCm1(K1M,ix)*DSm2(K2M,iy)
   do iz = 1, mshz
      DYCK = DCC12*DCm3(K3M,iz) - DSC12*DSm3(K3M,iz)
      DYSK = DSC12*DCm3(K3M,iz) + DCC12*DSm3(K3M,iz)
      m1 = mshxyz(ix,iy,iz)
      sumtmp(3) = sumtmp(3) + rho(m1)*DYCK
      sumtmp(4) = sumtmp(4) + rho(m1)*DYSK
   end do
   end do
   end do

   call gdsum(sumtmp,4,buff)
   dckm  = sumtmp(1)
   dskm  = sumtmp(2)
   dckmm = sumtmp(3)
   dskmm = sumtmp(4)

   recr    = 1.d0/recnrm(m)
   hepexp  = dexp(const2*recnrm(m))*recr

   hep  = const*hepexp
   hckm = hep*dckmm
   hskm = hep*dskmm
   dckmx = hckm*gx(m)
   dckmy = hckm*gy(m)
   dckmz = hckm*gz(m)
   dskmx = hskm*gx(m)
   dskmy = hskm*gy(m)
   dskmz = hskm*gz(m)

   do it = 1, ntype
   if( llking(it) ) then
   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)
      floc(1,ia) = floc(1,ia) + dckmx*dys(k) - dskmx*dyc(k)
      floc(2,ia) = floc(2,ia) + dckmy*dys(k) - dskmy*dyc(k)
      floc(3,ia) = floc(3,ia) + dckmz*dys(k) - dskmz*dyc(k)
   end do
   end if
   end do

!--- internal stress tensor
   if( lcstress ) then
       hep    = (recr - const2)*hep
       reeze  = dckm*dckmm + dskm*dskmm
       reezeh = hepexp*reeze
       hckm   = hep*reeze

       dckmx = hckm*gx(m)
       dckmy = hckm*gy(m)
       dckmz = hckm*gz(m)
       strdiag = strdiag + reezeh
       bufst(1) = bufst(1) + dckmx*gx(m)
       bufst(2) = bufst(2) + dckmy*gy(m)
       bufst(3) = bufst(3) + dckmz*gz(m)
       bufst(4) = bufst(4) + dckmy*gz(m)
       bufst(5) = bufst(5) + dckmz*gx(m)
       bufst(6) = bufst(6) + dckmx*gy(m)
   end if

end do


!--- constant term
if( lcstress ) then
    const3  = 2.d0*pi*dble(nel*nelv)/(volume*gamma*gamma)/rdelv
    consth =  4.d0/(volume*pi)
    bufst(1) = bufst(1)/pi + consth*strdiag - const3
    bufst(2) = bufst(2)/pi + consth*strdiag - const3
    bufst(3) = bufst(3)/pi + consth*strdiag - const3

    strloc(1,1) = strloc(1,1) + bufst(1) * rdelv
    strloc(2,2) = strloc(2,2) + bufst(2) * rdelv
    strloc(3,3) = strloc(3,3) + bufst(3) * rdelv
    strloc(2,3) = strloc(2,3) + bufst(4) * rdelv
    strloc(3,1) = strloc(3,1) + bufst(5) * rdelv
    strloc(1,2) = strloc(1,2) + bufst(6) * rdelv
    strloc(3,2) = strloc(2,3)
    strloc(1,3) = strloc(3,1)
    strloc(2,1) = strloc(1,2)
end if


return
end
#else
subroutine refloc( nfile, myid, nodes,  &
& floc, strloc, ntype, nhk1, nhk2, zv, nel, nelv, iatoit,  &
& iatmpt, natom, nion, llking, rdelv, lcstress, rho, nmt0,  &
& dyc, dys, mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, mshnod )
!-----------------------------------------------------------------------
!     atomic force by local pseudopotential in reciprocal space
!
!      for vector processors
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
real*8,  dimension(3,natom) :: floc
real*8,  dimension(3,3)  :: strloc
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(ntype) :: zv
integer, dimension(natom)     :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
logical, dimension(ntype) :: llking
real*8,  dimension(3*natom) :: dyc, dys

dimension      rho(nmt0)
dimension      mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)

dimension bufst(6)
logical lcstress


!--- set sine, cosine function for atoms -------------------------------
!cc      call stscfn

!--- set sine, cosine function for mesh points -------------------------
!cc      call stscms( mshx1, mshy1, mshz1, mshx, mshy, mshz )

!      ct  = timecnt()
!      ct0 = ct

do i = 1, 6
   bufst(i) = 0.d0
end do
const  = -8.d0/volume
const2 = -(pi/gamma)**2
strdiag = 0.d0
!      DO M = 1, numk2
!         K1M = nga(M)
!         K2M = ngb(M)
!         K3M = ngc(M)
!
!         sumtm1 = 0.d0
!         sumtm2 = 0.d0
!         do k = 1, nion
!            ia = iatmpt(k)
!            it = iatoit(ia)
!            DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
!            DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)
!
!            DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
!            DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
!            dyck = zv(it)*DYCK
!            dysk = zv(it)*DYSK
!            sumtm1 = sumtm1 + dyck
!            sumtm2 = sumtm2 + dysk
!         end do
!         sumtmp(M)       = sumtm1
!         sumtmp(M+numk2) = sumtm2
!ccc         call gdsum(sumtmp,2,buff)
!C****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!C****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )
!      end do

DO M = 1, numk2
   sumtmp(M)       = 0.d0
   sumtmp(M+numk2) = 0.d0
end do
do it = 1, ntype
if( llking(it) ) then
do k = nhk1(it), nhk2(it)
   DO M = 1, numk2
      K1M = nga(M)
      K2M = ngb(M)
      buff(2*M-1)  &
&     = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
      buff(2*M)  &
&     = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)
   end do

   DO M = 1, numk2
      K3M = ngc(M)
      DYCK = buff(2*M-1)*DCQ3(K3M,K) - buff(2*M)*DSQ3(K3M,K)
      DYSK = buff(2*M)*DCQ3(K3M,K) + buff(2*M-1)*DSQ3(K3M,K)
      dyck = zv(it)*DYCK
      dysk = zv(it)*DYSK
      sumtmp(M)       = sumtmp(M)       + dyck
      sumtmp(M+numk2) = sumtmp(M+numk2) + dysk
   end do
end do
end if
end do
!cc         call gdsum(sumtmp,2,buff)
!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-1:', ct - ct0
!      ct0 = ct

!--- This part is the most time-consuming. -----------------------------
! Benchmark for 98-atom GeSe2 system on 8-node VPP
!    refloc-1 : 3.6e-03 [sec]
!    refloc-2 : 5.7
!    refloc-3 : 1.3e-02
!    refloc-4 : 3.3e-04
!    refloc-5 : 5.6e-03
!    refloc-6 : 2.8e-04

DO M = 1, numk2
   sumtmp(M+2*numk2) = 0.d0
   sumtmp(M+3*numk2) = 0.d0
end do
do ix = 1, mshx
do iy = 1, mshy
   DO M = 1, numk2
      K1M = nga(M)
      K2M = ngb(M)
      buff(2*M-1)  &
&     = DCm1(K1M,ix)*DCm2(K2M,iy) - DSm1(K1M,ix)*DSm2(K2M,iy)
      buff(2*M)  &
&     = DSm1(K1M,ix)*DCm2(K2M,iy) + DCm1(K1M,ix)*DSm2(K2M,iy)
   end do
do iz = 1, mshz
   m1 = mshxyz(ix,iy,iz)
   DO M = 1, numk2
      K3M = ngc(M)
      DYCK = buff(2*M-1)*DCm3(K3M,iz) - buff(2*M)*DSm3(K3M,iz)
      DYSK = buff(2*M)*DCm3(K3M,iz) + buff(2*M-1)*DSm3(K3M,iz)
      sumtmp(M+2*numk2) = sumtmp(M+2*numk2) + rho(m1)*DYCK
      sumtmp(M+3*numk2) = sumtmp(M+3*numk2) + rho(m1)*DYSK
   end do
end do
end do
end do
!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-2:', ct - ct0
!      ct0 = ct
!-----------------------------------------------------------------------

call gdsum(sumtmp,4*numk2,buff)

!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-3:', ct - ct0
!      ct0 = ct


DO M = 1, numk2
   recr    = 1.d0/recnrm(m)
   hepexp  = dexp(const2*recnrm(m))*recr
   buff(2*m-1) = recr
   buff(2*m)   = hepexp
end do

!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-4:', ct - ct0
!      ct0 = ct

do k = 1, nion
   dys(3*k-2) = 0.d0
   dys(3*k-1) = 0.d0
   dys(3*k)   = 0.d0
end do
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)

   do it = 1, ntype
   if( llking(it) ) then
   do k = nhk1(it), nhk2(it)
      DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
      DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)

      DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
      DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
      dyck = zv(it)*DYCK
      dysk = zv(it)*DYSK
      dyc(2*k-1) = dyck
      dyc(2*k)   = dysk
   end do
   end if
   end do
   dckm  = sumtmp(M)
   dskm  = sumtmp(M+  numk2)
   dckmm = sumtmp(M+2*numk2)
   dskmm = sumtmp(M+3*numk2)

   recr   = buff(2*m-1)
   hepexp = buff(2*m)
!         recr    = 1.d0/recnrm(m)
!         hepexp  = dexp(const2*recnrm(m))*recr

   hep  = const*hepexp
   hckm = hep*dckmm
   hskm = hep*dskmm
   dckmx = hckm*gx(m)
   dckmy = hckm*gy(m)
   dckmz = hckm*gz(m)
   dskmx = hskm*gx(m)
   dskmy = hskm*gy(m)
   dskmz = hskm*gz(m)

   do it = 1, ntype
   if( llking(it) ) then
   do k = nhk1(it), nhk2(it)
      dys(3*k-2) = dys(3*k-2) + dckmx*dyc(2*k) - dskmx*dyc(2*k-1)
      dys(3*k-1) = dys(3*k-1) + dckmy*dyc(2*k) - dskmy*dyc(2*k-1)
      dys(3*k)   = dys(3*k)   + dckmz*dyc(2*k) - dskmz*dyc(2*k-1)
   end do
   end if
   end do
end do
do it = 1, ntype
if( llking(it) ) then
do k = nhk1(it), nhk2(it)
   ia = iatmpt(k)
   floc(1,ia) = floc(1,ia) + dys(3*k-2)
   floc(2,ia) = floc(2,ia) + dys(3*k-1)
   floc(3,ia) = floc(3,ia) + dys(3*k)
end do
end if
end do

!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-5:', ct - ct0
!      ct0 = ct

if( .not.lcstress ) return
!--- internal stress tensor
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)
   dckm  = sumtmp(M)
   dskm  = sumtmp(M+  numk2)
   dckmm = sumtmp(M+2*numk2)
   dskmm = sumtmp(M+3*numk2)

   recr   = buff(2*m-1)
   hepexp = buff(2*m)
!         recr    = 1.d0/recnrm(m)
!         hepexp  = dexp(const2*recnrm(m))*recr

   hep  = const*hepexp
   hckm = hep*dckmm
   hskm = hep*dskmm

   hep    = (recr - const2)*hep
   reeze  = dckm*dckmm + dskm*dskmm
   reezeh = hepexp*reeze
   hckm   = hep*reeze

   dckmx = hckm*gx(m)
   dckmy = hckm*gy(m)
   dckmz = hckm*gz(m)
   strdiag = strdiag + reezeh
   bufst(1) = bufst(1) + dckmx*gx(m)
   bufst(2) = bufst(2) + dckmy*gy(m)
   bufst(3) = bufst(3) + dckmz*gz(m)
   bufst(4) = bufst(4) + dckmy*gz(m)
   bufst(5) = bufst(5) + dckmz*gx(m)
   bufst(6) = bufst(6) + dckmx*gy(m)

end do

!      ct  = timecnt()
!      if( myid.eq.0 ) write(*,*) 'refloc-6:', ct - ct0
!      ct0 = ct

!--- constant term
const3  = 2.d0*pi*dble(nel*nelv)/(volume*gamma*gamma)/rdelv
consth =  4.d0/(volume*pi)
bufst(1) = bufst(1)/pi + consth*strdiag - const3
bufst(2) = bufst(2)/pi + consth*strdiag - const3
bufst(3) = bufst(3)/pi + consth*strdiag - const3

strloc(1,1) = strloc(1,1) + bufst(1) * rdelv
strloc(2,2) = strloc(2,2) + bufst(2) * rdelv
strloc(3,3) = strloc(3,3) + bufst(3) * rdelv
strloc(2,3) = strloc(2,3) + bufst(4) * rdelv
strloc(3,1) = strloc(3,1) + bufst(5) * rdelv
strloc(1,2) = strloc(1,2) + bufst(6) * rdelv
strloc(3,2) = strloc(2,3)
strloc(1,3) = strloc(3,1)
strloc(2,1) = strloc(1,2)


return
end
#endif




#ifndef VECTOR
subroutine frcrec( fclm, strclm, ntype, nhk1, nhk2, zv, nel,  &
& ecolmb, atom_ecolmb, iatoit, iatmpt, natom, nion,  &
& dyc, dys, lcstress )
!-----------------------------------------------------------------------
!     direct coulomb interaction in reciprocal space
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )

real*8,  dimension(3,natom) :: fclm
real*8,  dimension(3,3)  :: strclm
integer :: ntype
real*8  :: ecolmb
real*8,  dimension(natom) :: atom_ecolmb
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(natom) :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
real*8,  dimension(3*natom) :: dyc, dys

dimension bufst(6)
logical   lcstress


!--- set sine, cosine function for atoms -------------------------------
!cc      call stscfn


do i = 1, 6
   bufst(i) = 0.d0
end do

vol8   = 8.d0/volume
const2 = -(pi/gamma)**2
const  = 2.d0/(volume*pi)

recolmb = 0.d0
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)

   sumtmp(1) = 0.d0
   sumtmp(2) = 0.d0
   do k = 1, nion
      ia = iatmpt(k)
      it = iatoit(ia)
      DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
      DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)

      DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
      DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
      dyck = zv(it)*dyck
      dysk = zv(it)*dysk
      sumtmp(1) = sumtmp(1) + DYCK
      sumtmp(2) = sumtmp(2) + DYSK
      dyc(k) = dyck
      dys(k) = dysk
   end do
   call gdsum(sumtmp,2,buff)
   dckm = sumtmp(1)
   dskm = sumtmp(2)
!****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )

   recr = 1.d0/recnrm(m)
   hep  = ( dexp(const2*recnrm(m)) - recnrmex(m) )*recr
   hepdcs  = hep*( dckm*dckm + dskm*dskm )
   recolmb = recolmb + hepdcs

   hep8 = hep*vol8
   hckm = hep8*dckm
   hskm = hep8*dskm
   dckmx = hckm*gx(m)
   dckmy = hckm*gy(m)
   dckmz = hckm*gz(m)
   dskmx = hskm*gx(m)
   dskmy = hskm*gy(m)
   dskmz = hskm*gz(m)
   hepc = hep*const
   hckm = hepc*dckm
   hskm = hepc*dskm
   do k = 1, nion
      ia = iatmpt(k)
      fclm(1,ia) = fclm(1,ia) + dckmx*dys(k) - dskmx*dyc(k)
      fclm(2,ia) = fclm(2,ia) + dckmy*dys(k) - dskmy*dyc(k)
      fclm(3,ia) = fclm(3,ia) + dckmz*dys(k) - dskmz*dyc(k)
      it = iatoit(ia)
      atom_ecolmb(ia) = atom_ecolmb(ia)  &
&                     + dyc(k)*hckm + dys(k)*hskm
   end do

   if( lcstress ) then
       hckm  = (recr - const2)*hepdcs
       dckmx = hckm*gx(m)
       dckmy = hckm*gy(m)
       dckmz = hckm*gz(m)
       bufst(1) = bufst(1) + dckmx*gx(m)
       bufst(2) = bufst(2) + dckmy*gy(m)
       bufst(3) = bufst(3) + dckmz*gz(m)
       bufst(4) = bufst(4) + dckmy*gz(m)
       bufst(5) = bufst(5) + dckmz*gx(m)
       bufst(6) = bufst(6) + dckmx*gy(m)
   end if

end do
const  = 2.d0/(volume*pi)
constb = 4.d0/(volume*pi)
ecolmb = ecolmb + const*recolmb
if( lcstress ) then
    bufst(1) = bufst(1) * constb
    bufst(2) = bufst(2) * constb
    bufst(3) = bufst(3) * constb
    bufst(4) = bufst(4) * constb
    bufst(5) = bufst(5) * constb
    bufst(6) = bufst(6) * constb
end if


!--- constant term
ztot2 = 0.d0
do it = 1, ntype
   ztot2 = ztot2 + dble(nhk2(it) - nhk1(it) + 1)*zv(it)*zv(it)
end do
ecolmb = ecolmb  &
&      - dble(nel*nel)/volume* ( pi/(gamma*gamma) - recnrmex(0) )  &
&      - 2.d0*gamma*ztot2/sqrt(pi)
do k = 1, nion
   ia = iatmpt(k)
   it = iatoit(ia)
   atom_ecolmb(ia) = atom_ecolmb(ia)  &
&    - zv(it)*dble(nel)/volume* ( pi/(gamma*gamma) - recnrmex(0) )  &
&    - 2.d0*gamma*zv(it)*zv(it)/sqrt(pi)
end do

if( lcstress ) then
    constr = - const*recolmb  &
&            + pi*dble(nel*nel)/(volume*gamma*gamma)
    bufst(1) = bufst(1) + constr
    bufst(2) = bufst(2) + constr
    bufst(3) = bufst(3) + constr


    strclm(1,1) = strclm(1,1) + bufst(1)
    strclm(2,2) = strclm(2,2) + bufst(2)
    strclm(3,3) = strclm(3,3) + bufst(3)
    strclm(2,3) = strclm(2,3) + bufst(4)
    strclm(3,1) = strclm(3,1) + bufst(5)
    strclm(1,2) = strclm(1,2) + bufst(6)
    strclm(3,2) = strclm(2,3)
    strclm(1,3) = strclm(3,1)
    strclm(2,1) = strclm(1,2)
end if


return
end
#else
subroutine frcrec( fclm, strclm, ntype, nhk1, nhk2, zv, nel,  &
& ecolmb, atom_ecolmb, iatoit, iatmpt, natom, nion,  &
& dyc, dys, lcstress )
!-----------------------------------------------------------------------
!     direct coulomb interaction in reciprocal space
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )

real*8,  dimension(3,natom) :: fclm
real*8,  dimension(3,3)  :: strclm
integer :: ntype
real*8  :: ecolmb
real*8,  dimension(natom) :: atom_ecolmb
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(natom) :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
real*8,  dimension(3*natom) :: dyc, dys

dimension bufst(6)
logical   lcstress


!--- set sine, cosine function for atoms -------------------------------
!cc      call stscfn


do i = 1, 6
   bufst(i) = 0.d0
end do

vol8   = 8.d0/volume
const2 = -(pi/gamma)**2
const  = 2.d0/(volume*pi)


DO M = 1, 2*numk2
   sumtmp(m) = 0.d0
end do
do k = 1, nion
   ia = iatmpt(k)
   it = iatoit(ia)
   DO M = 1, numk2
      K1M = nga(M)
      K2M = ngb(M)
      buff(2*m-1) =DCQ1(K1M,K)*DCQ2(K2M,K)-DSQ1(K1M,K)*DSQ2(K2M,K)
      buff(2*m  ) =DSQ1(K1M,K)*DCQ2(K2M,K)+DCQ1(K1M,K)*DSQ2(K2M,K)
   end do
   DO M = 1, numk2
      K3M = ngc(M)
      DYCK = buff(2*m-1)*DCQ3(K3M,K) - buff(2*m  )*DSQ3(K3M,K)
      DYSK = buff(2*m  )*DCQ3(K3M,K) + buff(2*m-1)*DSQ3(K3M,K)
      dyck = zv(it)*dyck
      dysk = zv(it)*dysk
      sumtmp(m      ) = sumtmp(m      ) + DYCK
      sumtmp(m+numk2) = sumtmp(m+numk2) + DYSK
   end do
end do
call gdsum(sumtmp,2*numk2,buff)


DO M = 1, numk2
   recr = 1.d0/recnrm(m)
   hep  = ( dexp(const2*recnrm(m)) - recnrmex(m) )*recr
   sumtmp(m+2*numk2) = recr
   sumtmp(m+3*numk2) = hep
end do

recolmb = 0.d0
do k = 1, nion
   dys(3*k-2) = 0.d0
   dys(3*k-1) = 0.d0
   dys(3*k  ) = 0.d0
   dyc(3*k-2) = 0.d0
end do
DO M = 1, numk2
   K1M = nga(M)
   K2M = ngb(M)
   K3M = ngc(M)

   do k = 1, nion
      ia = iatmpt(k)
      it = iatoit(ia)
      DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) - DSQ1(K1M,K)*DSQ2(K2M,K)
      DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) + DCQ1(K1M,K)*DSQ2(K2M,K)

      DYCK = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
      DYSK = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
      dyck = zv(it)*dyck
      dysk = zv(it)*dysk
      dyc(3*k-1) = dyck
      dyc(3*k  ) = dysk
   end do
   dckm = sumtmp(m)
   dskm = sumtmp(m+numk2)
!****        DCKM -> C(KM) = SUM( ZV*COS(2*PAI*KM*R) )
!****        DSKM -> S(KM) = SUM( ZV*SIN(2*PAI*KM*R) )

   hep  = sumtmp(m+3*numk2)
   hepdcs  = hep*( dckm*dckm + dskm*dskm )
   recolmb = recolmb + hepdcs

   hep8 = hep*vol8
   hckm = hep8*dckm
   hskm = hep8*dskm
   dckmx = hckm*gx(m)
   dckmy = hckm*gy(m)
   dckmz = hckm*gz(m)
   dskmx = hskm*gx(m)
   dskmy = hskm*gy(m)
   dskmz = hskm*gz(m)
   hepc = hep*const
   hckm = hepc*dckm
   hskm = hepc*dskm
   do k = 1, nion
      dys(3*k-2) = dys(3*k-2) + dckmx*dyc(3*k) - dskmx*dyc(3*k-1)
      dys(3*k-1) = dys(3*k-1) + dckmy*dyc(3*k) - dskmy*dyc(3*k-1)
      dys(3*k  ) = dys(3*k  ) + dckmz*dyc(3*k) - dskmz*dyc(3*k-1)
      dyc(3*k-2) = dyc(3*k-2) + dyc(3*k-1)*hckm + dyc(3*k)*hskm
   end do

end do
do k = 1, nion
   ia = iatmpt(k)
   fclm(1,ia) = fclm(1,ia) + dys(3*k-2)
   fclm(2,ia) = fclm(2,ia) + dys(3*k-1)
   fclm(3,ia) = fclm(3,ia) + dys(3*k)
   it = iatoit(ia)
   atom_ecolmb(ia) = atom_ecolmb(ia) + dyc(3*k-2)
end do
const  = 2.d0/(volume*pi)
ecolmb = ecolmb + const*recolmb

!--- constant term
ztot2 = 0.d0
do it = 1, ntype
   ztot2 = ztot2 + dble(nhk2(it) - nhk1(it) + 1)*zv(it)*zv(it)
end do
ecolmb = ecolmb  &
&      - dble(nel*nel)/volume* ( pi/(gamma*gamma) - recnrmex(0) )  &
&      - 2.d0*gamma*ztot2/sqrt(pi)
do k = 1, nion
   ia = iatmpt(k)
   it = iatoit(ia)
   atom_ecolmb(ia) = atom_ecolmb(ia)  &
&    - zv(it)*dble(nel)/volume* ( pi/(gamma*gamma) - recnrmex(0) )  &
&    - 2.d0*gamma*zv(it)*zv(it)/sqrt(pi)
end do


if( lcstress ) then
    DO M = 1, numk2
       dckm = sumtmp(m)
       dskm = sumtmp(m+numk2)
       recr = sumtmp(m+2*numk2)
       hep  = sumtmp(m+3*numk2)
       hepdcs  = hep*( dckm*dckm + dskm*dskm )
       hckm  = (recr - const2)*hepdcs
       dckmx = hckm*gx(m)
       dckmy = hckm*gy(m)
       dckmz = hckm*gz(m)
       bufst(1) = bufst(1) + dckmx*gx(m)
       bufst(2) = bufst(2) + dckmy*gy(m)
       bufst(3) = bufst(3) + dckmz*gz(m)
       bufst(4) = bufst(4) + dckmy*gz(m)
       bufst(5) = bufst(5) + dckmz*gx(m)
       bufst(6) = bufst(6) + dckmx*gy(m)
    end do

    constb = 4.d0/(volume*pi)
    bufst(1) = bufst(1) * constb
    bufst(2) = bufst(2) * constb
    bufst(3) = bufst(3) * constb
    bufst(4) = bufst(4) * constb
    bufst(5) = bufst(5) * constb
    bufst(6) = bufst(6) * constb

!--- constant term
    constr = - const*recolmb  &
&            + pi*dble(nel*nel)/(volume*gamma*gamma)
    bufst(1) = bufst(1) + constr
    bufst(2) = bufst(2) + constr
    bufst(3) = bufst(3) + constr


    strclm(1,1) = strclm(1,1) + bufst(1)
    strclm(2,2) = strclm(2,2) + bufst(2)
    strclm(3,3) = strclm(3,3) + bufst(3)
    strclm(2,3) = strclm(2,3) + bufst(4)
    strclm(3,1) = strclm(3,1) + bufst(5)
    strclm(1,2) = strclm(1,2) + bufst(6)
    strclm(3,2) = strclm(2,3)
    strclm(1,3) = strclm(3,1)
    strclm(2,1) = strclm(1,2)
end if


return
end
#endif




subroutine stscfn(  &
& ntype, satm, iatoit, iatmpt, natom, nion,  &
& ldouble_grid_recip, nd1v, nd1vks )
!-----------------------------------------------------------------------
!    set sine, cosine function for atoms
!-----------------------------------------------------------------------
use para_ewald
implicit real*8 ( a-h, o-z )

integer :: ntype
real*8,  dimension(3,natom) :: satm
integer, dimension(natom)   :: iatoit
integer, dimension(nion) :: iatmpt
integer :: nion
logical :: ldouble_grid_recip
integer, dimension(3) :: nd1v, nd1vks


dpi = 2.d0*pi
if( .not.ldouble_grid_recip ) then
    fact1 = dpi
    fact2 = dpi
    fact3 = dpi
  else
    fact1 = dpi*dble(nd1v(1))/dble(nd1vks(1))
    fact2 = dpi*dble(nd1v(2))/dble(nd1vks(2))
    fact3 = dpi*dble(nd1v(3))/dble(nd1vks(3))
end if
do k = 1, nion
   ia = iatmpt(k)
   it = iatoit(ia)

!****  for p(1)
   D2P1  = satm(1,ia)*fact1
   DCSN1 = COS(D2P1)
   DSGN1 = SIN(D2P1)

   DCQ1(0,K)  = 1.0D0
   DSQ1(0,K)  = 0.d0
   DCQ1(1,K)  = DCSN1
   DSQ1(1,K)  = DSGN1
   DCQ1(-1,K) = DCSN1 
   DSQ1(-1,K) = -DSGN1

!****  for p(2)
   D2P2  = satm(2,ia)*fact2
   DCSN2 = COS(D2P2)
   DSGN2 = SIN(D2P2)

   DCQ2(0,K)  = 1.0D0
   DSQ2(0,K)  = 0.d0
   DCQ2(1,K)  = DCSN2
   DSQ2(1,K)  = DSGN2
   DCQ2(-1,K) = DCSN2 
   DSQ2(-1,K) = -DSGN2

!****  for p(3)
   D2P3  = satm(3,ia)*fact3
   DCSN3 = COS(D2P3)
   DSGN3 = SIN(D2P3)

   DCQ3(0,K)  = 1.0D0
   DSQ3(0,K)  = 0.d0
   DCQ3(1,K)  = DCSN3
   DSQ3(1,K)  = DSGN3
   DCQ3(-1,K) = DCSN3 
   DSQ3(-1,K) = -DSGN3

   DO KM = 2, KCSMA1
      KMM = KM - 1
      DCQ1(KM,K)  = DCSN1*DCQ1(KMM,K) - DSGN1*DSQ1(KMM,K)
      DSQ1(KM,K)  = DSGN1*DCQ1(KMM,K) + DCSN1*DSQ1(KMM,K)
      DCQ1(-KM,K) = DCQ1(KM,K)
      DSQ1(-KM,K) = -DSQ1(KM,K)
   end do
   DO KM = 2, KCSMA2
      KMM = KM - 1
      DCQ2(KM,K)  = DCSN2*DCQ2(KMM,K) - DSGN2*DSQ2(KMM,K)
      DSQ2(KM,K)  = DSGN2*DCQ2(KMM,K) + DCSN2*DSQ2(KMM,K)
      DCQ2(-KM,K) = DCQ2(KM,K)
      DSQ2(-KM,K) = -DSQ2(KM,K)
   end do
   DO KM = 2, KCSMA3
      KMM = KM - 1
      DCQ3(KM,K)  = DCSN3*DCQ3(KMM,K) - DSGN3*DSQ3(KMM,K)
      DSQ3(KM,K)  = DSGN3*DCQ3(KMM,K) + DCSN3*DSQ3(KMM,K)
      DCQ3(-KM,K) = DCQ3(KM,K)
      DSQ3(-KM,K) = -DSQ3(KM,K)
   end do
end do


return
end




subroutine stscms( nfile, myid, nodes,  &
& nd1v, mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!    set sine, cosine function for mesh points
!-----------------------------------------------------------------------
use para_ewald
implicit none
integer :: nfile(*), myid, nodes
integer :: nd1v(3)
integer :: mshx1, mshy1, mshz1, mshx, mshy, mshz

!-----declare local variables
integer :: k, ix, iy, iz, km, kmm
real*8  :: dpi, q1, q2, q3, D2P1, DCSN1, DSGN1
real*8  :: the_mem
integer :: status
!integer :: ifsetm = 0
!save ifsetm


dpi = 2.d0*pi

do k = 1, mshx
   ix = mshx1 + k - 1
   q1 = dble(ix-1)/dble(nd1v(1))
   D2P1  = q1*dpi
   DCSN1 = COS(D2P1)
   DSGN1 = SIN(D2P1)

   DCm1(0,K)  = 1.0D0
   DSm1(0,K)  = 0.0
   DCm1(1,K)  = DCSN1
   DSm1(1,K)  = DSGN1
   DCm1(-1,K) = DCSN1 
   DSm1(-1,K) = -DSGN1
   DO KM = 2, KCSMA1
      KMM = KM - 1
      DCm1(KM,K)  = DCSN1*DCm1(KMM,K) - DSGN1*DSm1(KMM,K)
      DSm1(KM,K)  = DSGN1*DCm1(KMM,K) + DCSN1*DSm1(KMM,K)
      DCm1(-KM,K) = DCm1(KM,K)
      DSm1(-KM,K) = -DSm1(KM,K)
   end do
end do

do k = 1, mshy
   iy = mshy1 + k - 1
   q2 = dble(iy-1)/dble(nd1v(2))
   D2P1  = q2*dpi
   DCSN1 = COS(D2P1)
   DSGN1 = SIN(D2P1)

   DCm2(0,K)  = 1.0D0
   DSm2(0,K)  = 0.0
   DCm2(1,K)  = DCSN1
   DSm2(1,K)  = DSGN1
   DCm2(-1,K) = DCSN1 
   DSm2(-1,K) = -DSGN1
   DO KM = 2, KCSMA2
      KMM = KM - 1
      DCm2(KM,K)  = DCSN1*DCm2(KMM,K) - DSGN1*DSm2(KMM,K)
      DSm2(KM,K)  = DSGN1*DCm2(KMM,K) + DCSN1*DSm2(KMM,K)
      DCm2(-KM,K) = DCm2(KM,K)
      DSm2(-KM,K) = -DSm2(KM,K)
   end do
end do

do k = 1, mshz
   iz = mshz1 + k - 1
   q3 = dble(iz-1)/dble(nd1v(3))
   D2P1  = q3*dpi
   DCSN1 = COS(D2P1)
   DSGN1 = SIN(D2P1)

   DCm3(0,K)  = 1.0D0
   DSm3(0,K)  = 0.0
   DCm3(1,K)  = DCSN1
   DSm3(1,K)  = DSGN1
   DCm3(-1,K) = DCSN1 
   DSm3(-1,K) = -DSGN1
   DO KM = 2, KCSMA3
      KMM = KM - 1
      DCm3(KM,K)  = DCSN1*DCm3(KMM,K) - DSGN1*DSm3(KMM,K)
      DSm3(KM,K)  = DSGN1*DCm3(KMM,K) + DCSN1*DSm3(KMM,K)
      DCm3(-KM,K) = DCm3(KM,K)
      DSm3(-KM,K) = -DSm3(KM,K)
   end do
end do


return
end




