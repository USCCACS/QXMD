



module nlppr_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkbr.f, nlkbppr.f & nlvandr.f
!-----------------------------------------------------------------------
implicit none

!-----variables for nonlocal pseudopotential
real*8, allocatable, dimension(:,:,:) :: tabnl, tabnla
real*8, allocatable, dimension(:,:,:) :: tbfnl, tbfnla
real*8, allocatable, dimension(:,:)   :: dltnlc, rmxnlc
real*8, allocatable, dimension(:,:)   :: clalp_ncr

real*8  :: rmxnlx   ! max. cutoff length for NCPP
real*8  :: rmxnlxv  ! max. cutoff length for USPP

save


end module




subroutine nlppr_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, mx1, lkbpp_r, lvand_r )
!-----------------------------------------------------------------------
!     allocate memory for variables for atoms and pseudopotentials
!-----------------------------------------------------------------------
use nlpp_parameters
use nlppr_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: ntype
integer :: mxl
integer :: mx1
logical :: lkbpp_r, lvand_r

!-----declare local variables
integer :: status
real*8  :: the_mem


!-----if already allocated, nothing is needed below
if( allocated(tabnl) ) return


!------allocate memory
allocate( tabnl(0:mx1,lnref,ntype), tabnla(0:mx1,lnref,ntype),  &
& tbfnl(0:mx1,lnref,ntype), tbfnla(0:mx1,lnref,ntype),  &
& dltnlc(lnref,ntype), rmxnlc(lnref,ntype), clalp_ncr(lnref,ntype),  &
& stat=status )

the_mem =  &
&  8.d0 * ( (mx1+1)*lnref*ntype*4 + lnref*ntype*3 )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlppr_variables_alloc', .true. )


if( lkbpp_r ) call nlkbppr_prealloc( nfile, myid, nodes,  &
& alloc_mem, lnref )

!if( lvand_r ) call nlvandr_prealloc( nfile, myid, nodes,  &
!& alloc_mem, lnref, mxl )

return
end subroutine




subroutine nlkbpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutsoft, lkbppi,  &
& lking, rking, gkgmax, gkgexct, ms, aname,  &
& mxl, mx1, vlocli ,xitgrd, rr, nvlcl )
!-----------------------------------------------------------------------
!---  Kleinmann and Bylander with Troullier and Martins  ---
!     non-local pseudopotential elements
!-----------------------------------------------------------------------
use outfile
use psvariables
use nlpp_variables
use nlppr_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
logical :: lrela
integer :: ntype
real*8,  dimension(ntype) :: zatom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lclno
logical, dimension(0:mxl,ntype) :: lchk
logical, dimension(ntype) :: lkbppi
logical, dimension(ntype) :: lking
real*8,  dimension(ntype) :: rking, gkgmax, gkgexct
character(1), dimension(7)   :: ms
character(2), dimension(103) :: aname

integer :: nvlcl
real*8, dimension(0:nvlcl)  :: vlocli ,xitgrd, rr
!-----declare local variables
real*8,  allocatable, dimension(:) :: tabwv
integer, parameter :: noofdt = 200
real*8,  allocatable, dimension(:) :: betaq
real*8 :: gdel
real*8,  allocatable, dimension(:) :: betar
real*8 :: rspdel
real*8,  allocatable, dimension(:,:) :: amatrx
real*8,  allocatable, dimension(:,:) :: vmatrx
integer, allocatable, dimension(:)   :: ivmatrx
integer :: status
real*8  :: alocmem, dealocmem

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 1
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip

character(50) :: fname
character(1)  :: dummy
character(1), dimension(0:9) :: num =  &
&                    (/ '0','1','2','3','4','5','6','7','8','9' /)


!------allocate memory
allocate( tabwv(0:mx1), betaq(0:noofdt), betar(0:noofdt),  &
& amatrx(0:noofdt,0:noofdt), vmatrx(noofdt,noofdt+1),  &
& ivmatrx(noofdt),  &
& stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in nlkbpp' )

alocmem = 8.d0 * ( mx1+1 + (noofdt+1)*2 + (noofdt+1)*(noofdt+1)  &
&                 + noofdt*(noofdt+1) )  &
& + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, allocated (nlkbpp)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, allocated (nlkbpp)'


!      do it=1, ntype
!      do l=0, lmax(it)
!      do I = 1, MESH
!         PSUORG(I,l,it) = PSUORG(I,l,it) - PLOCAL(I,it)
!      end do
!      end do
!      end do

ifpp = 1
if( lrela ) ifpp = ifpp + 3

call allocate_unit_number( jpl )


!-----zero clear
do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
       tabnl(0:mx1,1:lnref,it) = 0.d0
       tbfnl(0:mx1,1:lnref,it) = 0.d0
       rmxnlc(1:lnref,it) = 0.d0
       dltnlc(1:lnref,it) = 0.d0
   end if
end do

!-----------------------------------------------------------------------
LS = 5
N0 = 30
NOVERH = LS*2
NOVER  = NOVERH*2

typedo: do it=1, ntype
typeif: if( lkbppi(it) .and. lking(it) ) then
!=======================================================================

MESH =  meshi(it)
RMAX =  rmaxi(it)
DX   =  dxi(it)
DO I = 1, MESH
   R(I)   =  rdial(I,it)
end do

lmxx = 0
laxx = 0
ldo: do l=0, lmax(it)
if( lchk(l,it) .and. l.ne.lclno(it) ) then
jdo: do j = MJTAB1(l,it), MJTAB2(l,it)

!--------------------------------------------
!        set index
      do mj = 1, 2*l + 1
         lmm = mj - ( l + 1 )
         lmmx = lmxx + mj
!               itolm(lmmx,it) = l*(l+1)/2 + abs(lmm) + 1
         litoll(lmmx,it) = mod(l,2) .eq. 0
!               if(        lmm.ge.0 ) then
!                          itolmf(lmmx,it) = 1
!                 else if( mod(l+lmm,2).eq.0 ) then
!                          itolmf(lmmx,it) = 2
!                 else
!                          itolmf(lmmx,it) = 3
!               endif
         itolm( lmmx,it) = l
         itolmf(lmmx,it) = lmm
         itola( lmmx,it) = laxx + 1
      enddo
!--------------------------------------------
   lmm = laxx + 1

   ix = 0
   do I = MESH, 1, -1
      if( abs(PSUORG(I,j,it)).gt.1.0d-10 ) then
          ix = I + 1
          exit
      end if
   end do
   !--- error trap
   if( ix == 0 ) &
&      call fstop( nfile, myid, nodes, 'something wrong with nonlocal pp.' )
   !----------------
   rmxnlc(lmm,it) = R(ix)
   dltnlc(lmm,it) = R(ix)/dble(mx1-2)
!-----------------------------------------------------------------------
!---  read pseudo-wave-functions of reference state
!---       calculated by scalor relativistic atomic calculation.
   fname = (( ms(l+1)//num(rjtab(j)) )//'.nwp' )
   call openpp( nfile, myid, nodes,  &
& ifpp, jpl, zatom(it), fname, jgga, aname, ierror )

   do
      READ(JPL,'(a1)') DUMMY
      IF( DUMMY /= '#' ) exit
   end do
   backspace JPL
   i = 1
   do
      READ(JPL,*,iostat=istat) RI, PNL(i)
      if( istat /= 0 ) exit
      i = i + 1
      IF( i > MESH ) exit
   end do
   PNL(i:MESH) = 0.d0
   CLOSE(JPL)
!-----------------------------------------------------------------------
!     interpolation of wavefunctions  by  spline
   if( l.eq.0 ) then
       tabwv(0) = PNL(1) / R(1)
     else
       tabwv(0) = 0.d0
   end if
   tabwv(mx1-2) = PNL(ix) / R(ix)
   tabwv(mx1-1) = 0.d0
   tabwv(mx1  ) = 0.d0

LL   = LS
INI  = 1
ILSM = 1
xrr = 0.d0
do ir = 1, mx1 - 3
   xrr = xrr + dltnlc(lmm,it)
   IF( xrr.GT.r(ILSM) ) THEN
       do
          ILS  = MIN( INI + N0 - 1, MESH )
          IF( ILS-INI+1.LT.LL+2 ) LL = ILS - INI - 1
          IF( ILS.EQ.MESH ) THEN
              ILSM = ILS
            ELSE
              ILSM = ILS - NOVERH
          end if
          IF( R(ILSM) >= xrr ) exit
          INI = INI + NOVERH
       end do

       Nm = 0
       DO 670 LP = INI, ILS
          Nm = Nm + 1
          XS(Nm) = R(LP)
          YS(Nm) = PNL(LP) / R(LP)
670        CONTINUE
       CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,  &
&                   KQM, MM, MM1, NDV )
       INI  = ILS - NOVER
   end if

   CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
   CALL SPLINV( LL, Nm, xrr, ppfsv, QS, AS, WS, KQM, MM, NDV, 1 )
   tabwv(ir) = ppsv

end do
!-----------------------------------------------------------------------
!     interpolation of pseudopotentials  by  spline
   tabnl(0,lmm,it)     = PSUORG(1,j,it)
   tabnl(mx1-2,lmm,it) = PSUORG(ix,j,it)
   tabnl(mx1-1,lmm,it) = 0.d0
   tabnl(mx1  ,lmm,it) = 0.d0

LL   = LS
INI  = 1
ILSM = 1
xrr = 0.d0
do ir = 1, mx1 - 3
   xrr = xrr + dltnlc(lmm,it)
   IF( xrr.GT.r(ILSM) ) THEN
       do
          ILS  = MIN( INI + N0 - 1, MESH )
          IF( ILS-INI+1.LT.LL+2 ) LL = ILS - INI - 1
          IF( ILS.EQ.MESH ) THEN
              ILSM = ILS
            ELSE
              ILSM = ILS - NOVERH
          end if
          IF( R(ILSM) >= xrr ) exit
          INI = INI + NOVERH
       end do

       Nm = 0
       do LP = INI, ILS
          Nm = Nm + 1
          XS(Nm) = R(LP)
          YS(Nm) = PSUORG(LP,j,it)
       end do
       CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,  &
&                   KQM, MM, MM1, NDV )
       INI  = ILS - NOVER
   END IF

   CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
   CALL SPLINV( LL, Nm, xrr, ppfsv, QS, AS, WS, KQM, MM, NDV, 1 )
   tabnl(ir,lmm,it) = ppsv
   tbfnl(ir,lmm,it) = ppfsv

end do
tbfnl(0,lmm,it)     = 0.d0
tbfnl(mx1-2,lmm,it) = tbfnl(mx1-3,lmm,it)
tbfnl(mx1-1,lmm,it) = 0.d0
tbfnl(mx1,lmm,it)   = 0.d0
!-----------------------------------------------------------------------
clalpt = 0.d0
xrr = 0.d0
do ir = 1, mx1
   xrr = xrr + dltnlc(lmm,it)
   pnlt = tabwv(ir)*xrr
   clalpt = clalpt + tabnl(ir,lmm,it)*pnlt*pnlt
end do
clalp_ncr(lmm,it) = 1.d0/( clalpt*dltnlc(lmm,it) ) * factl(j,it)
!-----------------------------------------------------------------------
do ir = 0, mx1
   tabnl(ir,lmm,it) = tabnl(ir,lmm,it)*tabwv(ir)
end do
!      if( lking(it) ) then
!     --- King-Smith implementation ---
       lltk = l
       iitk = it
       xrr = 0.d0
       do ir = 0, mx1
          rr(ir) = xrr
          vlocli(ir) = tabnl(ir,lmm,it)*xrr*xrr
          xrr = xrr + dltnlc(lmm,it)
       end do
       dkmax4 = sqrt(ecutsoft) * gkgmax(it)
       dkmax  = sqrt(ecutsoft) * gkgexct(it)
       eg4   = dkmax4*dkmax4
       egmax = dkmax*dkmax
       rmxnlc(lmm,it) = rmxnlc(lmm,it) * rking(it)
!            -----------------------------------------------------------
       call kingsm( nfile, myid, nodes,  &
& lltk, iitk, dltnlc(lmm,it), egmax, eg4, .true., rmxnlc(lmm,it),  &
& vlocli ,xitgrd, rr, nvlcl, betaq, gdel, betar, rspdel, noofdt,  &
& amatrx, vmatrx, ivmatrx )
!            -----------------------------------------------------------
       dltnlc(lmm,it) = rmxnlc(lmm,it)/dble(mx1-2)
!        --- interpolation by  spline ---
       rskip = rmxnlc(lmm,it)/dble(mm-1)
       n  = 1
       nr = 0
       xs(n) = rspdel*dble(nr)
       ys(n) = betar(nr)
       xrr = rskip
       do nr = 1, noofdt
          xnrr = rspdel*dble(nr)
          if( xnrr.ge.xrr ) then
              n = n + 1
              xs(n) = xnrr
              ys(n) = betar(nr)
              xrr = xrr + rskip
          end if
       end do
       nr = noofdt
       xnrr = rspdel*dble(nr)
       xs(mm) = xnrr
       ys(mm) = betar(nr)

       nm = mm
       LL = LS
       CALL SPLINS( LL, nm, XS, YS, QS, BS, AS, WS, IP,  &
&                   KQM, MM, MM1, NDV )
       xrr = 0.d0
       do ir = 1, mx1 - 2
          xrr = xrr + dltnlc(lmm,it)
          CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS,  &
&                      KQM, MM, NDV, 0 )
          CALL SPLINV( LL, Nm, xrr, ppfsv, QS, AS, WS,  &
&                      KQM, MM, NDV, 1 )
          tabnl(ir,lmm,it) = ppsv
          tbfnl(ir,lmm,it) = ppfsv
       end do
       tabnl(0,lmm,it)     = betar(0)
       tabnl(mx1-2,lmm,it) = betar(noofdt)
       tabnl(mx1-1,lmm,it) = 0.d0
       tabnl(mx1  ,lmm,it) = 0.d0
       if( l.eq.1 ) then
           tbfnl(0,lmm,it) = tbfnl(1,lmm,it)
         else
           tbfnl(0,lmm,it) = 0.d0
       end if
       tbfnl(mx1-2,lmm,it) = tbfnl(mx1-3,lmm,it)
       tbfnl(mx1-1,lmm,it) = 0.d0
       tbfnl(mx1,lmm,it)   = 0.d0


do ir = 0, mx1 - 2
  tabnla(ir,lmm,it) = 2.D0*( tabnl(ir,lmm,it) + tabnl(ir+2,lmm,it)  &
&                                     - 2.d0*tabnl(ir+1,lmm,it) )
  tbfnla(ir,lmm,it) = 2.D0*( tbfnl(ir,lmm,it) + tbfnl(ir+2,lmm,it)  &
&                                     - 2.d0*tbfnl(ir+1,lmm,it) )
end do

!-----------------------------------------------------------------------
   lmxx = lmxx + 2*l + 1
   laxx = laxx + 1
end do jdo
end if
end do ldo
lmx(it) = lmxx
lax(it) = laxx

!=======================================================================
end if typeif
end do typedo


rmxnlx = 0.d0
do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
       do l = 1, lax(it)
          rmxnlx = max( rmxnlx, rmxnlc(l,it) )
       end do
   end if
end do




!------deallocate memory
deallocate( tabwv, betaq, betar, amatrx, vmatrx, ivmatrx,  &
& stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in nlkbpp' )

dealocmem = 8.d0 * ( mx1+1 + (noofdt+1)*2 + (noofdt+1)*(noofdt+1)  &
&                 + noofdt*(noofdt+1) )  &
& + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d3),' KB, deallocated (nlkbpp)'
if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d3),' KB, deallocated (nlkbpp)'

call deallocate_unit_number( jpl )


return
end




subroutine kingsm( nfile, myid, nodes,  &
& l, it, xrdel, egmax, eg4, llking, rking,  &
& vlocli ,xitgrd, rr, nvlcl, betaq, gdel, betar, rspdel, noofdt,  &
& amatrx, vmatrx, ivmatrx )
!-----------------------------------------------------------------------
!     filtering functions
!
!  input : vlocli : r * r * beta
!          xrdel  : integration step in r-space
!          egmax  : energy cutoff for wave function
!          eg4    : energy cutoff for modified function
!          llking : = .true.   King-Smith implementation is applied
!                     .false.  Gaussian filter function is used
!
! output : betaq  : modified beta function in q-space
!          betar  : modified beta function in r-space
!-----------------------------------------------------------------------
use outfile
#ifdef VECTOR
use param_for_Bessel
#endif
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile

character(1), dimension(7)   :: ms =  &
& (/ 'S', 'P', 'D', 'F', 'G', 'H', 'I' /)

integer :: nvlcl
real*8, dimension(0:nvlcl)  :: vlocli ,xitgrd, rr
integer :: noofdt
real*8, dimension(0:noofdt) :: betaq
real*8 :: gdel
real*8, dimension(0:noofdt) :: betar
real*8 :: rspdel
real*8,  dimension(0:noofdt,0:noofdt) :: amatrx
real*8,  dimension(noofdt,noofdt+1)   :: vmatrx
integer, dimension(noofdt)            :: ivmatrx
CHARACTER*50   FNAME
logical lout, llking
#ifdef VECTOR
logical lflag
save lflag
#endif
save lout
data lout / .false. /
#ifdef VECTOR
data lflag / .true. /

if( lflag ) call cnstts
lflag = .false.
#endif

if( lout .and. myid == 0 ) then
    call allocate_unit_number( jpl )
end if

gmax  = sqrt( egmax )
gmax4 = sqrt( eg4 )
gdel  = gmax4/dble(noofdt)
noofm = int( gmax/gdel )
noofmt = noofdt - noofm
do ig = 0, noofdt
   gq = gdel*dble(ig)
#ifndef VECTOR
   do ir = 0, nvlcl
      xitgrd(ir) = vlocli(ir) * sbessl( gq*rr(ir), l )
   end do
#else
 if( l.eq.0 ) then
   do ir = 0, nvlcl
      q = gq*rr(ir)
      IF( ABS(Q).LT.7.D-01 ) THEN
          XF12 = Q*Q
          BSL = 1.0D0 + XF12*( FJ02 + XF12*( FJ04 + XF12*( FJ06  &
&                + XF12*( FJ08 + XF12*FJ010 ) ) ) )
        ELSE
          BSL = SIN(Q)/Q
      end if
      xitgrd(ir) = vlocli(ir) * bsl
   end do
 else if( l.eq.1 ) then
   do ir = 0, nvlcl
      q = gq*rr(ir)
      IF( ABS(Q).LT.7.D-01 ) THEN
          XF12 = Q*Q
          BSL = Q*( FJ11 + XF12*( FJ13 + XF12*( FJ15  &
&                      + XF12*( FJ17 + XF12*FJ19 ) ) ) )
        ELSE
          BSL = ( SIN(Q)/Q - COS(Q) )/Q
      end if
      xitgrd(ir) = vlocli(ir) * bsl
   end do
 else if( l.eq.2 ) then
   do ir = 0, nvlcl
      q = gq*rr(ir)
      IF( ABS(Q).LT.7.D-01 ) THEN
          XF12 = Q*Q
          BSL = XF12*( FJ22 + XF12*( FJ24 + XF12*( FJ26  &
&                      + XF12*( FJ28 + XF12*FJ210 ) ) ) )
        ELSE
          BS0 = SIN(Q)/Q
          BS1 = ( BS0 - COS(Q) )/Q
          BSL = 3.D0*BS1/Q - BS0
      end if
      xitgrd(ir) = vlocli(ir) * bsl
   end do
 else if( l.eq.3 ) then
   do ir = 0, nvlcl
      q = gq*rr(ir)
      IF( ABS(Q).LT.7.D-01 ) THEN
          XF12 = Q*Q
          BSL = Q*XF12*( FJ33 + XF12*( FJ35 + XF12*( FJ37  &
&                           + XF12*( FJ39 + XF12*FJ311 ) ) ) )
        ELSE
          BS0 = SIN(Q)/Q
          BS1 = ( BS0 - COS(Q) )/Q
          BS2 = 3.D0*BS1/Q - BS0
          BSL = 5.D0*BS2/Q - BS1
      end if
      xitgrd(ir) = vlocli(ir) * bsl
   end do
 end if
#endif
   call INTGB3( nvlcl, xrdel, xitgrd, CL )
   betaq(ig) = CL
end do


if( lout .and. myid == 0 ) then
!cc#ifdef SCALAR
!--- output the original beta function in r-space ----------------------
FNAME = ((( 'Beta_r_'//'_' )//MS(l+1) )//'.original' )
write(*,*) 'output to file : ', FNAME
OPEN(JPL, FILE = FNAME, status='unknown' )
write(JPL,'(a30)') '# beta function in r-space    '
write(JPL,'(a30,f10.6)') '#    original cutoff length : ',  &
&                        xrdel*dble(nvlcl)
write(JPL,'(a30,f10.6)') '#    modified cutoff length : ',  &
&                        rking
nskip = nvlcl/50
do ir = 1, nvlcl, nskip
   write(JPL,'(f10.6,e14.6)') rr(ir), vlocli(ir)/rr(ir)/rr(ir)
end do
CLOSE(JPL)

!--- output the original beta function in q-space ----------------------
FNAME = ((( 'Beta_q_'//'_' )//MS(l+1) )//'.original' )
write(*,*) 'output to file : ', FNAME
OPEN(JPL, FILE = FNAME, status='unknown' )
write(JPL,'(a30)') '# beta function in q-space    '
write(JPL,'(a35,2f10.6)') '#    original energy/w.v. cutoff : ',  &
&                        egmax, gmax
write(JPL,'(a35,2f10.6)') '#    modified energy/w.v. cutoff : ',  &
&                        eg4,   gmax4
write(JPL,'(a35,i6)') '#  # of q-vectors .le. orig. Ecut: ',  &
&                        noofm
write(JPL,'(a35,i6)') '#  # of q-vectors .gt. orig. Ecut: ',  &
&                        noofmt
do ig = 0, noofdt
   write(JPL,'(f10.6,e14.6)') gdel*dble(ig), betaq(ig)
end do
CLOSE(JPL)
!cc#endif
!-----------------------------------------------------------------------
end if


pi     = acos(-1.d0)
pih    = 0.5d0*pi
if( llking ) then
!=== King-Smith implementention ========================================
!--- set A matrix : amatrx
!      if( myid.eq.0 ) then
!      write(*,*) 'set A matrix'
!      end if

r0  = rking
r02 = r0*r0
r03 = r02*r0

do ig1 = 0, noofdt
   q1  = gdel*dble(ig1)
   q12 = q1*q1
   sbl1 = sbessl(q1*r0,l)
   sbp1 = sbessl(q1*r0,l+1)
   amatrx(ig1,ig1) = q12*q12*0.5d0*r03*( sbl1*sbl1  &
&                                      - sbessl(q1*r0,l-1)*sbp1 )
   do ig2 = 0, ig1 - 1
      q2  = gdel*dble(ig2)
      q22 = q2*q2
      sbl2 = sbessl(q2*r0,l)
      sbp2 = sbessl(q2*r0,l+1)
      amatrx(ig1,ig2) = q12*q22*r02/(q12-q22) * ( q1*sbl2*sbp1  &
&                                               - q2*sbl1*sbp2 )
      amatrx(ig2,ig1) = amatrx(ig1,ig2)
   end do
end do


!--- set matrix equations
!      if( myid.eq.0 ) then
!      write(*,*) 'set matrix equations'
!      end if

nofmt1 = noofmt + 1
do m1 = 1, noofmt
   ntrue = noofm + m1

!   --- already known quantities ---
   vmatrx(m1, nofmt1) = 0.5d0*amatrx(ntrue,0)*betaq(0)
   do m2 = 1, noofm
      vmatrx(m1, nofmt1) = vmatrx(m1, nofmt1)  &
&                        + amatrx(ntrue,m2)*betaq(m2)
   end do
   vmatrx(m1, nofmt1) = gdel*vmatrx(m1, nofmt1)

!   --- coefficients ---
   q1  = gdel*dble(ntrue)
   q12 = q1*q1
   do m2 = 1, noofmt - 1
      vmatrx(m1, m2) = -gdel*amatrx(ntrue, noofm+m2)
   end do
   m2 = noofmt
   vmatrx(m1, m2) = -0.5d0*gdel*amatrx(ntrue, noofm+m2)
   vmatrx(m1, m1) = vmatrx(m1, m1) + pih*q12

end do


!--- solve matrix equations
!      if( myid.eq.0 ) then
!      write(*,*) 'solve matrix equations'
!      end if

INDER = 0
call GSSJOR( vmatrx, ivmatrx, noofmt, noofdt, noofdt+1, INDER )
if( INDER.ne.0 ) then
    if(loutfile(1)) write(nfile(1),*) 'error in kingsm (GSSJOR) : ', INDER
    if(loutfile(2)) write(nfile(2),*) 'error in kingsm (GSSJOR) : ', INDER
!-----Finalize the parallel environment
    call end_parallel(ierr)
    stop
end if

do m1 = 1, noofmt
   ntrue = noofm + m1
   betaq(ntrue) = vmatrx(m1, nofmt1)
end do
!=== end of King-Smith implementention =================================
else
!=== filtering with Gaussian function in k-space =======================
    cfact  = 18.d0/((gmax4-gmax)**2)
    do ig = 0, noofdt
       gq = gdel*dble(ig)
       if( gq.gt.gmax )  &
&          betaq(ig) = betaq(ig)*dexp(-cfact*(gq-gmax)**2)
    end do
!=== end of filtering with Gaussian function in k-space ================
end if


!--- modified beta function in r-space ---
rspdel = rking/dble(noofdt)
do ir = 0, noofdt
   r = rspdel*dble(ir)
   do ig1 = 0, noofdt
      q1  = gdel*dble(ig1)
      xitgrd(ig1) = q1*q1*betaq(ig1) * sbessl( q1*r, l )
   end do
   call INTGB3( noofdt, gdel, xitgrd, CL )
   betar(ir) = CL/pih
end do


if( .not.llking ) then
!=== filtering with Gaussian function in r-space =======================
    cfact  = 18.d0/((rking-xrdel*dble(nvlcl))**2)
    do ir = 0, noofdt
       r = rspdel*dble(ir)
       if( r.gt.xrdel*dble(nvlcl) )  &
&          betar(ir) =  &
&                 betar(ir)*dexp(-cfact*(r-xrdel*dble(nvlcl))**2)
    end do
!=== end of filtering with Gaussian function in r-space ================
end if


if( lout .and. myid.eq.0 ) then
!ccc#ifdef SCALAR
!--- output the modified beta function in r-space ----------------------
FNAME = ((( 'Beta_r_'//'_' )//MS(l+1) )//'.modified' )
write(*,*) 'output to file : ', FNAME
OPEN(JPL, FILE = FNAME, status='unknown' )
write(JPL,'(a30)') '# beta function in r-space    '
write(JPL,'(a30,f10.6)') '#    original cutoff length : ',  &
&                        xrdel*dble(nvlcl)
write(JPL,'(a30,f10.6)') '#    modified cutoff length : ',  &
&                        rking
do ir = 0, noofdt
   r = rspdel*dble(ir)
   write(JPL,'(f10.6,e14.6)') r, betar(ir)
end do
CLOSE(JPL)

!--- output the modified beta function in q-space ----------------------
FNAME = ((( 'Beta_q_'//'_' )//MS(l+1) )//'.modified' )
write(*,*) 'output to file : ', FNAME
OPEN(JPL, FILE = FNAME, status='unknown' )
write(JPL,'(a30)') '# beta function in q-space    '
write(JPL,'(a35,2f10.6)') '#    original energy/w.v. cutoff : ',  &
&                        egmax, gmax
write(JPL,'(a35,2f10.6)') '#    modified energy/w.v. cutoff : ',  &
&                        eg4,   gmax4
write(JPL,'(a35,i6)') '#  # of q-vectors l.e. orig. Ecut: ',  &
&                        noofm
write(JPL,'(a35,i6)') '#  # of q-vectors g.t. orig. Ecut: ',  &
&                        noofmt
do ig = 0, noofdt
   write(JPL,'(f10.6,e14.6)') gdel*dble(ig), betaq(ig)
end do
CLOSE(JPL)
!ccc#endif
!-----------------------------------------------------------------------
end if

if( lout .and. myid == 0 ) then
    call deallocate_unit_number( jpl )
end if

return
end
