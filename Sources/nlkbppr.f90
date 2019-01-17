



module nlkbppr_variables
!-----------------------------------------------------------------------
! type declaration of shared variables for nonlocal pp.
!-----------------------------------------------------------------------
implicit none

!integer :: nknod = 1   ! the number of k points
!logical :: lgamma
!logical :: lnoncollinear

!------- variables for nonlocal pp.
real*8,  allocatable, dimension(:,:) :: sumr, sumi
real*8,  allocatable, dimension(:,:) :: subr, subi

#ifndef VECTOR
real*8,  allocatable, dimension(:) :: vnlc
real*8,  allocatable, dimension(:) :: fvnlc
#endif

integer, allocatable, dimension(:) :: meshnl
real*8,  allocatable, dimension(:) ::  vnlcr

integer :: nionnew
integer, allocatable, dimension(:) :: iatmpt     ! atom index to calculate nonlocal pp
real*8,  allocatable, dimension(:,:) :: cratm
integer, allocatable, dimension(:,:) :: numsnl, npmsnl
integer, allocatable, dimension(:,:) :: npvnlc
integer,  allocatable, dimension(:,:) :: ncratm_disp
#ifdef VECTOR
integer, allocatable, dimension(:,:,:) :: ncornl, nwidnl
integer, allocatable, dimension(:,:,:) :: ndv4nl
#endif
real*8,  dimension(3) :: xc

save

end module




module force_nlc_variables
!-----------------------------------------------------------------------
! type declaration of shared variables for nonlocal pp forces
!-----------------------------------------------------------------------
implicit none

real*8, allocatable, dimension(:,:,:) :: ylmr, ylmi
#ifdef FAST_NONLOCAL_FORCE
real*8, allocatable, dimension(:) :: fvnlxr, fvnlyr, fvnlzr
real*8, allocatable, dimension(:,:) :: sumxr, sumyr, sumzr,  &
&                 sumxx, sumyy, sumzz, sumyz, sumzx, sumxy
#else
real*8, allocatable, dimension(:,:) :: sumxr, sumxx, sumxy
#endif
real*8, allocatable, dimension(:,:) :: sumxi, sumxxi, sumxyi

save

end module




subroutine nlkbppr_prealloc( nfile, myid, nodes,  &
& alloc_mem, lnref )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: lnref

!-----declare local variables
integer :: status
real*8  :: the_mem


#ifndef VECTOR
!------allocate memory
allocate( vnlc(lnref), fvnlc(lnref), stat=status )

the_mem =  &
& + 8.d0 * ( lnref*2 )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlkbppr_prealloc', .true. )
#endif


return
end subroutine




subroutine nlkbppr_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, mxl, nylmmx, nbnod, nspnmx )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: natom, mxl, nylmmx, nbnod
integer :: nspnmx

!-----declare local variables
logical :: lgamma, lnoncollinear
integer :: status
real*8  :: the_mem


!---get No. of k points
!call get_nknod( nknod )
call get_lgamma( lgamma )

call get_lnoncollinear( lnoncollinear )

!------allocate memory
allocate( sumr(nylmmx,natom), stat=status )

the_mem = 8.d0 * ( size(sumr) )

if( status == 0 ) then
    if( .not.lgamma .or. lnoncollinear ) then
        !------allocate memory
        allocate( sumi(nylmmx,natom), stat=status )

        the_mem = the_mem + 8.d0 * ( size(sumi) )

        sumi(:,:) = 0.d0
    end if
end if

if( status == 0 ) then
    !-----noncollinear magnetism
    if( lnoncollinear ) then
        !------allocate memory
        allocate( subr(nylmmx,natom), subi(nylmmx,natom), stat=status )

        the_mem = the_mem + 8.d0 * ( size(subr) + size(subi) )

        subr(:,:) = 0.d0
        subi(:,:) = 0.d0
    end if
end if


!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlkbppr_variables_alloc', .true. )


!-----allocate memory for variables for nonlocal pp forces
call force_nlc_alloc( nfile, myid, nodes,  &
& alloc_mem, mxl, nylmmx, natom )


return
end subroutine




subroutine vpp_ncl_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nhk1, nhk2, ntall, natom, lnref, lmmax,  &
& lkbppi, lking, lgamma )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
integer :: ntall
integer :: natom
integer :: lnref
integer :: lmmax
logical, dimension(ntype) :: lkbppi, lking
logical :: lgamma

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: i, it, ii


!------allocate memory
allocate( iatmpt(ntall), cratm(3,ntall),  &
& numsnl(lnref,ntall), npmsnl(lnref,ntall), npvnlc(lmmax,ntall),  &
#ifdef VECTOR
& ncornl(3,lnref,ntall), nwidnl(3,lnref,ntall),  &
& ndv4nl(3,lnref,ntype),  &
#endif
& stat=status )

the_mem =  &
& + 4.d0 * ( ntall )  &
& + 8.d0 * ( 3*ntall )  &
#ifdef VECTOR
& + 4.d0 * ( 3*lnref*ntall*2 + 3*lnref*ntype )  &
#endif
& + 4.d0 * ( lnref*ntall*2 + lmmax*ntall )

if( status == 0 .and. .not.lgamma ) then
    allocate( ncratm_disp(3,ntall),  &
& stat=status )
    the_mem = the_mem + 4.d0 * ( 3*ntall )
end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'vpp_ncl_alloc', .true. )


ii = 0
do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
       do i = nhk1(it), nhk2(it)
          ii = ii + 1
          iatmpt(ii) = i
       end do
   end if
end do


return
end subroutine




subroutine vpp_ncl_dealloc( nfile, myid, nodes,  &
& alloc_mem, ntype, ntall, natom, lnref, lmmax )
!-----------------------------------------------------------------------
!     deallocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer :: ntall
integer :: natom
integer :: lnref
integer :: lmmax

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: i


!------allocate memory
deallocate( iatmpt, cratm, numsnl, npmsnl, npvnlc,  &
#ifdef VECTOR
& ncornl, nwidnl, ndv4nl,  &
#endif
& stat=status )

the_mem =  &
& + 4.d0 * ( ntall )  &
& + 8.d0 * ( 3*ntall )  &
#ifdef VECTOR
& + 4.d0 * ( 3*lnref*ntall*2 + 3*lnref*ntype )  &
#endif
& + 4.d0 * ( lnref*ntall*2 + lmmax*ntall )

if( status == 0 .and. allocated(ncratm_disp) ) then
    deallocate( ncratm_disp,  &
& stat=status )
    the_mem = the_mem + 4.d0 * ( 3*ntall )
end if

!------error trap
call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'vpp_ncl_dealloc', .true. )


return
end subroutine




subroutine nlcmx_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1, nlcmx2 )
!-----------------------------------------------------------------------
!     allocate memory for variables for related to nlcmx
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nlcmx1, nlcmx2

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( meshnl(nlcmx2), vnlcr(nlcmx1), stat=status )

the_mem = 4.d0 * nlcmx2 + 8.d0 * nlcmx1

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlcmx_alloc', .true. )


!-----allocate memory for variables for nonlocal pp forces
call force_nlcmx_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1 )


return
end subroutine




subroutine nlcmx_dealloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1, nlcmx2 )
!-----------------------------------------------------------------------
!     deallocate memory for variables for related to nlcmx
!-----------------------------------------------------------------------
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nlcmx1, nlcmx2

!-----declare local variables
integer :: status
real*8  :: the_mem


!------deallocate memory
deallocate( meshnl, vnlcr, stat=status )

the_mem = 4.d0 * nlcmx2 + 8.d0 * nlcmx1

!------error trap
call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlcmx_dealloc', .true. )


!-----allocate memory for variables for nonlocal pp forces
call force_nlcmx_dealloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1 )


return
end subroutine




subroutine force_nlc_alloc( nfile, myid, nodes,  &
& alloc_mem, mxl, nylmmx_, natom )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp forces
!-----------------------------------------------------------------------
use nlpp_variables
use force_nlc_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: mxl
integer :: nylmmx_
integer :: natom

!-----declare local variables
!logical :: lgamma
integer :: status
real*8  :: the_mem


!call get_lgamma( lgamma )

!------allocate memory
allocate( ylmr(3,-mxl:mxl,0:mxl), ylmi(3,-mxl:mxl,0:mxl),  &
#ifdef FAST_NONLOCAL_FORCE
& sumxr(nylmmx,natom), sumyr(nylmmx,natom), sumzr(nylmmx,natom),  &
& sumxx(nylmmx,natom), sumyy(nylmmx,natom), sumzz(nylmmx,natom),  &
& sumyz(nylmmx,natom), sumzx(nylmmx,natom), sumxy(nylmmx,natom),  &
#else
& sumxr(nylmmx,natom), sumxx(nylmmx,natom), sumxy(nylmmx,natom),  &
#endif
& stat=status )

the_mem =  &
&  8.d0 * ( 3*(2*mxl+1)*(mxl+1)*2  &
#ifdef FAST_NONLOCAL_FORCE
&         + nylmmx*natom*9 )
#else
&         + nylmmx*natom*3 )
#endif

if( status == 0 ) then
    if( .not.lgamma .or. lnoncollinear ) then
        !------allocate memory
        allocate(  &
& sumxi(nylmmx,natom), sumxxi(nylmmx,natom), sumxyi(nylmmx,natom),  &
&  stat=status )

        the_mem = the_mem + 8.d0 * ( size(sumxi) + size(sumxxi) + size(sumxyi) )
    end if
end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'force_nlc_alloc', .true. )


return
end subroutine




subroutine force_nlcmx_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx )
!-----------------------------------------------------------------------
!     allocate memory for variables related to nlcmx
!-----------------------------------------------------------------------
use force_nlc_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nlcmx

!-----declare local variables
integer :: status
real*8  :: the_mem


#ifdef FAST_NONLOCAL_FORCE
!------allocate memory
allocate( fvnlxr(nlcmx), fvnlyr(nlcmx), fvnlzr(nlcmx),  &
& stat=status )

the_mem = 8.d0 * nlcmx*3

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'force_nlcmx_alloc', .true. )
#endif


return
end subroutine




subroutine force_nlcmx_dealloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx )
!-----------------------------------------------------------------------
!     deallocate memory for variables related to nlcmx
!-----------------------------------------------------------------------
use force_nlc_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nlcmx

!-----declare local variables
integer :: status
real*8  :: the_mem


#ifdef FAST_NONLOCAL_FORCE
!------allocate memory
deallocate( fvnlxr, fvnlyr, fvnlzr, stat=status )

the_mem = 8.d0 * nlcmx*3

!------error trap
call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'force_nlcmx_dealloc', .true. )
#endif


return
end subroutine




subroutine vlocal( nfile, myid, nodes,  &
& lclust, lorthrhmbc, nd1v, ntype, nhk1, nhk2, nhk1_nod, nhk2_nod,  &
& zv, nel, ratm,  &
& iatoit, iatmpt, nion, natom, llking, mx1,  &
& hcell, rdel, lhfull, vext, tablc, tablca, dltlc, rmxlc,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, tmpk, ltimecnt )
!-----------------------------------------------------------------------
!    local pseudopotential : vext
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
logical :: lclust
logical :: lorthrhmbc
integer, dimension(3) :: nd1v
integer :: ntype
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(ntype) :: nhk1_nod, nhk2_nod
real*8,  dimension(3,natom) :: ratm
integer, dimension(natom)   :: iatoit
integer, dimension(nion)    :: iatmpt
integer :: nion
logical, dimension(ntype) :: llking
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3) :: rdel
logical :: lhfull
real*8, dimension(0:mx1,ntype) :: tablc, tablca
real*8, dimension(ntype)       :: dltlc, rmxlc

!--- for bulk calculation ----------------------------------------
real*8, dimension(8) :: sbl1 =  &
&             (/ 1.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0 /)
real*8, dimension(8) :: sbl2 =  &
&             (/ 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0 /)
real*8, dimension(8) :: sbl3 =  &
&             (/ 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0 /)
dimension      mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
!-----------------------------------------------------------------
dimension vext(mshnod)
logical        ltimecnt

dimension      mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension tmpk(*)
save sbl1, sbl2, sbl3


if( ltimecnt ) then
    ct0 = timecnt()
end if
clustif: if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
#ifdef VECTOR
   do it = 1, ntype
      zvit2  = 2.d0*zv(it)
      rmxlc2 = rmxlc(it)*rmxlc(it)
   do i = nhk1(it), nhk2(it)
#endif
do m1 = 1, mshnod
   ix = mshnx(m1)
   ix = mshx1 + ix - 1
   iy = mshny(m1)
   iy = mshy1 + iy - 1
   iz = mshnz(m1)
   iz = mshz1 + iz - 1
   xm = dble(ix-1)*rdel(1)
   ym = dble(iy-1)*rdel(2)
   zm = dble(iz-1)*rdel(3)
#ifndef VECTOR
   do it = 1, ntype
      zvit2  = 2.d0*zv(it)
      rmxlc2 = rmxlc(it)*rmxlc(it)
   do i = nhk1(it), nhk2(it)
#endif
      xx = xm - ratm(1,i)
      yy = ym - ratm(2,i)
      zz = zm - ratm(3,i)
      r2 = xx*xx + yy*yy + zz*zz
      if( r2.le.rmxlc2 ) then
          m = r2/(2.d0*dltlc(it))
          m = 2*m
          d = 0.5d0*( r2/dltlc(it) - dble(m) )
          vloc = d*( (d-1.d0)*tablca(m,it) + tablc(m+2,it)  &
&                          - tablc(m,it) ) + tablc(m,it)
          vext(m1) = vext(m1) + vloc
      else
          r = sqrt(r2)
          vext(m1) = vext(m1) - zvit2/r
      end if
   end do
   end do
end do
else clustif
!--- for bulk calculation ----------------------------------------------
if( lhfull ) then
    ihroop = 8
    dhalf  = 0.d0
  else
    ihroop = 1
    dhalf  = 0.5d0
end if
orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
#ifdef VECTOR
   typedo: do it = 1, ntype
   typeif: if( llking(it) ) then
      rmxlc2 = rmxlc(it)*rmxlc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
   do m1 = 1, mshnod
      ix = mshnx(m1)
      ix = mshx1 + ix - 1
      iy = mshny(m1)
      iy = mshy1 + iy - 1
      iz = mshnz(m1)
      iz = mshz1 + iz - 1
      xm = dble(ix-1)/dble(nd1v(1))
      ym = dble(iy-1)/dble(nd1v(2))
      zm = dble(iz-1)/dble(nd1v(3))
#ifndef VECTOR
      typedo: do it = 1, ntype
      typeif: if( llking(it) ) then
         rmxlc2 = rmxlc(it)*rmxlc(it)
      do i = nhk1(it), nhk2(it)
      do lloop = 1, ihroop
#endif
         q1 = xm - ratm(1,i)
         q2 = ym - ratm(2,i)
         q3 = zm - ratm(3,i)
         if( abs(q1).ge.dhalf ) q1 = q1 - sign(sbl1(lloop),q1)
         if( abs(q2).ge.dhalf ) q2 = q2 - sign(sbl2(lloop),q2)
         if( abs(q3).ge.dhalf ) q3 = q3 - sign(sbl3(lloop),q3)
         xx = hcell(1,1)*q1
         yy = hcell(2,2)*q2
         zz = hcell(3,3)*q3
         r2 = xx*xx + yy*yy + zz*zz
!#ifdef VECTOR
!               tmpk(m1) = r2
!         end do
!         do m1 = 1, mshnod
!               r2 = tmpk(m1)
!#endif
         if( r2.le.rmxlc2 ) then
             m = r2/(2.d0*dltlc(it))
             m = 2*m
             d = 0.5d0*( r2/dltlc(it) - dble(m) )
             vloc = d*( (d-1.d0)*tablca(m,it) + tablc(m+2,it)  &
&                             - tablc(m,it) ) + tablc(m,it)
             vext(m1) = vext(m1) + vloc
         end if
#ifndef VECTOR
      end do
      end do
      end if typeif
      end do typedo
#endif
   end do
#ifdef VECTOR
   end do
   end do
   end if typeif
   end do typedo
#endif

else orthoif

!------non-orthorhombic super cell -------------------------------
#ifdef VECTOR
   typedo2: do it = 1, ntype
   typeif2: if( llking(it) ) then
      rmxlc2 = rmxlc(it)*rmxlc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
   do m1 = 1, mshnod
      ix = mshnx(m1)
      ix = mshx1 + ix - 1
      iy = mshny(m1)
      iy = mshy1 + iy - 1
      iz = mshnz(m1)
      iz = mshz1 + iz - 1
      xm = dble(ix-1)/dble(nd1v(1))
      ym = dble(iy-1)/dble(nd1v(2))
      zm = dble(iz-1)/dble(nd1v(3))
#ifndef VECTOR
      typedo2: do it = 1, ntype
      typeif2: if( llking(it) ) then
         rmxlc2 = rmxlc(it)*rmxlc(it)
      do i = nhk1(it), nhk2(it)
      do lloop = 1, ihroop
#endif
         q1 = xm - ratm(1,i)
         q2 = ym - ratm(2,i)
         q3 = zm - ratm(3,i)
         if( abs(q1).ge.dhalf ) q1 = q1 - sign(sbl1(lloop),q1)
         if( abs(q2).ge.dhalf ) q2 = q2 - sign(sbl2(lloop),q2)
         if( abs(q3).ge.dhalf ) q3 = q3 - sign(sbl3(lloop),q3)
         xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
         yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
         zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
         r2 = xx*xx + yy*yy + zz*zz
!#ifdef VECTOR
!               tmpk(m1) = r2
!         end do
!         do m1 = 1, mshnod
!               r2 = tmpk(m1)
!#endif
         if( r2.le.rmxlc2 ) then
             m = r2/(2.d0*dltlc(it))
             m = 2*m
             d = 0.5d0*( r2/dltlc(it) - dble(m) )
             vloc = d*( (d-1.d0)*tablca(m,it) + tablc(m+2,it)  &
&                             - tablc(m,it) ) + tablc(m,it)
             vext(m1) = vext(m1) + vloc
         end if
#ifndef VECTOR
      end do
      end do
      end if typeif2
      end do typedo2
#endif
   end do
#ifdef VECTOR
   end do
   end do
   end if typeif2
   end do typedo2
#endif

end if orthoif
!--- end of local pseudopotential in real space  -----------------------
end if clustif
if( ltimecnt ) then
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '         real space : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '         real space : cpu-time :', ct - ct0
    ct0 = ct
end if

if( .not.lclust ) then
!          call stscfn(
!     & ntype, ratm, iatoit, iatmpt, natom, nion )
!          if( ltimecnt ) then
!              call gsync
!              ct = timecnt()
!              ct1 = ct - ct0
!              ct0 = ct
!          end if
!--- set sine, cosine function for mesh points
    call stscms( nfile, myid, nodes,  &
& nd1v, mshx1, mshy1, mshz1, mshx, mshy, mshz )
    if( ltimecnt ) then
        call gsync
        ct = timecnt()
        ct2 = ct - ct0
        ct0 = ct
    end if
!--- reciprocal space contribution
    nelv = 0
    do it = 1, ntype
       if( llking(it) ) then
           nelv = nelv + nint(zv(it))*( nhk2(it) - nhk1(it) + 1 )
       end if
    end do
    call revloc(  &
& ntype, nhk1_nod, nhk2_nod, zv, nelv, iatoit, iatmpt, natom, nion,  &
& llking, vext, mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, mshnod )
    if( ltimecnt ) then
        call gsync
        ct = timecnt()
        do ii = 1, 2
        if( loutfile(ii) ) then
!              write(nfile(ii),*) '             stscfn :          :', ct1
            write(nfile(ii),*) '             stscms :          :', ct2
            write(nfile(ii),*) '             revloc :          :', ct - ct0
        end if
        end do
        ct0 = ct
    end if
end if


return
end




subroutine cnanpp( nfile, myid, nodes,  &
& ntall, lclust, lvacuum, ntype, nhk1, nhk2, ratm, lkbppi, lking,  &
& rmxnlx, hcell, b, xc )
!-----------------------------------------------------------------------
!     count the number of atoms to calculate nonlocal pp. in real space
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: ntall
logical :: lclust
logical, dimension(3) :: lvacuum
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,*) :: ratm
logical, dimension(ntype) :: lkbppi, lking
real*8 :: rmxnlx
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3,3) :: b
real*8,  dimension(3) :: xc

!-----declare local variables
integer :: i, j, ix, iy, iz, npp, it, mltx, mlty, mltz
real*8 :: expand = 1.1d0    ! 


ntall = 0
do it = 1, ntype
   if( lkbppi(it) .and. lking(it) )  &
&      ntall = ntall + nhk2(it) - nhk1(it) + 1
end do

!---call rciprl( hcell, b, v )

do i = 1, 3
   xc(i) = sqrt( b(1,i)*b(1,i) + b(2,i)*b(2,i) + b(3,i)*b(3,i) )  &
&     * rmxnlx
   xc(i) = xc(i) * expand
end do
if( lclust ) return


!if( xc(1) >= 1.d0 .or. xc(2) >= 1.d0 .or. xc(3) >= 1.d0 ) then
!    ntall = ntall  &
!&       * ( 2*xc(1) + 1.d0 )*( 2*xc(2) + 1.d0 )*( 2*xc(3) + 1.d0 )
!    return
!end if


typedo: do it = 1, ntype
typeif: if( lkbppi(it) .and. lking(it) ) then

do ix = 1, 3
do mltx = 0, int(xc(ix))
   if( .not.lvacuum(ix) ) then
   do j = nhk1(it), nhk2(it)
      if(        ratm(ix,j) + dble(mltx) <= xc(ix) ) ntall = ntall + 1
      if( 1.d0 - ratm(ix,j) + dble(mltx) <= xc(ix) ) ntall = ntall + 1
   end do
   end if
end do
end do

do ix = 1, 2
do mltx = 0, int(xc(ix))
do iy = ix + 1, 3
do mlty = 0, int(xc(iy))
   if( .not.lvacuum(ix) .and. .not.lvacuum(iy) ) then
   do j = nhk1(it), nhk2(it)
      if(        ratm(ix,j) + dble(mltx) <= xc(ix) ) then
         if(        ratm(iy,j) + dble(mlty) <= xc(iy) ) ntall = ntall + 1
         if( 1.d0 - ratm(iy,j) + dble(mlty) <= xc(iy) ) ntall = ntall + 1
      end if
      if( 1.d0 - ratm(ix,j) + dble(mltx) <= xc(ix) ) then
         if(        ratm(iy,j) + dble(mlty) <= xc(iy) ) ntall = ntall + 1
         if( 1.d0 - ratm(iy,j) + dble(mlty) <= xc(iy) ) ntall = ntall + 1
      end if
   end do
   end if
end do
end do
end do
end do

ix = 1
iy = 2
iz = 3
do mltx = 0, int(xc(ix))
do mlty = 0, int(xc(iy))
do mltz = 0, int(xc(iz))
   if( .not.lvacuum(ix) .and. .not.lvacuum(iy) .and.  &
&      .not.lvacuum(iz) ) then
   do j = nhk1(it), nhk2(it)
      if(        ratm(ix,j) + dble(mltx) <= xc(ix) ) then
         if(        ratm(iy,j) + dble(mlty) <= xc(iy) ) then
            if(        ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
            if( 1.d0 - ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
         end if
         if( 1.d0 - ratm(iy,j) + dble(mlty) <= xc(iy) ) then
            if(        ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
            if( 1.d0 - ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
         end if
      end if
      if( 1.d0 - ratm(ix,j) + dble(mltx) <= xc(ix) ) then
         if(        ratm(iy,j) + dble(mlty) <= xc(iy) ) then
            if(        ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
            if( 1.d0 - ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
         end if
         if( 1.d0 - ratm(iy,j) + dble(mlty) <= xc(iy) ) then
            if(        ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
            if( 1.d0 - ratm(iz,j) + dble(mltz) <= xc(iz) ) ntall = ntall + 1
         end if
      end if
   end do
   end if
end do
end do
end do

end if typeif
end do typedo


return
end subroutine




#ifdef VECTOR
module nlcvector
!-----------------------------------------------------------------------
! type declaration of work variables for nonlocal pp, only for VECTOR
!-----------------------------------------------------------------------
implicit none

real*8, parameter :: pi = 3.141592653589793d0,  &
&                  pi00 = 2.820947917738781d-01,  &
&                  pi10 = 4.886025119029199d-01,  &
&                  pi11 = 3.454941494713355d-01,  &
&                  pi20 = 3.153915652525200d-01,  &
&                  pi21 = 7.725484040463791d-01,  &
&                  pi22 = 3.862742020231896d-01,  &
&                  pi30 = 3.731763325901154d-01,  &
&                  pi31 = 3.231801841141507d-01,  &
&                  pi32 = 1.021985476433282d0  ,  &
&                  pi33 = 4.172238236327841d-01

real*8,  allocatable, dimension(:) ::  &
& xmr, ymr, zmr, g1r, vnlcrr, xxnl,  &
& qqa, qqb, qqc, qqd, qqe, qqf, qqg, qqh, qqi,  &
& qqj, qqk, qql, qqm, qqn, qqo, qqp, qqq, qqr, qqs, qqt, qqu

end module




subroutine nlcvector_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx3 )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlcvector
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nlcmx3

!-----declare local variables
integer :: status
real*8  :: the_mem


!-----if already allocated, deallocate arrays
if( allocated(xmr) ) then
    the_mem = 8.d0 * size(xmr) * 27.d0

    !------deallocate memory
    deallocate(  &
& xmr, ymr, zmr, g1r, vnlcrr, xxnl,  &
& qqa, qqb, qqc, qqd, qqe, qqf, qqg, qqh, qqi,  &
& qqj, qqk, qql, qqm, qqn, qqo, qqp, qqq, qqr, qqs, qqt, qqu,  &
& stat=status )

    !------error trap
    call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlcvector_alloc', .true. )
end if


!------allocate memory
allocate(  &
& xmr(nlcmx3), ymr(nlcmx3), zmr(nlcmx3), g1r(nlcmx3), vnlcrr(nlcmx3), xxnl(nlcmx3),  &
& qqa(nlcmx3), qqb(nlcmx3), qqc(nlcmx3), qqd(nlcmx3), qqe(nlcmx3),  &
& qqf(nlcmx3), qqg(nlcmx3), qqh(nlcmx3), qqi(nlcmx3), qqj(nlcmx3),  &
& qqk(nlcmx3), qql(nlcmx3), qqm(nlcmx3), qqn(nlcmx3), qqo(nlcmx3),  &
& qqp(nlcmx3), qqq(nlcmx3), qqr(nlcmx3), qqs(nlcmx3), qqt(nlcmx3), qqu(nlcmx3),  &
& stat=status )

the_mem = 8.d0 * size(xmr) * 27.d0

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlcvector_alloc', .true. )


return
end subroutine
#endif




subroutine setnlc( nfile, myid, nodes,  &
& lvacuum, lclust, ntype, nhk1, nhk2, lmax, lclno, lchk, lkbppi, lking,  &
& iatoit, rdel, rdelg, hcell, hci, lorthrhmbc, ratm,  &
& mshglb, kfft1, kfft2, kfft3,  &
& ylmr, ylmi, natom, mxl, nylmmx_, mx1, alloc_mem )
!-----------------------------------------------------------------------
!    set parameters for calculation of nonlocal pseudopotentials
!-----------------------------------------------------------------------
use outfile
#ifdef VECTOR
use nlcvector
#endif
use nlpp_variables
use nlppr_variables
use nlkbppr_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
logical, dimension(3) :: lvacuum
logical :: lclust
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lclno
integer :: natom, mxl, nylmmx_, mx1
logical, dimension(0:mxl,ntype) :: lchk
logical, dimension(ntype) :: lkbppi, lking
integer, dimension(natom) :: iatoit
real*8,  dimension(3)   :: rdel
real*8,  dimension(3,3) :: rdelg
logical :: lorthrhmbc
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3,3) :: hci
real*8,  dimension(3,natom)  :: ratm

integer :: kfft1, kfft2, kfft3
integer, dimension(kfft1,kfft2,kfft3) :: mshglb

real*8,  dimension(-mxl:mxl,0:mxl) :: ylmr, ylmi
real*8  :: alloc_mem

!-----declare local variables
integer, dimension(3) :: ndv3nl
integer :: nlcmx1 = 0, nlcmx2 = 0, nlcmx3 = 0, ntall = 0
save nlcmx1, nlcmx2, nlcmx3, ntall


!-----count the number of atoms to calculate nonlocal pp. in real space
call cnanpp( nfile, myid, nodes,  &
& ntanew, lclust, lvacuum, ntype, nhk1, nhk2, ratm, lkbppi, lking,  &
& rmxnlx, hcell, hci, xc )

!-----check array sizes
if( ntanew > ntall ) then

   !-----if already allocated, deallocate arrays
   if( allocated(iatmpt) ) then
       call vpp_ncl_dealloc( nfile, myid, nodes,  &
& alloc_mem, ntype, ntall, natom, lnref, lmmax )
   end if

   !-----allocate arrays
   ntall = ntanew * 1.05d0
   call vpp_ncl_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nhk1, nhk2, ntall, natom, lnref, lmmax,  &
& lkbppi, lking, lgamma )

end if


mshx1 = 1
mshy1 = 1
mshz1 = 1
mshx  = kfft1
mshy  = kfft2
mshz  = kfft3
!if( lgamma ) then
    call baset( nfile, myid, nodes,  &
& lvacuum, mshx1, mshy1, mshz1, mshx, mshy, mshz, lclust,  &
& ntype, nhk1, nhk2, iatmpt, nionnew, rdelg, hcell, ratm, cratm,  &
& lkbppi, lking, rmxnlx, xc, ntall )
!else
!    !---for k-point sampling
!    call baset_k( nfile, myid, nodes,  &
!& lvacuum, mshx1, mshy1, mshz1, mshx, mshy, mshz, lclust,  &
!& ntype, nhk1, nhk2, iatmpt, nionnew, rdelg, hcell, ratm, cratm,  &
!& lkbppi, lking, rmxnlx, xc, ntall, ncratm_disp )
!end if


!      lmxxxx = 0
!      do it = 1, ntype
!         lmxx = 0
!         do l = 0, lmax(it)
!         if( lchk(l,it) .and. l.ne.lclno(it) ) then
!             lmxx = lmxx + 2*l + 1
!         end if
!         end do
!         lmxxxx = max( lmxxxx, lmxx )
!      end do


rmxnl2 = rmxnlx*rmxnlx
do i = 1, 3
!   ndv3nl(i) = xc(i)*hcell(i,i)/rdel(i) + 1.d0
   ndv3nl(i) = xc(i)*vecratio(hcell(1,i),rdelg(1,i)) + 1.d0
end do
#ifdef VECTOR
do it = 1, ntype
if( lkbppi(it) .and. lking(it) ) then
    do la = 1, lax(it)
       do i = 1, 3
!          ndv4nl(i,l,it) = xc(i)*hcell(i,i)/rdel(i)  &
          ndv4nl(i,la,it) = xc(i)*vecratio(hcell(1,i),rdelg(1,i))  &
&                        * rmxnlc(la,it)/rmxnlx + 1.d0
       end do
    end do
end if
end do
#endif

small=1.0d-08
sqrtwo = sqrt(2.d0)

!dn1v1 = hcell(1,1)/rdel(1)
!dn1v2 = hcell(2,2)/rdel(2)
!dn1v3 = hcell(3,3)/rdel(3)
dn1v1 = vecratio(hcell(1,1),rdelg(1,1))
dn1v2 = vecratio(hcell(1,2),rdelg(1,2))
dn1v3 = vecratio(hcell(1,3),rdelg(1,3))
!--- pre-calculation to count the number of mesh points around each atom
do i = 1, nionnew
do la = 1, lnref
   numsnl(la,i) = 0
end do
end do
do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)
#ifdef VECTOR
   do la = 1, lax(it)
#endif
   xx = cratm(1,i_nod)
   yy = cratm(2,i_nod)
   zz = cratm(3,i_nod)
!         nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
!         nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
!         nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
   nxx = ( hci(1,1)*xx + hci(2,1)*yy + hci(3,1)*zz )*dn1v1 + 0.5d0
   nyy = ( hci(1,2)*xx + hci(2,2)*yy + hci(3,2)*zz )*dn1v2 + 0.5d0
   nzz = ( hci(1,3)*xx + hci(2,3)*yy + hci(3,3)*zz )*dn1v3 + 0.5d0
   nxx = nxx + 1
   nyy = nyy + 1
   nzz = nzz + 1
#ifdef VECTOR
   ncornl(1,la,i_nod) = max(1,nxx-ndv4nl(1,la,it))
   ncornl(2,la,i_nod) = max(1,nyy-ndv4nl(2,la,it))
   ncornl(3,la,i_nod) = max(1,nzz-ndv4nl(3,la,it))
   nwidnl(1,la,i_nod) = min(mshx,nxx+ndv4nl(1,la,it))
   nwidnl(2,la,i_nod) = min(mshy,nyy+ndv4nl(2,la,it))
   nwidnl(3,la,i_nod) = min(mshz,nzz+ndv4nl(3,la,it))
#ifndef SR8000
   numsnl(la,i_nod) =  &
&              max( nwidnl(1,la,i_nod) - ncornl(1,la,i_nod) + 1, 0 )  &
&            * max( nwidnl(2,la,i_nod) - ncornl(2,la,i_nod) + 1, 0 )  &
&            * max( nwidnl(3,la,i_nod) - ncornl(3,la,i_nod) + 1, 0 )
#else
   numl = 0
   do iz = ncornl(3,la,i_nod), nwidnl(3,la,i_nod) 
   do iy = ncornl(2,la,i_nod), nwidnl(2,la,i_nod) 
   do ix = ncornl(1,la,i_nod), nwidnl(1,la,i_nod) 
!         m1 = ix + ( (iy-1) + (iz-1)*mshy )*mshx
!c         m1 = mshglb(ix,iy,iz)
!c         if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
!             if( g1.le.rmxnl2 ) then

           g1   = sqrt( g1 )

          if( g1.le.rmxnlc(la,it) )then
              numl = numl + 1
          end if

!             end if
!         end if
   end do
   end do
   end do
   numsnl(la,i_nod) = numl
#endif
#else
   do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
   do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
   do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
   m1 = mshglb(ix,iy,iz)
   if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
       if( g1.le.rmxnl2 ) then

       g1   = sqrt( g1 )

       do la = 1, lax(it)
          if( g1.le.rmxnlc(la,it) )then
              numsnl(la,i_nod) = numsnl(la,i_nod) + 1
          end if
       end do

       end if
   end if
   end do
   end do
   end do
#endif
#ifdef VECTOR
   end do
#endif
end do


!-----get the total number of mesh points : nlcmx
numwei = 0
numwe2 = 0
numwe3 = 0
do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)
   do la = 1, lax(it)
      do l = 1, lmx(it)
         if( itola(l,it) == la ) then
             ll = itolm(l,it)
             exit
         end if
      end do
      numwei = numwei + numsnl(la,i_nod)*(2*ll + 1)
      numwe2 = numwe2 + numsnl(la,i_nod)
      numwe3 = max( numwe3, numsnl(la,i_nod) )
   end do
end do

!      call gimax(numwei)
!      call gimax(numwe2)
if(loutfile(1)) write(nfile(1),*) ' *** nlcmx1&2 :', numwei, numwe2, numwe3
if(loutfile(2)) write(nfile(2),*) ' *** nlcmx1&2 :', numwei, numwe2, numwe3

!-----check array sizes
checkif: if( numwei > nlcmx1 .or. numwe2 > nlcmx2 ) then

   !-----if already allocated, deallocate arrays
   if( allocated(meshnl) ) then
       call nlcmx_dealloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1, nlcmx2 )
   end if

   !-----allocate arrays
   nlcmx1 = numwei * 1.05d0
   nlcmx2 = numwe2 * 1.05d0
   call nlcmx_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx1, nlcmx2 )

end if checkif

#ifdef VECTOR
if( numwe3 > nlcmx3 ) then
    !-----allocate arrays
    nlcmx3 = numwe3 * 1.05d0
    call nlcvector_alloc( nfile, myid, nodes,  &
& alloc_mem, nlcmx3 )
end if
#endif


!--- calculate pointer for each atom
mshsum = 0
mshsu2 = 0
do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)
   do l = 1, lmx(it)
      la = itola(l,it)
      npvnlc(l,i_nod) = mshsum
      mshsum = mshsum + numsnl(la,i_nod)
   end do
   do la = 1, lax(it)
      npmsnl(la,i_nod) = mshsu2
      mshsu2 = mshsu2 + numsnl(la,i_nod)
   end do
end do




!--- post-calculation for nonlocal term
#ifndef VECTOR
!==================== for SCALAR processors ============================
do i = 1, nionnew
do la = 1, lnref
   numsnl(la,i) = 0
end do
end do
do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)
   xx = cratm(1,i_nod)
   yy = cratm(2,i_nod)
   zz = cratm(3,i_nod)
!         nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
!         nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
!         nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
   nxx = ( hci(1,1)*xx + hci(2,1)*yy + hci(3,1)*zz )*dn1v1 + 0.5d0
   nyy = ( hci(1,2)*xx + hci(2,2)*yy + hci(3,2)*zz )*dn1v2 + 0.5d0
   nzz = ( hci(1,3)*xx + hci(2,3)*yy + hci(3,3)*zz )*dn1v3 + 0.5d0
   nxx = nxx + 1
   nyy = nyy + 1
   nzz = nzz + 1
   do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
   do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
   do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
   m1 = mshglb(ix,iy,iz)
   if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
       if( g1.le.rmxnl2 ) then

       g1xy = sqrt( g1xy )
       g1   = sqrt( g1 )
       if( g1.gt.small ) then
           cos1 = zm/g1
           sin1 = g1xy/g1
           if( g1xy.gt.small ) then
               cos2 = xm/g1xy
               sin2 = ym/g1xy
             else
               cos2 = 1.d0
               sin2 = 0.d0
           end if
         else
           cos1 = 1.d0
           sin1 = 0.d0
           cos2 = 1.d0
           sin2 = 0.d0
       end if
!            --------------------------------------
!            Spherical harmonics
       call sphrcl( lmax(it), cos1, sin1, cos2, sin2,  &
&                   ylmr, ylmi, mxl )

       do la = 1, lax(it)
          if( g1.le.rmxnlc(la,it) )then
              numl = numsnl(la,i_nod) + 1
              numsnl(la,i_nod) = numl
              meshnl(numl+npmsnl(la,i_nod)) = m1
!                    meshnl(numl+npmsnl(l,i_nod)) = 
!     &                                ix + ( (iy-1) + (iz-1)*mshy )*mshx

              m = g1/(2.d0*dltnlc(la,it))
              m = 2*m
              d = 0.5d0*( g1/dltnlc(la,it) - dble(m) )
              vnlc(la) = d*( (d-1.d0)*tabnla(m,la,it) + tabnl(m+2,la,it)  &
&                          - tabnl(m,la,it) ) + tabnl(m,la,it)
          end if
       end do

       do l = 1, lmx(it)
          la = itola(l,it)
          if( g1.le.rmxnlc(la,it) )then
              numl = numsnl(la,i_nod)
              ll  = itolm( l,it)
              lmm = itolmf(l,it)
              if( lmm == 0 ) then
                  vnlcr(numl+npvnlc(l,i_nod)) = vnlc(la)*ylmr(lmm,ll)
              else if( lmm > 0 ) then
                  vnls2  = vnlc(la)*sqrtwo
                  vnlcr(numl+npvnlc(l,i_nod)) = vnls2*ylmr(lmm,ll)
              else
                  vnls2  = vnlc(la)*sqrtwo
                  vnlcr(numl+npvnlc(l,i_nod)) = vnls2*ylmi(-lmm,ll)
              end if
          end if
       end do

       end if
   end if
   end do
   end do
   end do
end do
#else
!==================== for vector processors ============================
do i = 1, nionnew
do la = 1, lnref
   numsnl(la,i) = 0
end do
end do
do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)

   do la = 1, lax(it)

!         xx = cratm(1,i_nod) - xmn
!         yy = cratm(2,i_nod) - ymn
!         zz = cratm(3,i_nod) - zmn
!         nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
!         nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
!         nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
!---
   numl = 0
   do iz = ncornl(3,la,i_nod), nwidnl(3,la,i_nod) 
   do iy = ncornl(2,la,i_nod), nwidnl(2,la,i_nod) 
   do ix = ncornl(1,la,i_nod), nwidnl(1,la,i_nod) 
   m1 = ix + ( (iy-1) + (iz-1)*mshy )*mshx
!         m1 = mshglb(ix,iy,iz)
!         if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
!             if( g1.le.rmxnl2 ) then

           g1   = sqrt( g1 )

#ifdef SR8000
          if( g1.le.rmxnlc(la,it) )then
#endif
              numl = numl + 1
              meshnl(numl+npmsnl(la,i_nod)) = m1

              xmr(numl) = xm
              ymr(numl) = ym
              zmr(numl) = zm
              g1r(numl) = g1
#ifdef SR8000
          end if
#endif

!             end if
!         end if
   end do
   end do
   end do
   numsnl(la,i_nod) = numl


!---
   do numl = 1, numsnl(la,i_nod)
      g1 = g1r(numl)
      if( g1.le.rmxnlc(la,it) )then

       m = g1/(2.d0*dltnlc(la,it))
       m = 2*m
       d = 0.5d0*( g1/dltnlc(la,it) - dble(m) )
       vnlc = d*( (d-1.d0)*tabnla(m,la,it) + tabnl(m+2,la,it)  &
&                       - tabnl(m,la,it) ) + tabnl(m,la,it)

      vnlcrr(numl) = vnlc
      else
      vnlcrr(numl) = 0.d0
      end if
   end do

   do numl = 1, numsnl(la,i_nod)
      xm = xmr(numl)
      ym = ymr(numl)
      zm = zmr(numl)
      g1 = g1r(numl)

       g1xy = xm*xm + ym*ym
       g1xy = sqrt( g1xy )

       if( g1.gt.small ) then
           cos1 = zm/g1
           sin1 = g1xy/g1
           if( g1xy.gt.small ) then
               cos2 = xm/g1xy
               sin2 = ym/g1xy
             else
               cos2 = 1.d0
               sin2 = 0.d0
           end if
         else
           cos1 = 1.d0
           sin1 = 0.d0
           cos2 = 1.d0
           sin2 = 0.d0
       end if

      xmr(numl) = cos1
      ymr(numl) = sin1
      zmr(numl) = cos2
      g1r(numl) = sin2
   end do

   do lmxx = 1, lmx(it)
      ll  = itolm( lmxx,it)
      lmm = itolmf(lmxx,it)
   if( itola(lmxx,it) == la .and. lmm == -ll ) then
   if( ll.eq.0 ) then
       npvnl1 = npvnlc(lmxx,i_nod)
       do numl = 1, numsnl(la,i_nod)
!                c1 = xmr(numl)
!                s1 = ymr(numl)
!                c2 = zmr(numl)
!                s2 = g1r(numl)
          vnlc = vnlcrr(numl)

           vnls2  = vnlc*sqrtwo

!            --- ylmr(0,0) = sqrt( 1.d0/(4.d0*pi) )
           ylmr00 = pi00
              vnlcr(numl+npvnl1) = vnlc*ylmr00
       end do
   else if( ll.eq.1 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       do numl = 1, numsnl(la,i_nod)
          c1 = xmr(numl)
          s1 = ymr(numl)
          c2 = zmr(numl)
          s2 = g1r(numl)
          vnlc = vnlcrr(numl)

           vnls2  = vnlc*sqrtwo

!            --- ylmr(1, 0) = sqrt( 3.d0/(4.d0*pi) )*c1
           ylmr01 = pi10*c1
!            --- a1 = sqrt( 3.d0/(8.d0*pi) )*s1
           a1 = pi11*s1
           ylmr11 = a1*c2
           ylmi11 = a1*s2
              qqa(numl) = vnlc*ylmr01
              qqb(numl) = vnls2*ylmr11
              qqc(numl) = vnls2*ylmi11
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1) = qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2) = qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3) = qqc(numl)
       end do
   else if( ll.eq.2 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       npvnl4 = npvnlc(lmxx + 3,i_nod)
       npvnl5 = npvnlc(lmxx + 4,i_nod)
       do numl = 1, numsnl(la,i_nod)
          c1 = xmr(numl)
          s1 = ymr(numl)
          c2 = zmr(numl)
          s2 = g1r(numl)
          vnlc = vnlcrr(numl)

           vnls2  = vnlc*sqrtwo

           c12 = c1*c1
           s12 = s1*s1
           c22p = c2*c2 - s2*s2
           s22p = s2*c2 + c2*s2
!            --- ylmr(0, 2) = sqrt( 5.d0/(16.d0*pi) )*(3.d0*c12 - 1.d0)
           ylmr02 = pi20*(3.d0*c12 - 1.d0)
!            --- a1 = sqrt( 15.d0/(8.d0*pi) )*s1*c1
           a1 = pi21*s1*c1
           ylmr12 = a1*c2
           ylmi12 = a1*s2
!            --- a1 = sqrt( 15.d0/(32.d0*pi) )*s12
           a1 = pi22*s12
           ylmr22 = a1*c22p
           ylmi22 = a1*s22p
              qqa(numl) = vnlc*ylmr02
              qqb(numl) = vnls2*ylmr12
              qqc(numl) = vnls2*ylmi12
              qqd(numl) = vnls2*ylmr22
              qqe(numl) = vnls2*ylmi22
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1) = qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2) = qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl4) = qqc(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3) = qqd(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl5) = qqe(numl)
       end do
   else if( ll.eq.3 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       npvnl4 = npvnlc(lmxx + 3,i_nod)
       npvnl5 = npvnlc(lmxx + 4,i_nod)
       npvnl6 = npvnlc(lmxx + 5,i_nod)
       npvnl7 = npvnlc(lmxx + 6,i_nod)
       do numl = 1, numsnl(la,i_nod)
          c1 = xmr(numl)
          s1 = ymr(numl)
          c2 = zmr(numl)
          s2 = g1r(numl)
          vnlc = vnlcrr(numl)

           vnls2  = vnlc*sqrtwo

           c12 = c1*c1
           s12 = s1*s1
           c22p = c2*c2 - s2*s2
           s22p = s2*c2 + c2*s2
           c13 = c12*c1
           s13 = s12*s1
           c23p = c22p*c2 - s22p*s2
           s23p = s22p*c2 + c22p*s2
!            --- ylmr(0, 3) = sqrt( 7.d0/(16.d0*pi) )*c1*(5.d0*c12 - 3.d0)
           ylmr03 = pi30*c1*(5.d0*c12 - 3.d0)
!            --- a1 = sqrt( 21.d0/(64.d0*pi) )*s1*(5.d0*c12 - 1.d0)
           a1 = pi31*s1*(5.d0*c12 - 1.d0)
           ylmr13 = a1*c2
           ylmi13 = a1*s2
!            --- a1 = sqrt( 105.d0/(32.d0*pi) )*s12*c1
           a1 = pi32*s12*c1
           ylmr23 = a1*c22p
           ylmi23 = a1*s22p
!            --- a1 = sqrt( 35.d0/(64.d0*pi) )*s13
           a1 = pi33*s13
           ylmr33 = a1*c23p
           ylmi33 = a1*s23p
              qqa(numl) = vnlc*ylmr03
              qqb(numl) = vnls2*ylmr13
              qqc(numl) = vnls2*ylmi13
              qqd(numl) = vnls2*ylmr23
              qqe(numl) = vnls2*ylmi23
              qqf(numl) = vnls2*ylmr33
              qqg(numl) = vnls2*ylmi33
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1) = qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2) = qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl5) = qqc(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3) = qqd(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl6) = qqe(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl4) = qqf(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl7) = qqg(numl)
       end do
   end if

   end if
   end do
   end do

end do
!==================== end ==============================================
#endif

!check
!      do i_nod = 1, nionnew
!         i = iatmpt(i_nod)
!         it = iatoit(i)
!         lmxx = 0
!         l = 0
!c         do l = 0, lmax(it)
!c         if( lchk(l,it) .and. l.ne.lclno(it) ) then
!         tsum = 0.d0
!         j   = lmxx + 1
!         do numl = 1, numsnl(l,i_nod)
!            tsum = tsum + vnlcr(numl+npvnlc(j,i_nod))**2
!         end do
!         if( myid.eq.0 ) write(*,*) i_nod, l, j, tsum
!c         lmxx = lmxx + 2*l + 1
!c         end if
!c         end do
!      end do
!c-----Finalize the parallel environment
!      call end_parallel(ierr)
!
!      stop

return
end




subroutine calnlc( nfile, myid, nodes,  &
& ntype, lmax, mxl, lchk, lclno, rk, eigv, nhk1, nhk2,  &
& lkbppi, lking, iatoit, rdelv, nylmmx_, nd1vks  )
!-----------------------------------------------------------------------
!     calculation of nonlocal pseudopotentials
!-----------------------------------------------------------------------
#ifdef VECTOR
use nlcvector
#endif
use nlpp_variables
use nlppr_variables
use nlkbppr_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ntype, nhk1(*), nhk2(*), iatoit(*)
integer :: lmax(*), mxl, lclno(*)
logical :: lchk(0:mxl,*)
real*8  :: rk(*), eigv(0:*)
logical :: lkbppi(ntype), lking(ntype)
real*8  :: rdelv
integer :: nylmmx_, nd1vks(3)

!------declare local variables
integer :: kfft1, kfft2, kfft3, i, j, ia, it, l, la, mm, m1
integer :: ix, iy, iz


#ifdef VECTOR
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)
#endif

sumr(1:nylmmx,1:nhk2(ntype)) = 0.d0

do i = 1, nionnew
   ia = iatmpt(i)
   it = iatoit(ia)
#ifndef VECTOR
   do l = 1, lmx(it)
      la = itola(l,it)
      do mm = 1, numsnl(la,i)
         m1 = meshnl(mm+npmsnl(la,i))
         sumr(l,ia) = sumr(l,ia) + vnlcr(mm+npvnlc(l,i))*eigv(m1)
      end do
!cc            sumr(j,ia) = sumr(j,ia) * rdelv * clalp_ncr(la,it)
   end do
#else
   do la = 1, lax(it)
      do mm = 1, numsnl(la,i)
         m1 = meshnl(mm+npmsnl(la,i))
         xxnl(mm) = eigv(m1)
      end do
   do l = 1, lmx(it)
      if( itola(l,it) == la ) then
          do mm = 1, numsnl(la,i)
             sumr(l,ia) = sumr(l,ia) + vnlcr(mm+npvnlc(l,i))*xxnl(mm)
          end do
!cc            sumr(j,ia) = sumr(j,ia) * rdelv * clalp_ncr(la,it)
      end if
   end do
   end do
#endif
end do


do it = 1, ntype
if( lkbppi(it) .and. lking(it) ) then
    do ia = nhk1(it), nhk2(it)
       do l = 1, lmx(it)
          la = itola(l,it)
          sumr(l,ia) = sumr(l,ia) * rdelv * clalp_ncr(la,it)
       end do
    end do
end if
end do


!cc      call bacopy( sumr, lmxxxx, t_comm )


do i = 1, nionnew
   ia = iatmpt(i)
   it = iatoit(ia)
#ifndef VECTOR
   do l = 1, lmx(it)
      la = itola(l,it)
      do mm = 1, numsnl(la,i)
         m1 = meshnl(mm+npmsnl(la,i))
         rk(m1) = rk(m1) + vnlcr(mm+npvnlc(l,i))*sumr(l,ia)
      end do
   end do
#else
   do la = 1, lax(it)
      do mm = 1, numsnl(la,i)
         xxnl(mm) = 0.d0
      end do
   do l = 1, lmx(it)
      if( itola(l,it) == la ) then
          do mm = 1, numsnl(la,i)
             xxnl(mm) = xxnl(mm) + vnlcr(mm+npvnlc(l,i))*sumr(l,ia)
          end do
      end if
   end do
#ifdef SR8000
      do mm = 1, numsnl(la,i)
         m1 = meshnl(mm+npmsnl(la,i))
         rk(m1)= rk(m1) + xxnl(mm)
      end do
#else
   mm = 0
   do iz = ncornl(3,la,i), nwidnl(3,la,i) 
   do iy = ncornl(2,la,i), nwidnl(2,la,i) 
   do ix = ncornl(1,la,i), nwidnl(1,la,i) 
      m1 = ix + ((iy-1) + (iz-1)*kfft2 )*kfft1
      mm = mm + 1
      rk(m1)= rk(m1) + xxnl(mm)
   end do
   end do
   end do
#endif
   end do
#endif
end do


return
end




subroutine savesumr( nfile, myid, nodes,  &
& i, ntype, nhk1, nhk2, natom, lkbppi, lking, nylmmx_, nspin )
!-----------------------------------------------------------------------
!     store sumr
!-----------------------------------------------------------------------
use nlpp_variables
use nlkbppr_variables
use nlkbpp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: i
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lkbppi, lking
integer :: nylmmx_, nspin
!------declare local variables
integer :: it, nk


call get_kvec( nk )

do it = 1, ntype
if( lkbppi(it) .and. lking(it) ) then
    svsumr(1:nylmmx,nhk1(it):nhk2(it),i,nspin,nk) = sumr(1:nylmmx,nhk1(it):nhk2(it))
    if( .not.lgamma .or. lnoncollinear ) then
        svsumi(1:nylmmx,nhk1(it):nhk2(it),i,nspin,nk) = sumi(1:nylmmx,nhk1(it):nhk2(it))
    end if

    if( lnoncollinear ) then
        svsubr(1:nylmmx,nhk1(it):nhk2(it),i,nspin,nk) = subr(1:nylmmx,nhk1(it):nhk2(it))
        svsubi(1:nylmmx,nhk1(it):nhk2(it),i,nspin,nk) = subi(1:nylmmx,nhk1(it):nhk2(it))
    end if
end if
end do


return
end




subroutine fnnlcl( nfile, myid, nodes,  &
& fnlc, bufst, ltimecnt, lcstress, t_comm,  &
& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
& lkbppi, lking, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& mshglb, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& nspin, dbuf, dbufr, natom, mxl, eigr, eigi )
!-----------------------------------------------------------------------
!   fnlc : force by nonlocal pseudopotential
!
!   This routine fast but needs big memory.
!-----------------------------------------------------------------------
use ncmagne_variables
use nlpp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
real*8  :: fnlc(3,natom)
real(8) :: bufst(6)
logical :: ltimecnt, lcstress
real(8) :: t_comm
integer :: ntype, nhk1(ntype), nhk2(ntype)
integer :: lmax(ntype), lclno(ntype), iatoit(*)
integer :: mxl
logical :: lchk(0:mxl,ntype)
logical :: lkbppi(ntype), lking(ntype)
real(8) :: rdel(3), rdelg(3,3), rdelv, hcell(3,3), hci(3,3), rvol
logical :: lorthrhmbc
integer :: kfft1, kfft2, kfft3
integer :: mshglb(kfft1,kfft2,kfft3)
integer :: mshgnx(*), mshgny(*), mshgnz(*)
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
real*8  :: occ(nband)
real(8) :: cgjr(nplwex,*)
integer :: iod(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real(8) :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh, ijkgd(-ngenh:ngenh)
integer :: nspin
real*8  :: dbuf(3*natom), dbufr(3*natom)
real*8  :: eigr(0:*), eigi(0:*)


!if( .not.lnoncollinear ) then
    call fnnlcl2( nfile, myid, nodes,  &
& fnlc, bufst, ltimecnt, lcstress, t_comm,  &
& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
& lkbppi, lking, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& mshglb, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& nspin, dbuf, dbufr, natom, mxl )
!else
!    !-----noncollinear magnetism
!    call ncfnnlcl2( nfile, myid, nodes,  &
!& fnlc, bufst, ltimecnt, lcstress, t_comm,  &
!& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
!& lkbppi, lking, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
!& mshglb, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
!& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
!& fft3x, fft3y, fftwork, ncijkg,  &
!& nspin, dbuf, dbufr, natom, mxl, eigr, eigi )
!end if


return
end




subroutine fnnlcl2( nfile, myid, nodes,  &
& fnlc, bufst, ltimecnt, lcstress, t_comm,  &
& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
& lkbppi, lking, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& mshglb, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& nspin, dbuf, dbufr, natom, mxl )
!-----------------------------------------------------------------------
!   fnlc : force by nonlocal pseudopotential
!-----------------------------------------------------------------------
use outfile
#ifdef VECTOR
use nlcvector
#endif
use nlpp_variables
use nlppr_variables
use nlkbppr_variables
use nlkbpp_variables
use force_nlc_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
real*8,  dimension(3,natom) :: fnlc
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lkbppi, lking
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lclno
logical, dimension(0:mxl,ntype) :: lchk
real*8,  dimension(nband) :: occ
real*8,  dimension(3*natom) :: dbuf, dbufr
real*8,  dimension(3,3) :: rdelg
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3,3) :: hci
logical :: lorthrhmbc

dimension iatoit(*), rdel(3)
dimension mshglb(kfft1,kfft2,kfft3)
dimension mshgnx(*), mshgny(*), mshgnz(*)

dimension ndv3nl(3)
dimension bufst(6)
logical   ltimecnt, lcstress

dimension cgjr(nplwex,*)
dimension iod(*)
dimension mfd2ft(*), nd1vks(*)
dimension fft3x(*), fft3y(*)
complex*16 fftwork(*)
dimension ijkgd(-ngenh:ngenh)
!dimension gdcrsv(*), wegud(*)
!logical   lspin, lcgjsv


call get_kvec( nk )

mshx1 = 1
mshy1 = 1
mshz1 = 1
mshx  = kfft1
mshy  = kfft2
mshz  = kfft3
if( ltimecnt ) then
    ct0 = timecnt()
    t_com0 = t_comm
end if

!do i = 1, 6
!   bufst(i) = 0.d0
!end do


rmxnl2 = rmxnlx*rmxnlx
do i = 1, 3
!   ndv3nl(i) = xc(i)*hcell(i,i)/rdel(i) + 1.d0
   ndv3nl(i) = xc(i)*vecratio(hcell(1,i),rdelg(1,i)) + 1.d0
end do

small=1.0d-08
sqrtwo = sqrt(2.d0)

!dn1v1 = hcell(1,1)/rdel(1)
!dn1v2 = hcell(2,2)/rdel(2)
!dn1v3 = hcell(3,3)/rdel(3)
dn1v1 = vecratio(hcell(1,1),rdelg(1,1))
dn1v2 = vecratio(hcell(1,2),rdelg(1,2))
dn1v3 = vecratio(hcell(1,3),rdelg(1,3))
do igxxx = 1, 3
!=======================================================================
!   if igxxx = 1, force x component is calculated
!   if igxxx = 2, force y component is calculated
!   if igxxx = 3, force z component is calculated
!=======================================================================
!      ct0 = timecnt()
#ifndef VECTOR
!==================== for SCALAR processors ============================
do i = 1, nionnew
do la = 1, lnref
   numsnl(la,i) = 0
end do
end do

do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)
   xx = cratm(1,i_nod)
   yy = cratm(2,i_nod)
   zz = cratm(3,i_nod)
!         nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
!         nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
!         nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
   nxx = ( hci(1,1)*xx + hci(2,1)*yy + hci(3,1)*zz )*dn1v1 + 0.5d0
   nyy = ( hci(1,2)*xx + hci(2,2)*yy + hci(3,2)*zz )*dn1v2 + 0.5d0
   nzz = ( hci(1,3)*xx + hci(2,3)*yy + hci(3,3)*zz )*dn1v3 + 0.5d0
   nxx = nxx + 1
   nyy = nyy + 1
   nzz = nzz + 1
   do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
   do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
   do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
   m1 = mshglb(ix,iy,iz)
   if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
       if( g1.le.rmxnl2 ) then

       g1xy = sqrt( g1xy )
       g1   = sqrt( g1 )
       if( g1.gt.small ) then
           cos1 = zm/g1
           sin1 = g1xy/g1
           if( g1xy.gt.small ) then
               cos2 = xm/g1xy
               sin2 = ym/g1xy
             else
               cos2 = 1.d0
               sin2 = 0.d0
           end if
         else
           cos1 = 1.d0
           sin1 = 0.d0
           cos2 = 1.d0
           sin2 = 0.d0
       end if
!            --------------------------------------
!            Spherical harmonics
       call sphrb2( lmax(it), cos1, sin1, cos2, sin2,  &
&                   ylmr, ylmi, mxl )
!            --- unit vector ---
       erx = sin1*cos2
       ery = sin1*sin2
       erz = cos1
       esx = cos1*cos2
       esy = cos1*sin2
       esz = -sin1
       efx = -sin2
       efy =  cos2
       efz = 0.d0

       do la = 1, lax(it)
          if( g1.le.rmxnlc(la,it) )then
              numl = numsnl(la,i_nod) + 1
              numsnl(la,i_nod) = numl
              meshnl(numl+npmsnl(la,i_nod)) = m1
              m = g1/(2.d0*dltnlc(la,it))
              m = 2*m
              d = 0.5d0*( g1/dltnlc(la,it) - dble(m) )
              vnlc(la) = d*( (d-1.d0)*tabnla(m,la,it) + tabnl(m+2,la,it)  &
&                          - tabnl(m,la,it) ) + tabnl(m,la,it)
              fvnlc(la)= d*( (d-1.d0)*tbfnla(m,la,it) + tbfnl(m+2,la,it)  &
&                          - tbfnl(m,la,it) ) + tbfnl(m,la,it)
          end if
       end do

       do l = 1, lmx(it)
          la = itola(l,it)
          if( g1.le.rmxnlc(la,it) )then
              numl = numsnl(la,i_nod)
              numlj = numl+npvnlc(l,i_nod)
              ll  = itolm( l,it)
              lmm = itolmf(l,it)
              if( lmm == 0 ) then
!cc                   vnlcr(numlj) = vnlc*ylmr(1,lmm,l)
                  conrr = fvnlc(la)*ylmr(1,lmm,ll)
                  if( g1.gt.small ) then
                      vnlcg1 = vnlc(la)/g1
                  else
                      vnlcg1 = fvnlc(la)
                  end if
                  consr = vnlcg1*ylmr(2,lmm,ll)
                  confr = vnlcg1*ylmr(3,lmm,ll)
                  if( igxxx.eq.1 ) then
                      vnlcr(numlj)=conrr*erx + consr*esx + confr*efx
                  else if( igxxx.eq.2 ) then
                      vnlcr(numlj)=conrr*ery + consr*esy + confr*efy
                  else
                      vnlcr(numlj)=conrr*erz + consr*esz + confr*efz
                  end if
              else
                   vnls2 =  vnlc(la)*sqrtwo
                  fvnls2 = fvnlc(la)*sqrtwo
                  if( g1.gt.small ) then
                      vnlcg1 = vnls2/g1
                  else
                      vnlcg1 = fvnls2
                  end if
                  if( lmm > 0 ) then
!cc                   vnlcr(numlj) = vnls2*ylmr(1,lmm,l)
                      conrr = fvnls2*ylmr(1,lmm,ll)
                      consr = vnlcg1*ylmr(2,lmm,ll)
                      confr = vnlcg1*ylmr(3,lmm,ll)
                      if( igxxx.eq.1 ) then
                          vnlcr(numlj)=conrr*erx + consr*esx + confr*efx
                      else if( igxxx.eq.2 ) then
                          vnlcr(numlj)=conrr*ery + consr*esy + confr*efy
                      else
                          vnlcr(numlj)=conrr*erz + consr*esz + confr*efz
                      end if
                  else
!cc                   vnlcr(numlj)=vnls2*ylmi(1,lmm,l)
                      conri = fvnls2*ylmi(1,-lmm,ll)
                      consi = vnlcg1*ylmi(2,-lmm,ll)
                      confi = vnlcg1*ylmi(3,-lmm,ll)
                      if( igxxx.eq.1 ) then
                          vnlcr(numlj)=conri*erx + consi*esx + confi*efx
                      else if( igxxx.eq.2 ) then
                          vnlcr(numlj)=conri*ery + consi*esy + confi*efy
                      else
                          vnlcr(numlj)=conri*erz + consi*esz + confi*efz
                      end if
                  end if
              end if
          end if
       end do
       end if
   end if
   end do
   end do
   end do
end do
#else
!==================== for vector processors ============================
do i = 1, nionnew
do la = 1, lnref
   numsnl(la,i) = 0
end do
end do

do i_nod = 1, nionnew
   i = iatmpt(i_nod)
   it = iatoit(i)

   do la = 1, lax(it)

!         xx = cratm(1,i_nod)
!         yy = cratm(2,i_nod)
!         zz = cratm(3,i_nod)
!         nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
!         nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
!         nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
   numl = 0
   do iz = ncornl(3,la,i_nod), nwidnl(3,la,i_nod) 
   do iy = ncornl(2,la,i_nod), nwidnl(2,la,i_nod) 
   do ix = ncornl(1,la,i_nod), nwidnl(1,la,i_nod) 
!         m1 = ix + ( (iy-1) + (iz-1)*mshy )*mshx
   m1 = mshglb(ix,iy,iz)
   if( m1.gt.0 ) then
!             ixx = mshx1 + ix - 1
!             xm = dble(ixx-1)*rdel(1)
!             iyy = mshy1 + iy - 1
!             ym = dble(iyy-1)*rdel(2)
!             izz = mshz1 + iz - 1
!             zm = dble(izz-1)*rdel(3)
       d1 = ix - 1
       d2 = iy - 1
       d3 = iz - 1
       xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
       ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
       zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

       xm = xm - cratm(1,i_nod)
       ym = ym - cratm(2,i_nod)
       zm = zm - cratm(3,i_nod)
       g1xy = xm*xm + ym*ym
       g1   = g1xy + zm*zm
!             if( g1.le.rmxnl2 ) then

           g1   = sqrt( g1 )

          if( g1.le.rmxnlc(la,it) )then
              numl = numl + 1
              meshnl(numl+npmsnl(la,i_nod)) = m1

              xmr(numl) = xm
              ymr(numl) = ym
              zmr(numl) = zm
              g1r(numl) = g1
          end if

!             end if
   end if
   end do
   end do
   end do
   numsnl(la,i_nod) = numl


   do numl = 1, numsnl(la,i_nod)
      g1 = g1r(numl)
!            if( g1.le.rmxnlc(la,it) )then

       m = g1/(2.d0*dltnlc(la,it))
       m = 2*m
       d = 0.5d0*( g1/dltnlc(la,it) - dble(m) )
       vnlc = d*( (d-1.d0)*tabnla(m,la,it) + tabnl(m+2,la,it)  &
&                       - tabnl(m,la,it) ) + tabnl(m,la,it)
       fvnlc= d*( (d-1.d0)*tbfnla(m,la,it) + tbfnl(m+2,la,it)  &
&                       - tbfnl(m,la,it) ) + tbfnl(m,la,it)

       if( g1.gt.small ) then
           vnlcg1 = vnlc/g1
         else
           vnlcg1 = fvnlc
       end if
      vnlcrr(numl) = vnlcg1
      xxnl(numl)   = fvnlc
!            else
!            vnlcrr(numl) = 0.d0
!            xxnl(numl)   = 0.d0
!            end if
   end do

   do numl = 1, numsnl(la,i_nod)
      xm = xmr(numl)
      ym = ymr(numl)
      zm = zmr(numl)
      g1 = g1r(numl)

       g1xy = xm*xm + ym*ym
       g1xy = sqrt( g1xy )

       if( g1.gt.small ) then
           cos1 = zm/g1
           sin1 = g1xy/g1
           if( g1xy.gt.small ) then
               cos2 = xm/g1xy
               sin2 = ym/g1xy
             else
               cos2 = 1.d0
               sin2 = 0.d0
           end if
         else
           cos1 = 1.d0
           sin1 = 0.d0
           cos2 = 1.d0
           sin2 = 0.d0
       end if

      xmr(numl) = cos1
      ymr(numl) = sin1
      zmr(numl) = cos2
      g1r(numl) = sin2
   end do


   do lmxx = 1, lmx(it)
      ll  = itolm( lmxx,it)
      lmm = itolmf(lmxx,it)
   if( itola(lmxx,it) == la .and. lmm == -ll ) then
   if( ll.eq.0 ) then
       npvnl1 = npvnlc(lmxx,i_nod)
   if( igxxx.eq.1 ) then
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
          sin1 = ymr(numl)
          cos2 = zmr(numl)
          sin2 = g1r(numl)
!                vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
          erx = sin1*cos2
!                ery = sin1*sin2
!                erz = cos1

!           --- ylmr(1,0,0) = sqrt( 1.d0/(4.d0*pi) )
          ylmr100 = pi00
             conrr =  fvnlc*ylmr100
                 vnlcr(numl+npvnl1)=conrr*erx
       end do
   else if( igxxx.eq.2 ) then
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
          sin1 = ymr(numl)
          cos2 = zmr(numl)
          sin2 = g1r(numl)
!                vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
!                erx = sin1*cos2
          ery = sin1*sin2
!                erz = cos1

!           --- ylmr(1,0,0) = sqrt( 1.d0/(4.d0*pi) )
          ylmr100 = pi00
             conrr =  fvnlc*ylmr100
                 vnlcr(numl+npvnl1)=conrr*ery
       end do
   else
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
!                sin1 = ymr(numl)
!                cos2 = zmr(numl)
!                sin2 = g1r(numl)
!                vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
!                erx = sin1*cos2
!                ery = sin1*sin2
          erz = cos1

!           --- ylmr(1,0,0) = sqrt( 1.d0/(4.d0*pi) )
          ylmr100 = pi00
             conrr =  fvnlc*ylmr100
                 vnlcr(numl+npvnl1)=conrr*erz
       end do
   end if

   else if( ll.eq.1 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
          sin1 = ymr(numl)
          cos2 = zmr(numl)
          sin2 = g1r(numl)
          vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
          erx = sin1*cos2
          ery = sin1*sin2
          erz = cos1
          esx = cos1*cos2
          esy = cos1*sin2
          esz = -sin1
          efx = -sin2
          efy =  cos2
          efz = 0.d0

          c1 = cos1
          s1 = sin1
          c2 = cos2
          s2 = sin2

!           --- ylmr(1,1, 0) = sqrt( 3.d0/(4.d0*pi) )*c1
          ylmr101 = pi10*c1
          ylmr201 = -pi10*s1
          ylmr301 = 0.d0
!           --- a1 = sqrt( 3.d0/(8.d0*pi) )*s1
          a1  = pi11*s1
          ac1 = pi11*c1
          ylmr111 = a1*c2
          ylmi111 = a1*s2
          ylmr211 = ac1*c2
          ylmi211 = ac1*s2
          ylmr311 = -pi11*s2
          ylmi311 =  pi11*c2

             numlj = numl+npvnl1
             conrr =  fvnlc*ylmr101
             consr = vnlcg1*ylmr201
             confr = vnlcg1*ylmr301

          fvnls2 = fvnlc*sqrtwo
          vnlcg2 = vnlcg1*sqrtwo
             c1nrr = fvnls2*ylmr111
             c1nri = fvnls2*ylmi111
             c1nsr = vnlcg2*ylmr211
             c1nsi = vnlcg2*ylmi211
             c1nfr = vnlcg2*ylmr311
             c1nfi = vnlcg2*ylmi311

             if( igxxx.eq.1 ) then
              qqa(numl)=conrr*erx + consr*esx + confr*efx
              qqb(numl)=c1nrr*erx + c1nsr*esx + c1nfr*efx
              qqc(numl)=c1nri*erx + c1nsi*esx + c1nfi*efx
             else if( igxxx.eq.2 ) then
              qqa(numl)=conrr*ery + consr*esy + confr*efy
              qqb(numl)=c1nrr*ery + c1nsr*esy + c1nfr*efy
              qqc(numl)=c1nri*ery + c1nsi*esy + c1nfi*efy
             else
              qqa(numl)=conrr*erz + consr*esz + confr*efz
              qqb(numl)=c1nrr*erz + c1nsr*esz + c1nfr*efz
              qqc(numl)=c1nri*erz + c1nsi*esz + c1nfi*efz
             end if
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1)=qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2)=qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3)=qqc(numl)
       end do

   else if( ll.eq.2 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       npvnl4 = npvnlc(lmxx + 3,i_nod)
       npvnl5 = npvnlc(lmxx + 4,i_nod)
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
          sin1 = ymr(numl)
          cos2 = zmr(numl)
          sin2 = g1r(numl)
          vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
          erx = sin1*cos2
          ery = sin1*sin2
          erz = cos1
          esx = cos1*cos2
          esy = cos1*sin2
          esz = -sin1
          efx = -sin2
          efy =  cos2
          efz = 0.d0

          c1 = cos1
          s1 = sin1
          c2 = cos2
          s2 = sin2

          c12 = c1*c1
          s12 = s1*s1
          c22p = c2*c2 - s2*s2
          s22p = s2*c2 + c2*s2
!           --- ylmr(1, 0, 2) = sqrt( 5.d0/(16.d0*pi) )*(3.d0*c12 - 1.d0)
          ylmr102 = pi20*(3.d0*c12 - 1.d0)
          ylmr202 = -pi20*6.d0*c1*s1
          ylmr302 = 0.d0
!           --- a1 = sqrt( 15.d0/(8.d0*pi) )*s1*c1
          a01 = pi21*c1
          a1  =  a01*s1
          ac1 = pi21*( c12 - s12 )
          ylmr112 = a1*c2
          ylmi112 = a1*s2
          ylmr212 = ac1*c2
          ylmi212 = ac1*s2
          ylmr312 = -a01*s2
          ylmi312 =  a01*c2
!           --- a1 = sqrt( 15.d0/(32.d0*pi) )*s12
          a1  = pi22*s12
          a01 = 2.d0*pi22*s1
          ac1 =  a01*c1
          ylmr122 = a1*c22p
          ylmi122 = a1*s22p
          ylmr222 = ac1*c22p
          ylmi222 = ac1*s22p
          ylmr322 = -a01*s22p
          ylmi322 =  a01*c22p

             conrr =  fvnlc*ylmr102
             consr = vnlcg1*ylmr202
             confr = vnlcg1*ylmr302

          fvnls2 = fvnlc*sqrtwo
          vnlcg2 = vnlcg1*sqrtwo
             c1nrr = fvnls2*ylmr112
             c1nri = fvnls2*ylmi112
             c1nsr = vnlcg2*ylmr212
             c1nsi = vnlcg2*ylmi212
             c1nfr = vnlcg2*ylmr312
             c1nfi = vnlcg2*ylmi312
             c2nrr = fvnls2*ylmr122
             c2nri = fvnls2*ylmi122
             c2nsr = vnlcg2*ylmr222
             c2nsi = vnlcg2*ylmi222
             c2nfr = vnlcg2*ylmr322
             c2nfi = vnlcg2*ylmi322

             if( igxxx.eq.1 ) then
              qqa(numl)=conrr*erx + consr*esx + confr*efx
              qqb(numl)=c1nrr*erx + c1nsr*esx + c1nfr*efx
              qqc(numl)=c1nri*erx + c1nsi*esx + c1nfi*efx
              qqd(numl)=c2nrr*erx + c2nsr*esx + c2nfr*efx
              qqe(numl)=c2nri*erx + c2nsi*esx + c2nfi*efx
             else if( igxxx.eq.2 ) then
              qqa(numl)=conrr*ery + consr*esy + confr*efy
              qqb(numl)=c1nrr*ery + c1nsr*esy + c1nfr*efy
              qqc(numl)=c1nri*ery + c1nsi*esy + c1nfi*efy
              qqd(numl)=c2nrr*ery + c2nsr*esy + c2nfr*efy
              qqe(numl)=c2nri*ery + c2nsi*esy + c2nfi*efy
             else
              qqa(numl)=conrr*erz + consr*esz + confr*efz
              qqb(numl)=c1nrr*erz + c1nsr*esz + c1nfr*efz
              qqc(numl)=c1nri*erz + c1nsi*esz + c1nfi*efz
              qqd(numl)=c2nrr*erz + c2nsr*esz + c2nfr*efz
              qqe(numl)=c2nri*erz + c2nsi*esz + c2nfi*efz
             end if
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1)=qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2)=qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl4)=qqc(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3)=qqd(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl5)=qqe(numl)
       end do

   else if( ll.eq.3 ) then
       npvnl1 = npvnlc(lmxx + 0,i_nod)
       npvnl2 = npvnlc(lmxx + 1,i_nod)
       npvnl3 = npvnlc(lmxx + 2,i_nod)
       npvnl4 = npvnlc(lmxx + 3,i_nod)
       npvnl5 = npvnlc(lmxx + 4,i_nod)
       npvnl6 = npvnlc(lmxx + 5,i_nod)
       npvnl7 = npvnlc(lmxx + 6,i_nod)
       do numl = 1, numsnl(la,i_nod)
          cos1 = xmr(numl)
          sin1 = ymr(numl)
          cos2 = zmr(numl)
          sin2 = g1r(numl)
          vnlcg1 = vnlcrr(numl)
          fvnlc  = xxnl(numl)

!            --- unit vector ---
          erx = sin1*cos2
          ery = sin1*sin2
          erz = cos1
          esx = cos1*cos2
          esy = cos1*sin2
          esz = -sin1
          efx = -sin2
          efy =  cos2
          efz = 0.d0

          c1 = cos1
          s1 = sin1
          c2 = cos2
          s2 = sin2

          c12 = c1*c1
          s12 = s1*s1
          c22p = c2*c2 - s2*s2
          s22p = s2*c2 + c2*s2
          c13 = c12*c1
          s13 = s12*s1
          c23p = c22p*c2 - s22p*s2
          s23p = s22p*c2 + c22p*s2
!           --- ylmr(1,0, 3) = sqrt( 7.d0/(16.d0*pi) )*c1*(5.d0*c12 - 3.d0)
          ylmr103 = pi30*c1*( 5.d0*c12 - 3.d0)
          ylmr203 =-pi30*s1*(15.d0*c12 - 3.d0)
          ylmr303 = 0d0
!           --- a1 = sqrt( 21.d0/(64.d0*pi) )*s1*(5.d0*c12 - 1.d0)
          a01 = pi31*(5.d0*c12 - 1.d0)
          a1  = a01*s1
          ac1 = pi31*c1*(5.d0*c12 - 10d0*s12 - 1.d0)
          ylmr113 = a1*c2
          ylmi113 = a1*s2
          ylmr213 = ac1*c2
          ylmi213 = ac1*s2
          ylmr313 = -a01*s2
          ylmi313 =  a01*c2
!           --- a1 = sqrt( 105.d0/(32.d0*pi) )*s12*c1
          a1  = pi32*s12*c1
          ac1 = pi32*(2d0*c12*s1 - s13)
          a01 = 2d0*pi32*s1*c1
          ylmr123 = a1*c22p
          ylmi123 = a1*s22p
          ylmr223 = ac1*c22p
          ylmi223 = ac1*s22p
          ylmr323 = -a01*s22p
          ylmi323 =  a01*c22p
!           --- a1 = sqrt( 35.d0/(64.d0*pi) )*s13
          a1  = pi33*s13
          a01 = 3d0*pi33*s12
          ac1 = a01*c1
          ylmr133 = a1*c23p
          ylmi133 = a1*s23p
          ylmr233 = ac1*c23p
          ylmi233 = ac1*s23p
          ylmr333 = -a01*s23p
          ylmi333 =  a01*c23p

             conrr =  fvnlc*ylmr103
             consr = vnlcg1*ylmr203
             confr = vnlcg1*ylmr303

          fvnls2 = fvnlc*sqrtwo
          vnlcg2 = vnlcg1*sqrtwo
             c1nrr = fvnls2*ylmr113
             c1nri = fvnls2*ylmi113
             c1nsr = vnlcg2*ylmr213
             c1nsi = vnlcg2*ylmi213
             c1nfr = vnlcg2*ylmr313
             c1nfi = vnlcg2*ylmi313
             c2nrr = fvnls2*ylmr123
             c2nri = fvnls2*ylmi123
             c2nsr = vnlcg2*ylmr223
             c2nsi = vnlcg2*ylmi223
             c2nfr = vnlcg2*ylmr323
             c2nfi = vnlcg2*ylmi323
             c3nrr = fvnls2*ylmr133
             c3nri = fvnls2*ylmi133
             c3nsr = vnlcg2*ylmr233
             c3nsi = vnlcg2*ylmi233
             c3nfr = vnlcg2*ylmr333
             c3nfi = vnlcg2*ylmi333

             if( igxxx.eq.1 ) then
              qqa(numl)=conrr*erx + consr*esx + confr*efx
              qqb(numl)=c1nrr*erx + c1nsr*esx + c1nfr*efx
              qqc(numl)=c1nri*erx + c1nsi*esx + c1nfi*efx
              qqd(numl)=c2nrr*erx + c2nsr*esx + c2nfr*efx
              qqe(numl)=c2nri*erx + c2nsi*esx + c2nfi*efx
              qqf(numl)=c3nrr*erx + c3nsr*esx + c3nfr*efx
              qqg(numl)=c3nri*erx + c3nsi*esx + c3nfi*efx
             else if( igxxx.eq.2 ) then
              qqa(numl)=conrr*ery + consr*esy + confr*efy
              qqb(numl)=c1nrr*ery + c1nsr*esy + c1nfr*efy
              qqc(numl)=c1nri*ery + c1nsi*esy + c1nfi*efy
              qqd(numl)=c2nrr*ery + c2nsr*esy + c2nfr*efy
              qqe(numl)=c2nri*ery + c2nsi*esy + c2nfi*efy
              qqf(numl)=c3nrr*ery + c3nsr*esy + c3nfr*efy
              qqg(numl)=c3nri*ery + c3nsi*esy + c3nfi*efy
             else
              qqa(numl)=conrr*erz + consr*esz + confr*efz
              qqb(numl)=c1nrr*erz + c1nsr*esz + c1nfr*efz
              qqc(numl)=c1nri*erz + c1nsi*esz + c1nfi*efz
              qqd(numl)=c2nrr*erz + c2nsr*esz + c2nfr*efz
              qqe(numl)=c2nri*erz + c2nsi*esz + c2nfi*efz
              qqf(numl)=c3nrr*erz + c3nsr*esz + c3nfr*efz
              qqg(numl)=c3nri*erz + c3nsi*esz + c3nfi*efz
             end if
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl1)=qqa(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl2)=qqb(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl5)=qqc(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl3)=qqd(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl6)=qqe(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl4)=qqf(numl)
       end do
       do numl = 1, numsnl(la,i_nod)
              vnlcr(numl+npvnl7)=qqg(numl)
       end do
   end if

   end if
   end do
   end do

end do
!==================== end ==============================================
#endif
!      ct = timecnt()
!      if( myid.eq.0 ) write(*,*) 'fnnlcl2-1:', ct - ct0
!      ct0 = ct


!do nspin = 1, nspnmx
!   if( lspin ) then
!       call stspud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
!       call ldocc( nspin, occ, wegud, nband )
!   end if

do ibb = 1, nbnod
   ib = iod(ibb)
   if( abs(occ(ib)).gt.1.d-15 ) then
!    --- FFT to obtain wavefunctions in r-space
   call wvg2r( nfile, myid, nodes,  &
& cgjr(1,ibb), nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, .true. )

   do i = 1, nhk2(ntype)
   do j = 1, nylmmx
      sumr(j,i) = svsumr(j,i,ibb,nspin,nk)
      sumxr(j,i) = 0.d0
      sumxx(j,i) = 0.d0
      sumxy(j,i) = 0.d0
   end do
   end do

   do i = 1, nionnew
      ia = iatmpt(i)
      it = iatoit(ia)
#ifndef VECTOR
      do l = 1, lmx(it)
         la = itola(l,it)
         do mm = 1, numsnl(la,i)
            m1 = meshnl(mm+npmsnl(la,i))
            eigfvx = vnlcr(mm+npvnlc(l,i))*fft3y(m1)
            sumxr(l,ia) = sumxr(l,ia) + eigfvx
         end do

         if( lcstress ) then
         if( igxxx.eq.1 ) then
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                xm = xm - cratm(1,i)
                ym = ym - cratm(2,i)
                eigfvx = vnlcr(mm+npvnlc(l,i))*fft3y(m1)
                sumxx(l,ia) = sumxx(l,ia) + eigfvx*xm
                sumxy(l,ia) = sumxy(l,ia) + eigfvx*ym
             end do
         else if( igxxx.eq.2 ) then
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3
                ym = ym - cratm(2,i)
                zm = zm - cratm(3,i)
                eigfvy = vnlcr(mm+npvnlc(l,i))*fft3y(m1)
                sumxx(l,ia) = sumxx(l,ia) + eigfvy*ym
                sumxy(l,ia) = sumxy(l,ia) + eigfvy*zm
             end do
         else
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                zm = zm - cratm(3,i)
                xm = xm - cratm(1,i)
                eigfvz = vnlcr(mm+npvnlc(l,i))*fft3y(m1)
                sumxx(l,ia) = sumxx(l,ia) + eigfvz*zm
                sumxy(l,ia) = sumxy(l,ia) + eigfvz*xm
             end do
         end if
         end if
      end do
#else
      do la = 1, lax(it)
         do mm = 1, numsnl(la,i)
            m1 = meshnl(mm+npmsnl(la,i))
            xxnl(mm) = fft3y(m1)
         end do
         if( lcstress ) then
         if( igxxx.eq.1 ) then
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                xmr(mm) = xm - cratm(1,i)
                ymr(mm) = ym - cratm(2,i)
             end do
         else if( igxxx.eq.2 ) then
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3
                ymr(mm) = ym - cratm(2,i)
                zmr(mm) = zm - cratm(3,i)
             end do
         else
             do mm = 1, numsnl(la,i)
                m1 = meshnl(mm+npmsnl(la,i))
                ix = mshgnx(m1)
                iy = mshgny(m1)
                iz = mshgnz(m1)
                d1 = ix - 1
                d2 = iy - 1
                d3 = iz - 1
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                zmr(mm) = zm - cratm(3,i)
                xmr(mm) = xm - cratm(1,i)
             end do
         end if
         end if

      do l = 1, lmx(it)
      if( itola(l,it) == la ) then
         do mm = 1, numsnl(la,i)
            eigfvx = vnlcr(mm+npvnlc(l,i))*xxnl(mm)
            sumxr(l,ia) = sumxr(l,ia) + eigfvx
         end do

         if( lcstress ) then
         if( igxxx.eq.1 ) then
             do mm = 1, numsnl(la,i)
                xm = xmr(mm)
                ym = ymr(mm)
                eigfvx = vnlcr(mm+npvnlc(l,i))*xxnl(mm)
                sumxx(l,ia) = sumxx(l,ia) + eigfvx*xm
                sumxy(l,ia) = sumxy(l,ia) + eigfvx*ym
             end do
         else if( igxxx.eq.2 ) then
             do mm = 1, numsnl(la,i)
                ym = ymr(mm)
                zm = zmr(mm)
                eigfvy = vnlcr(mm+npvnlc(l,i))*xxnl(mm)
                sumxx(l,ia) = sumxx(l,ia) + eigfvy*ym
                sumxy(l,ia) = sumxy(l,ia) + eigfvy*zm
             end do
         else
             do mm = 1, numsnl(la,i)
                zm = zmr(mm)
                xm = xmr(mm)
                eigfvz = vnlcr(mm+npvnlc(l,i))*xxnl(mm)
                sumxx(l,ia) = sumxx(l,ia) + eigfvz*zm
                sumxy(l,ia) = sumxy(l,ia) + eigfvz*xm
             end do
         end if
         end if
      end if
      end do
      end do
#endif
   end do

   do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
   do ia = nhk1(it), nhk2(it)
       do l = 1, lmx(it)
          sumxr(l,ia) = sumxr(l,ia) * rdelv
       end do
       if( lcstress ) then
           do l = 1, lmx(it)
              sumxx(l,ia) = sumxx(l,ia) * rdelv
              sumxy(l,ia) = sumxy(l,ia) * rdelv
           end do
       end if
   end do
   end if
   end do


   occib2   = occ(ib)*2.d0
   occib2bv = occib2 / rvol
   do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
   do ia = nhk1(it), nhk2(it)
      fcx = 0.d0
      stxx = 0.d0
      stxy = 0.d0
      do l = 1, lmx(it)
         fcx = fcx + sumr(l,ia)*sumxr(l,ia)
      end do
      if( lcstress ) then
          do l = 1, lmx(it)
             stxx = stxx + sumr(l,ia)*sumxx(l,ia)
             stxy = stxy + sumr(l,ia)*sumxy(l,ia)
          end do
      end if

      fnlc(igxxx,ia) = fnlc(igxxx,ia) + fcx*occib2bv
      if( lcstress ) then
          if( igxxx.eq.1 ) then
              bufst(1)   = bufst(1)   + stxx*occib2
              bufst(6)   = bufst(6)   + stxy*occib2
          else if( igxxx.eq.2 ) then
              bufst(2)   = bufst(2)   + stxx*occib2
              bufst(4)   = bufst(4)   + stxy*occib2
          else
              bufst(3)   = bufst(3)   + stxx*occib2
              bufst(5)   = bufst(5)   + stxy*occib2
          end if
      end if
   end do
   end if
   end do

   end if
end do

!end do
!      ct = timecnt()
!      if( myid.eq.0 ) write(*,*) 'fnnlcl2-2:', ct - ct0
!      ct0 = ct
!=======================================================================
end do

!do it = 1, ntype
!if( lkbppi(it) .and. lking(it) ) then
!do i = nhk1(it), nhk2(it)
!   fnlc(1,i) = fnlc(1,i)/rvol
!   fnlc(2,i) = fnlc(2,i)/rvol
!   fnlc(3,i) = fnlc(3,i)/rvol
!end do
!end if
!end do


!      ii = 0
!      do it = 1, ntype
!      if( lkbppi(it) .and. lking(it) ) then
!      do i = nhk1(it), nhk2(it)
!         ii = ii + 1
!         dbuf(3*ii-2) = fnlc(1,i)
!         dbuf(3*ii-1) = fnlc(2,i)
!         dbuf(3*ii-0) = fnlc(3,i)
!      end do
!      end if
!      end do
!
!      call gdsum(dbuf,3*ii,dbufr)
!
!      ii = 0
!      do it = 1, ntype
!      if( lkbppi(it) .and. lking(it) ) then
!      do i = nhk1(it), nhk2(it)
!         ii = ii + 1
!         fnlc(1,i) = dbuf(3*ii-2)
!         fnlc(2,i) = dbuf(3*ii-1)
!         fnlc(3,i) = dbuf(3*ii-0)
!      end do
!      end if
!      end do


    if( ltimecnt ) then
        ct = timecnt()
        do ii = 1, 2
        if( loutfile(ii) ) then
           write(nfile(ii),*) '                     nonlocal pp ',  &
&                             ':          :', ct - ct0 - t_comm
           write(nfile(ii),*) '            comm. in nonlocal pp ',  &
&                             ':          :', t_comm
        end if
        end do
        ct0 = ct
    end if

return
end




subroutine sphrb1( l2, c1, s1, c2, s2, ylmr, ylmi, mxl )
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
ylmr(3, 1, 1) = -a1*s2
ylmi(3, 1, 1) =  a1*c2
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
a1  = pi21*s1*c1
ac1 = pi21*( c12 - s12 )
ylmr(1, 1, 2) = a1*c2
ylmi(1, 1, 2) = a1*s2
ylmr(2, 1, 2) = ac1*c2
ylmi(2, 1, 2) = ac1*s2
ylmr(3, 1, 2) = -a1*s2
ylmi(3, 1, 2) =  a1*c2
ylmr(1,-1, 2) = -ylmr(1, 1, 2)
ylmi(1,-1, 2) =  ylmi(1, 1, 2)
ylmr(2,-1, 2) = -ylmr(2, 1, 2)
ylmi(2,-1, 2) =  ylmi(2, 1, 2)
ylmr(3,-1, 2) = -ylmr(3, 1, 2)
ylmi(3,-1, 2) =  ylmi(3, 1, 2)
!      a1 = sqrt( 15.d0/(32.d0*pi) )*s12
a1  = pi22*s12
ac1 = pi22*2.d0*s1*c1
ylmr(1, 2, 2) = a1*c22p
ylmi(1, 2, 2) = a1*s22p
ylmr(2, 2, 2) = ac1*c22p
ylmi(2, 2, 2) = ac1*s22p
ylmr(3, 2, 2) = -a1*2.d0*s22p
ylmi(3, 2, 2) =  a1*2.d0*c22p
ylmr(1,-2, 2) =  ylmr(1, 2, 2)
ylmi(1,-2, 2) = -ylmi(1, 2, 2)
ylmr(2,-2, 2) =  ylmr(2, 2, 2)
ylmi(2,-2, 2) = -ylmi(2, 2, 2)
ylmr(3,-2, 2) =  ylmr(3, 2, 2)
ylmi(3,-2, 2) = -ylmi(3, 2, 2)
if( l2.le.2 ) return


write(*,*) ' ***  l=3 is not supported yet .'
c13 = c12*c1
s13 = s12*s1
c23p = c22p*c2 - s22p*s2
s23p = s22p*c2 + c22p*s2
!      ylmr(1,0, 3) = sqrt( 7.d0/(16.d0*pi) )*c1*(5.d0*c12 - 3.d0)
ylmr(1,0, 3) = pi30*c1*(5.d0*c12 - 3.d0)
ylmi(1, 0, 3) = 0.d0
!      a1 = sqrt( 21.d0/(64.d0*pi) )*s1*(5.d0*c12 - 1.d0)
a1 = pi31*s1*(5.d0*c12 - 1.d0)
ylmr(1, 1, 3) = a1*c2
ylmi(1, 1, 3) = a1*s2
ylmr(1,-1, 3) = -ylmr(1,1,3)
ylmi(1,-1,3) =  ylmi(1, 1,3)
!      a1 = sqrt( 105.d0/(32.d0*pi) )*s12*c1
a1 = pi32*s12*c1
ylmr(1,2, 3) = a1*c22p
ylmi(1, 2, 3) = a1*s22p
ylmr(1,-2,3) =  ylmr(1,2,3)
ylmi(1, -2,3) = -ylmi(1, 2,3)
!      a1 = sqrt( 35.d0/(64.d0*pi) )*s13
a1 = pi33*s13
ylmr(1,3, 3) = a1*c23p
ylmi(1, 3, 3) = a1*s23p
ylmr(1,-3,3) = -ylmr(1,3,3)
ylmi(1, -3,3) =  ylmi(1, 3,3)
!cc      if( l2.le.3 ) return

return
end
