



module psvariables
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkb.f
!-----------------------------------------------------------------------
implicit none

real*8, parameter :: pi = 3.141592653589793d0

integer, parameter :: MSH = 421*5

integer :: NAE
integer, allocatable, dimension(:,:) :: MJTAB1, MJTAB2
real*8,  allocatable, dimension(:,:) :: factl
integer, parameter ::  jtab(7) = (/ 1, 3, 5, 7, 0, 0, 0 /)
integer, parameter :: rjtab(7) = (/ 1, 1, 3, 3, 5, 5, 7 /)

real*8  :: RMAX, DX
integer :: MESH
real*8,  allocatable, dimension(:) :: rmaxi, dxi
integer, allocatable, dimension(:) :: meshi

real*8, allocatable, dimension(:,:) :: rdial

real*8, allocatable, dimension(:,:) :: PLOCAL

real*8, allocatable, dimension(:)   :: PNL

real*8, allocatable, dimension(:,:) :: svrhops
real*8, allocatable, dimension(:,:) :: SGCAE

real*8, allocatable, dimension(:)   :: rad99

real*8, allocatable, dimension(:)     :: R
real*8, allocatable, dimension(:,:,:) :: PSUORG

save

end module




module dftCvariables
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkb.f
!-----------------------------------------------------------------------
implicit none

!-----for an empirical correction to non-local pp. (DFT+C)
logical :: lplusC = .false.
logical, allocatable, dimension(:) :: lplusC_at  ! on/off DFT+U for each atom
real*8,  allocatable, dimension(:,:) :: plusC_r    ! cutoff radius
real*8,  allocatable, dimension(:,:) :: plusC_e    ! energy shift

save

end module




module lopp_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkb.f & nlkbpp.f
!-----------------------------------------------------------------------
implicit none


real*8,  allocatable, dimension(:,:) :: svco, svcop
real*8,  allocatable, dimension(:)   :: rhmur
real*8,  allocatable, dimension(:)   :: svvext
save


end module




module nlpp_parameters
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkb.f, ppvd.f, nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
implicit none

integer :: nylmmx, nqlmmx
integer :: nrefxx, ncmxx, ncmtxx
integer :: lmmax, lpmax, lptmax, lfmax, lftmax
integer :: lnref

!-----for Q functions in ultrasoft pp
integer :: MXL2, MXL21
integer :: nqylmx
integer :: nqabm1, nqabm2
integer :: nqabmx
save


end module




module nlpp_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in ppkb.f, ppvd.f, nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
use nlpp_parameters
implicit none


integer :: nknod = 1   ! the number of k points
logical :: lgamma = .true.

logical :: ltddft, ltddft_fssh

logical :: lnoncollinear = .false.
integer :: ncscale = 1

integer :: nqemax, nqfmax
integer :: nqesum, nqfsum

integer, allocatable, dimension(:,:) :: lft2lf

integer, allocatable, dimension(:) :: lmx, lpmx, lptmx
integer, allocatable, dimension(:) :: lfmx, lftmx

integer, allocatable, dimension(:,:) :: iptv1, iptv2
integer, allocatable, dimension(:,:) :: ipv1,  ipv2
integer, allocatable, dimension(:,:) :: lpt2lp
real*8,  allocatable, dimension(:,:) :: clalp
real*8,  allocatable, dimension(:,:) :: qablm


integer, parameter :: ipking = 1
integer, parameter :: mxqlm = 550  ! dimension for table of nonlocal term
real*8,  allocatable, dimension(:,:,:) :: tbqlm, tbqlma
real*8,  allocatable, dimension(:)     :: dlqlm
real*8,  allocatable, dimension(:,:,:) :: tbqld, tbqlda
real*8,  allocatable, dimension(:)     :: dlqld

integer, allocatable, dimension(:) :: lax

real*8,  allocatable, dimension(:,:,:,:) :: qlmr, qlmi
real*8,  allocatable, dimension(:,:,:,:) :: qlmdgr, qlmdgi

integer, allocatable, dimension(:,:) :: itolm, itolmf, itola

integer, allocatable, dimension(:) :: lxxx
real*8,  allocatable, dimension(:) :: g1max

integer :: ikng1 = 0, ikng2 = 0
logical :: lyreal = .true.

logical, allocatable, dimension(:,:) :: litoll

real*8,  allocatable, dimension(:,:) :: bzk


!-----table for Q-functions
integer, parameter :: mxqab = 1100   ! dimension for table of Q-functions

real*8,  allocatable, dimension(:,:,:) :: tbqab, tbqaba
real*8,  allocatable, dimension(:)     :: dlqab
real*8,  allocatable, dimension(:,:,:) :: tbqad, tbqada
real*8,  allocatable, dimension(:)     :: dlqad
integer, allocatable, dimension(:)     :: nqabx

real*8,  allocatable, dimension(:,:)   :: qylmr, qylmi
real*8,  allocatable, dimension(:,:,:) :: qylmdr, qylmdi

integer, allocatable, dimension(:,:,:) :: ijtolm, ijtolf, ijtola
integer, allocatable, dimension(:,:)   :: ijtonm

real*8,  allocatable, dimension(:,:,:) :: cgijto

logical, allocatable, dimension(:,:,:) :: lijtll

integer :: lxxqab
real*8  :: g1maxq

real*8, allocatable, dimension(:,:) :: abgvec

save


end module




subroutine psvariables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype )
!-----------------------------------------------------------------------
!     allocate memory for pseudopotentials
!-----------------------------------------------------------------------
use psvariables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( rmaxi(ntype), dxi(ntype), meshi(ntype),  &
& rdial(MSH,ntype), PLOCAL(MSH,ntype), PNL(MSH),  &
& svrhops(MSH,ntype), SGCAE(MSH,ntype), rad99(ntype),  &
& stat=status )

the_mem =  &
&  8.d0 * ntype*2 + 4.d0 * ntype  &
& + 8.d0 * ( MSH*ntype*4 + MSH + ntype )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'psvariables_alloc', .true. )


return
end subroutine




subroutine psvariables_mxl_alloc( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     allocate memory for pseudopotentials
!-----------------------------------------------------------------------
use param
use param_atom
use psvariables
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: it, l, NNAE
real*8  :: dl
integer :: status
real*8  :: the_mem


NAE    = 2*(mxl + 1) - 1

!------allocate memory
allocate( R(MSH), PSUORG(MSH,NAE,ntype),  &
& MJTAB1(0:MXL,ntype), MJTAB2(0:MXL,ntype), factl(NAE,ntype),  &
& stat=status )

the_mem =  &
&   8d0 * ( size(R) + size(PSUORG) )  &
& + 4d0 * ( size(MJTAB1) + size(MJTAB2) )  &
& + 8d0 * size(factl)

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'psvariables_mxl_alloc', .true. )


do it = 1, ntype

    do L = 0, mxl
       MJTAB1(L,it) = 2 * L
       MJTAB2(L,it) = 2 * L + 1
    end do
    MJTAB1(0,it) = 1

    if( lsomax(it) >= 0 ) then
        !---relativistic data set
        do L = 0, mxl
           if( L == 0 ) then
               dl = L + 1
             else
               dl = L
           end if
           do NNAE = MJTAB1(L,it), MJTAB2(L,it)
              factl(NNAE,it) = dl/dble( 2*L + 1 )
              dl = dl + 1.d0
           end do
        end do
    else
        !---non-relativistic data set
        do L = 0, mxl
           MJTAB1(L,it) = MJTAB2(L,it)
           NNAE = MJTAB1(L,it)
           factl(NNAE,it) = 1.d0
        end do
    end if

end do


return
end subroutine




subroutine psvariables_dealloc( nfile, myid, nodes,  &
& alloc_mem, mxl, ntype )
!-----------------------------------------------------------------------
!     deallocate memory for pseudopotentials
!-----------------------------------------------------------------------
use psvariables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: mxl
integer :: ntype

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
deallocate( rmaxi, dxi, meshi, rdial, PLOCAL, PNL, svrhops,  &
& SGCAE, rad99, R, PSUORG, stat=status )

the_mem =  &
&  8.d0 * ntype*2 + 4.d0 * ntype  &
& + 8.d0 * ( MSH*ntype*4 + MSH + ntype )  &
& + 8.d0 * ( MSH + MSH*(mxl+1)*ntype )

!------error trap
call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'psvariables_dealloc', .true. )


return
end subroutine




subroutine lopp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nplw5, lstress, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for local pseudopotentials
!-----------------------------------------------------------------------
use lopp_variables
use planewave_decomp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer :: nplw5
logical :: lstress
real*8  :: pwscale

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: ngenh = 0
integer :: ngenhnod = 0
save ngenh, ngenhnod


if( nplw5+1  > ngenh     .or. &
&   nplw5nod > ngenhnod ) then

    !-----if already allocated, deallocate arrays
    if( allocated(svco) ) then

        the_mem = 8.d0 * ( size(svco) + size(rhmur) )

        !------deallocate memory
        deallocate( svco, rhmur, stat=status )

        if( allocated(svcop) ) then
            the_mem = the_mem + 8.d0 * size(svcop)
            deallocate( svcop, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'lopp_variables_alloc', .true. )

    end if

    ngenh    = (nplw5+1)* pwscale
    ngenhnod = nplw5nod * pwscale
    !------allocate memory
    allocate( svco(ngenhnod,ntype), rhmur(2*ngenh), stat=status )

    if( status == 0 .and. lstress )  &
&       allocate( svcop(ngenhnod,ntype), stat=status )

    the_mem = 8.d0 * ( size(svco) + size(rhmur) )
    if( lstress ) the_mem = the_mem + 8.d0 * size(svcop)

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'lopp_variables_alloc', .true. )

end if


return
end subroutine




subroutine svvext_alloc( nfile, myid, nodes,  &
& alloc_mem, llclpp_r, llclpp_g, lstress, mshnod, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for local pseudopotentials
!-----------------------------------------------------------------------
use lopp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
logical :: llclpp_r, llclpp_g, lstress
integer :: mshnod
real*8  :: pwscale

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: lplmshnod = 0
integer :: mshnodx
save lplmshnod


if( .not.lstress ) return

if( llclpp_r .and. llclpp_g ) then

    mshnodx = mshnod
    call gimax(mshnodx)

    !-----check array sizes
    if( mshnodx > lplmshnod ) then

        !-----if already allocated, deallocate arrays
        if( allocated(svvext) ) then
            deallocate( svvext, stat=status )
            the_mem = 8.d0 * ( lplmshnod )
            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'svvext_alloc', .true. )
        end if

        !------allocate memory
        lplmshnod = mshnodx * pwscale
        allocate( svvext(lplmshnod), stat=status )

        the_mem = 8.d0 * ( lplmshnod )

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'svvext_alloc', .true. )

    end if

end if


return
end subroutine




subroutine save_vext( nfile, myid, nodes, vext, mshnod )
!-----------------------------------------------------------------------
!     save local pp, if necessary
!-----------------------------------------------------------------------
use lopp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8,  dimension(mshnod) :: vext
integer :: i


if( allocated(svvext) ) then
    do i = 1, mshnod
       svvext(i) = vext(i)
    end do
end if


return
end subroutine




subroutine calc_vext( nfile, myid, nodes,  &
& elclv, elcl, rho, mshnod, rdelv )
!-----------------------------------------------------------------------
!     calculate local pp energy, if necessary
!-----------------------------------------------------------------------
use lopp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8  :: elclv, elcl
real*8,  dimension(mshnod) :: rho
real*8  :: rdelv
integer :: i
real*8  :: dbuf1r


if( allocated(svvext) ) then
    elclv = 0.d0
    do i = 1, mshnod
       elclv = elclv + svvext(i)*rho(i)
    end do
    call gdsum(elclv,1,dbuf1r)
    elclv = elclv * rdelv
  else
    elclv = elcl
end if


return
end subroutine




!subroutine set_param_for_kbpp( nfile, myid, nodes, mxl, mxref_ )
!!-----------------------------------------------------------------------
!!     allocate memory for shared work arrays in ppvd.f
!!-----------------------------------------------------------------------
!use nlpp_variables
!implicit none
!integer :: nfile(*), myid, nodes
!integer :: mxl
!integer :: mxref_
!
!
!!-----set constants
!nylmmx = mxl*(mxl+2) + 1
!nqlmmx = (mxl + 1)*(mxl + 2)/2
!nrefxx = mxref_
!ncmxx  = nrefxx*(nrefxx-1)/2 + nrefxx
!ncmtxx = nrefxx*nrefxx
!lmmax  = nrefxx*nylmmx
!lpmax  = ncmxx *nylmmx
!lptmax = ncmtxx*nylmmx
!lfmax  = lmmax*(lmmax + 1)/2
!lftmax = lmmax*lmmax
!lnref  = nrefxx*(mxl + 1)
!
!!-----Q functions are not needed for kbpp
!!      MXL2  = 2*mxl
!!      MXL21 = MXL2 + 1
!!      nqylmx = (MXL21 + 1)*(MXL21 + 2)/2
!!      nqabm1 = ncmtxx*MXL*(MXL + 1)*(MXL + 2)/6
!!      nqabm2 = ncmxx*(MXL + 1)*(MXL + 2)/2
!!      nqabmx = nqabm1 + nqabm2
!MXL2  = 0
!MXL21 = 0
!nqylmx = 0
!nqabm1 = 0
!nqabm2 = 0
!nqabmx = 0
!
!
!return
!end subroutine




subroutine nlpp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, nplw3, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for ultrasoft pseudopotentials
!-----------------------------------------------------------------------
use nlpp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer :: mxl
integer :: nplw3
real*8  :: pwscale

!-----declare local variables
integer :: i
integer :: status
real*8  :: the_mem
integer :: lpl3h = 0
save lpl3h


if( nplw3 > lpl3h .and. nqylmx > 0 ) then

    !-----if already allocated, deallocate arrays
    if( allocated(qylmr) ) then

        !------deallocate memory
        deallocate( qylmr, qylmi, qylmdr, qylmdi, stat=status )

        the_mem = 8.d0 * (lpl3h+1)*nqylmx*( 2 + 3*2 )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlpp_variables_alloc', .true. )

    end if

    lpl3h = nplw3 * pwscale
    !------allocate memory
    allocate(  &
& qylmr(0:lpl3h,nqylmx), qylmi(0:lpl3h,nqylmx),  &
& qylmdr(0:lpl3h,nqylmx,3), qylmdi(0:lpl3h,nqylmx,3),  &
& stat=status )

    the_mem = 8.d0 * (lpl3h+1)*nqylmx*( 2 + 3*2 )

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlpp_variables_alloc', .true. )

end if


!-----if already allocated, nothing is needed below
if( allocated(lft2lf) ) return


!------allocate memory
allocate( lft2lf(lftmax,ntype),  &
& lmx(ntype), lpmx(ntype), lptmx(ntype), lfmx(ntype), lftmx(ntype),  &
& iptv1(lftmax,ntype), iptv2(lftmax,ntype),  &
& ipv1(lfmax,ntype), ipv2(lfmax,ntype), lpt2lp(lptmax,ntype),  &
& clalp(lptmax,ntype), qablm(lptmax,ntype),  &
& lax(ntype),  &
& itolm(lmmax,ntype), itolmf(lmmax,ntype), itola(lmmax,ntype),  &
& litoll(lmmax,ntype),  &
& stat=status )

if( status == 0 .and. nqylmx > 0 ) allocate(  &
& tbqab(0:mxqab,nqabmx,ntype), tbqaba(0:mxqab,nqabmx,ntype),  &
& dlqab(ntype), tbqad(0:mxqab,nqabmx,ntype),  &
& tbqada(0:mxqab,nqabmx,ntype), dlqad(ntype), nqabx(ntype),  &
& ijtolm(MXL21,lfmax,ntype), ijtolf(MXL21,lfmax,ntype),  &
& ijtola(MXL21,lfmax,ntype), ijtonm(lfmax,ntype),  &
& cgijto(MXL21,lfmax,ntype), lijtll(MXL21,lfmax,ntype),  &
& stat=status )

the_mem =  &
&  4.d0 * ( lftmax*ntype + ntype*5 + lftmax*ntype*2  &
&         + lfmax*ntype*2 + lptmax*ntype )  &
& + 8.d0 * ( lptmax*ntype*2 )  &
& + 4.d0 * ( ntype )  &
& + 4.d0 * ( lmmax*ntype*3 )  &
& + 1.d0 * ( lmmax*ntype )

if( nqylmx > 0 ) then
    the_mem = the_mem  &
& + 8.d0 * ( (mxqab+1)*nqabmx*ntype*2 + ntype )*2  &
& + 4.d0 * ( ntype )  &
& + 4.d0 * ( MXL21*lfmax*ntype*3 + lfmax*ntype )  &
& + 8.d0 * ( MXL21*lfmax*ntype )  &
& + 1.d0 * ( MXL21*lfmax*ntype )
end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlpp_variables_alloc(2)', .true. )


!-----first-order SO effects by purturbation calculation
!call nlppSOp_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype )


return
end subroutine




subroutine nlppg_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, nplw, recnrm, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for ultrasoft pseudopotentials
!-----------------------------------------------------------------------
use nlpp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer :: mxl
integer :: nplw
real*8,  dimension(0:nplw) :: recnrm
real*8  :: pwscale

!-----declare local variables
integer :: i
integer :: status
real*8  :: the_mem
integer :: lplh = 0
real*8  :: pkscale
save lplh


!---for k-point sampling
call get_pkscale( pkscale )

if( nplw * pkscale > lplh ) then

    !-----if already allocated, deallocate arrays
    if( allocated(qlmr) ) then

        !------deallocate memory
        deallocate( qlmr, qlmi, qlmdgr, qlmdgi, abgvec,  &
& stat=status )

        the_mem =  &
& + 8.d0 * ( (lplh+1)*nqlmmx*nknod*( ipking*2 + 3*2 ) )  &
& + 8.d0 * ( (lplh+1)*nknod )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlppg_variables_alloc', .true. )

    end if

    lplh = nplw * pwscale * pkscale
    !------allocate memory
    allocate(  &
& qlmr(0:lplh,nqlmmx,nknod,ipking),  &
& qlmi(0:lplh,nqlmmx,nknod,ipking),  &
& qlmdgr(0:lplh,nqlmmx,3,nknod), qlmdgi(0:lplh,nqlmmx,3,nknod),  &
& abgvec(0:lplh,nknod),  &
& stat=status )

    the_mem =  &
& + 8.d0 * ( (lplh+1)*nqlmmx*nknod*( ipking*2 + 3*2 ) )  &
& + 8.d0 * ( (lplh+1)*nknod )

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlppg_variables_alloc', .true. )

end if


!-----if not allocated, allocate arrays
if( .not.allocated(lxxx) ) then

    !------allocate memory
    allocate(  &
& lxxx(ipking), g1max(ipking), bzk(3,nknod),  &
& tbqlm(0:mxqlm,lnref,ntype), tbqlma(0:mxqlm,lnref,ntype),  &
& dlqlm(ntype), tbqld(0:mxqlm,lnref,ntype),  &
& tbqlda(0:mxqlm,lnref,ntype), dlqld(ntype),  &
& stat=status )

    the_mem =  &
& + 4.d0 * ( ipking )  &
& + 8.d0 * ( ipking )  &
& + 8.d0 * ( 3*nknod )  &
& + 8.d0 * ( (mxqlm+1)*lnref*ntype*2 + ntype )*2

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlppg_variables_alloc(2)', .true. )

end if


do i = 1, ipking
   g1max(i) = 0.d0
end do
do i = 1, nknod
   bzk(1,i) = 0.d0
   bzk(2,i) = 0.d0
   bzk(3,i) = 0.d0
end do

do i = 0, nplw
   abgvec(i,1)  = sqrt(recnrm(i))
end do

!call get_bzk_and_norm( bzk, abgvec, lplh )


return
end subroutine




subroutine nlpp_settddft
!-----------------------------------------------------------------------
!     set parameters for TDDFT
!-----------------------------------------------------------------------
use nlpp_variables
implicit none

!---get ltddft
call get_ltddft( ltddft, ltddft_fssh )

return
end subroutine




subroutine setpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens, ecutsoft,  &
& lpaw, lpaw_sym, lkbpp, lvand, lvfull, llocl, lkbppi, lvandi,  &
& lvflag, llocli,  &
& lnlpp_r, lnlpp_g, lkbpp_r, lkbpp_g, lvand_r, lvand_g,  &
& llclpp_r, llclpp_g, llking,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& lking, rking, gkgmax, gkgexct, rctflc,  &
& lpcc_r, lpcc_g, lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, frac_rcomp_it,  &
& ms, aname, mxl, mx1, mx1loc, mxref, vlocli ,xitgrd, rr, nvlcl,  &
& nd1vks, rvol, nplw5, nplw, gx, gy, gz, recnrm,  &
& nplw3, gboxx, gboxy, gboxz, kfft1b, kfft2b, kfft3b, kfft0b,  &
& alloc_mem, lstress, ldouble_grid_recip, dgalpha, recnrmex,  &
& pwscale, lvshape, dkgx, dkgy, dkgz, dkrecnrm )
!-----------------------------------------------------------------------
!      local and non-local pseudopotential elements 
!                        by
!           relativistic pseudopotential 
!
!  based upon the work by 
!  N. Troullier and J. L. Martins :  Phys. Rev. B43(1991)1993.
!   and
!  G. B. Bachelet and M. Schluter :  Phys. Rev. B25(1982)2103.
!-----------------------------------------------------------------------
use nlpp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: jgga
logical :: lrela
integer :: ntype
integer :: mxl
integer :: mx1, mx1loc
integer :: mxref
real*8,  dimension(ntype) :: zatom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lclno
logical, dimension(0:mxl,ntype) :: lchk
real*8 :: ecutdens, ecutsoft
logical :: lpaw, lpaw_sym
logical ::                   lkbpp,  lvand,  lvfull, llocl
logical, dimension(ntype) :: lkbppi, lvandi, lvflag, llocli
logical :: lnlpp_r, lnlpp_g, lkbpp_r, lkbpp_g, lvand_r, lvand_g
logical :: llclpp_r, llclpp_g
logical, dimension(ntype) :: llking
real*8,  dimension(0:mx1loc,ntype) :: tablc, tablca
real*8,  dimension(ntype)          :: dltlc, rmxlc
real*8,  dimension(0:mx1loc,ntype) :: tbflc, tbflca
real*8,  dimension(ntype)          :: dltflc, rmxflc
logical, dimension(ntype) :: lking
real*8,  dimension(ntype) :: rking, gkgmax, gkgexct
real*8  :: rctflc
logical :: lpcc_r, lpcc_g
real*8,  dimension(ntype) :: rpcc
logical, dimension(ntype) :: lpcci
logical, dimension(ntype) :: lpking
real*8,  dimension(ntype) :: rpking, gpkgmax, gpkgexct
real*8,  dimension(ntype) :: frac_rcomp_it
character(1), dimension(7)   :: ms
character(2), dimension(103) :: aname
integer :: nvlcl
real*8,  dimension(0:nvlcl)  :: vlocli ,xitgrd, rr
integer, dimension(3) :: nd1vks
real*8  :: rvol
integer :: nplw
integer :: nplw5
real*8,  dimension(0:nplw5) :: gx, gy, gz
real*8,  dimension(0:nplw5) :: recnrm
integer :: nplw3
real*8,  dimension(0:nplw3) :: gboxx, gboxy, gboxz
integer :: kfft1b, kfft2b, kfft3b, kfft0b
real*8 :: alloc_mem
logical :: lstress
logical :: ldouble_grid_recip
real*8  :: dgalpha
real*8,  dimension(0:nplw5) :: recnrmex
real*8  :: pwscale
logical :: lvshape
real*8  :: dkgx(*), dkgy(*), dkgz(*), dkrecnrm(*)

!-----declare local variables
integer :: nplw_dummy(1)
integer :: ierror
logical :: licall = .true.
save licall


!-----set parameters for k-point sampling
!call nlpp_setkp

!-----set parameters for TDDFT
call nlpp_settddft

!-----set parameters for noncollinear magnetism
!call nlpp_setnoncollinear


!-----allocate memory for ultrasoft pseudopotential
!if( lvand ) then

    !-----allocate memory
!    call ppvd_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, mxl, mxref )

!  else
!
!    mxref = 1
!    call set_param_for_kbpp( nfile, myid, nodes, mxl, mxref )

!end if


!-----allocate memory for nonlocal pseudopotential
if( lnlpp_g .or. lnlpp_r .or. lvand ) then

    call nlpp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, nplw3, pwscale )

end if


!-----allocate memory for the calculation in real space
if( lnlpp_r ) then

    call nlppr_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, mx1, lkbpp_r, lvand_r )

end if


!-----allocate memory for the calculation in reciprocal space
if( lnlpp_g ) then

    if( .not.lnoncollinear ) then
        call nlppg_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, nplw, recnrm, pwscale )
    else
        call nlppg_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, nplw, dkrecnrm, pwscale )
    end if

end if


!-----load norm-conseving pseudopotentials
if( lkbpp .and. licall ) then
    call rdtmpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens, lkbppi,  &
& lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, ms, aname,  &
& mxl, mx1, vlocli, xitgrd, rr, mx1 )
end if


!if( lvand .and. licall ) then
!    if( lpaw ) then
!
!        !-----allocate memory for the PAW method
!        call paw_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, mxl, lpaw_sym )
!
!        !-----load data sets for the PAW method
!        call rdpaw( nfile, myid, nodes,  &
!& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens,  &
!& lvandi, lvflag, lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, frac_rcomp_it,  &
!& ms, aname, mxl, mx1, vlocli, xitgrd, rr, mx1 )
!
!      else
!
!        !-----load ultrasoft pseudopotentials
!        call rdvand( nfile, myid, nodes,  &
!& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens,  &
!& lvandi, lvflag, lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct,  &
!& ms, aname, mxl, mx1, vlocli, xitgrd, rr, mx1 )
!
!    end if
!end if


!-----set recnrmex
!----- = 4*pi/recnrm in the usual method
!----- or
!----- = 4*pi/recnrm * ( 1 - exp(-recnrm/(4*a*a)) )
!----- in the double-grid method
call set_recnrmex( nfile, myid, nodes,  &
& recnrmex, recnrm, nplw5, ldouble_grid_recip, dgalpha )



!-----set local pseudopotentials & partial core correction
if( lkbpp .or. lvand .or. llocl ) then

    !-----set cutoff length for local pp,
    !-----and tables for the calculations in real space
    !---if( llclpp_r ) then
    if( licall ) then
        call lotmpp( nfile, myid, nodes,  &
& ntype, zv, rctflc, llking, tablc, tablca, dltlc, rmxlc,  &
& tbflc, tbflca, dltflc, rmxflc, mx1loc )
    end if
    !---end if

    !-----in reciprocal space
    if( llclpp_g ) then

        !-----allocate memory
        call lopp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nplw5, lstress, pwscale )

        call lotmpp_g( nfile, myid, nodes,  &
& ntype, zv, llking, dltlc, rmxlc,  &
& recnrm, recnrmex, nplw5, vlocli, xitgrd, rr, nvlcl, lstress,  &
& lvshape )
    end if


    !-----partial core correction in reciprocal space
    if( lpcc_g ) then

        !-----allocate memory
        call pcc_variables_alloc2( nfile, myid, nodes,  &
& alloc_mem, ntype, nplw5, lstress, pwscale )

        call subpcc_g( nfile, myid, nodes,  &
& ntype, lpcci, lpking,  &
& recnrm, recnrmex, nplw5, vlocli, xitgrd, rr, nvlcl, lstress,  &
& lvshape )
    end if
end if


!-----set norm-conseving pseudopotentials in KB form
if( lkbpp ) then

    !-----in real space
    if( lkbpp_r .and. licall ) then
        call nlkbpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutsoft, lkbppi,  &
& lking, rking, gkgmax, gkgexct, ms, aname,  &
& mxl, mx1, vlocli, xitgrd, rr, mx1 )
    end if

    !-----in reciprocal space
    if( lkbpp_g ) then
!        if( lgamma ) then
            nplw_dummy(1) = nplw ! just to avoid an attribute error
!            if( .not.lnoncollinear ) then
                call nlkbpp_g( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, lkbppi, lking,  &
& ms, aname, mxl, nplw_dummy, gx, gy, gz, nplw,  &
& vlocli, xitgrd, rr, mx1, ierror )
!            else
!                call nlkbpp_g( nfile, myid, nodes,  &
!& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, lkbppi, lking,  &
!& ms, aname, mxl, nplw_dummy, dkgx, dkgy, dkgz, nplw,  &
!& vlocli, xitgrd, rr, mx1, ierror )
!            end if
!        else
!            call nlkbpp_g_k( nfile, myid, nodes,  &
!& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, lkbppi, lking,  &
!& ms, aname, mxl, vlocli, xitgrd, rr, mx1,  &
!& ierror )
!        end if
    end if

    !-----first-order SO effects by purturbation calculation
!    call nlkbppSOp_set( nfile, myid, nodes,  &
!& ntype, lmax, lclno, lchk, lkbppi, lking, mxl )

end if


!-----set ultrasoft pseudopotentials
!if( lvand ) then
!
!    if( .not.lnoncollinear ) then
!        call nlvand( nfile, myid, nodes,  &
!& lpaw, ntype, lmax, lchk, lvandi, lvflag, lvand_r, lvand_g,  &
!& lking, rking, gkgmax, gkgexct, ecutsoft,  &
!& mxl, nd1vks, rvol, nplw, gx, gy, gz,  &
!& nplw3, gboxx, gboxy, gboxz, kfft1b, kfft2b, kfft3b, kfft0b,  &
!& mx1, vlocli, xitgrd, rr, mx1 )
!    else
!        call nlvand( nfile, myid, nodes,  &
!& lpaw, ntype, lmax, lchk, lvandi, lvflag, lvand_r, lvand_g,  &
!& lking, rking, gkgmax, gkgexct, ecutsoft,  &
!& mxl, nd1vks, rvol, nplw, dkgx, dkgy, dkgz,  &
!& nplw3, gboxx, gboxy, gboxz, kfft1b, kfft2b, kfft3b, kfft0b,  &
!& mx1, vlocli, xitgrd, rr, mx1 )
!    end if
!
!    !-----first-order SO effects by purturbation calculation
!    call nlvandSOp_set( nfile, myid, nodes,  &
!& ntype, lmax, lchk, lvandi, mxl )
!
!    !-----deallocate memory
!          call ppvd_variables_dealloc( nfile, myid, nodes,
!     & alloc_mem, ntype, mxl )
!end if


licall = .false.

return
end




subroutine pptmpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, lkbppi, zatom, zv, lmax, lsomax, lclno,  &
& lchk, mxlx, aname )
!-----------------------------------------------------------------------
!     Troullier and Martins pseudopotential
!-----------------------------------------------------------------------
use psvariables
implicit none
integer :: nfile(*), myid, nodes
integer :: jgga, ntype
logical :: lrela
integer :: mxlx
logical, dimension(ntype) :: lkbppi
real*8,  dimension(ntype) :: zatom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lsomax
integer, dimension(ntype) :: lclno
logical, dimension(0:mxlx,ntype) :: lchk
character(2), dimension(103) :: aname

!-----declare local variables
integer :: j, l, ldd, jmax
character(50) :: fname
integer :: jpl, ifpp
integer :: ierror, istat


call allocate_unit_number( jpl )


ifpp = 1
if( lrela ) ifpp = ifpp + 3

typedo: do j = 1, ntype
typeif: if( lkbppi(j) ) then


   fname = 'A.rvs'
   call openpp( nfile, myid, nodes,  &
& ifpp, jpl, zatom(j), fname, jgga, aname, ierror )

   blockif_a1: if( ierror == 0 ) then

!----  zv     : core charge
!----  lmax   :
!----  lsomax :
!----  lclno  : l for local potential
      read(JPL,*,iostat=istat) zv(j), lmax(j), lsomax(j), lclno(j)
      !----- error trap
      if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                      'error 0001 in pptmpp' )

      !-----error trap
      if( lmax(j) > mxlx ) call fstop( nfile, myid, nodes,  &
&                       'lmax(j) > mxlx in subroutine pptmpp' )

      !-----error trap : check LDA
      if( jgga == 1 .and. abs(lsomax(j)) /= 1 )  &
&                          call fstop( nfile, myid, nodes,  &
&                       'error-LDA in subroutine pptmpp' )

      !-----error trap : check PBE
      if( jgga == 2 .and. abs(lsomax(j)) /= 4 )  &
&                          call fstop( nfile, myid, nodes,  &
&                       'error-PBE in subroutine pptmpp' )

      !-----error trap : check RPBE
      if( jgga == 3 .and. abs(lsomax(j)) /= 5 )  &
&                          call fstop( nfile, myid, nodes,  &
&                       'error-RPBE in subroutine pptmpp' )

      read(JPL,*,iostat=istat) meshi(j), rmaxi(j), dxi(j)
      if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                      'error 0002 in pptmpp' )

      do l = 0, lmax(j)
         read(JPL,*,iostat=istat) ldd, jmax
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error 0003 in pptmpp' )

         lchk(l,j) = jmax .ge. 1
      end do

   else blockif_a1

     !----- error trap
     call fstop( nfile, myid, nodes, 'open file error in pptmpp' )

   end if blockif_a1
   CLOSE(JPL)

end if typeif
end do typedo

call deallocate_unit_number( jpl )


return
end




subroutine rdtmpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens, lkbppi,  &
& lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, ms, aname,  &
& mxl, mx1, vlocli ,xitgrd, rr, nvlcl )
!-----------------------------------------------------------------------
!     read pseudopotentials
!-----------------------------------------------------------------------
use outfile
use psvariables
use dftCvariables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
logical :: lrela
real*8,  dimension(ntype) :: zatom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: lmax
integer, dimension(ntype) :: lclno
logical, dimension(0:mxl,ntype) :: lchk
logical, dimension(ntype) :: lkbppi
character(1), dimension(7)   :: ms
character(2), dimension(103) :: aname
integer :: nvlcl
real*8,  dimension(0:nvlcl)  :: vlocli ,xitgrd, rr

dimension rpcc(*)
logical   lpcci(*)
logical   lpking(*)
dimension rpking(*), gpkgmax(*), gpkgexct(*)

!-----declare local variables
real*8,  allocatable, dimension(:) :: rhoc, rhops, VR, work,  &
&                                     work2, work3, work4
integer, parameter :: noofdt = 200
real*8,  allocatable, dimension(:) :: betaq
real*8,  allocatable, dimension(:) :: betar
real*8,  allocatable, dimension(:,:) :: amatrx
real*8,  allocatable, dimension(:,:) :: vmatrx
integer, allocatable, dimension(:)   :: ivmatrx
integer :: status
real*8  :: alocmem
character(50) :: FNAME
character(1)  :: DUMMY
character(1), dimension(0:9) :: NUM =  &
&             (/ '0','1','2','3','4','5','6','7','8','9' /)
integer :: jpl


!------allocate memory for local variables
allocate( rhoc(MSH), rhops(MSH), VR(MSH), work(MSH), work2(MSH),  &
& work3(MSH), work4(MSH),  &
& betaq(0:noofdt), betar(0:noofdt),  &
& amatrx(0:noofdt,0:noofdt), vmatrx(noofdt,noofdt+1),  &
& ivmatrx(noofdt), stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in rdtmpp' )

alocmem = 8.d0 * msh*7  &
&       + 8.d0 * ( noofdt + 1 ) * ( 2 + noofdt + 1 + noofdt )  &
&       + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, allocated (rdtmpp)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, allocated (rdtmpp)'


call allocate_unit_number( jpl )

PAI = ACOS(-1.D0)
typedo: do it = 1, ntype
typeif: if( lkbppi(it) ) then
!=======================================================================

!--- prepare radial meshes ---------------------------------------------
MESH =  meshi(it)
RMAX =  rmaxi(it)
DX   =  dxi(it)

D    =  EXP( DX )
CR   =  RMAX / D ** MESH
do I = 1, MESH
   CR     =  CR * D
   R(I)   =  CR
   rdial(I,it)   =  CR
end do
!-----------------------------------------------------------------------

ifpp = 1
if( lrela ) ifpp = ifpp + 3

do l=0, lmax(it)
if( lchk(l,it) ) then
do j = MJTAB1(l,it), MJTAB2(l,it)
    fname = (( ms(l+1)//num(rjtab(j)) )//'.rvs' )
    call openpp( nfile, myid, nodes,  &
& ifpp, jpl, zatom(it), fname, jgga, aname, ierror )

   do
      READ(JPL,'(a1)') DUMMY
      IF( DUMMY /= '#' ) exit
   end do
   backspace JPL
   i = 1
!   READ(JPL,2000) RI, ps1, ps2
!   PSUORG(I,j,it) = ( DBLE(l)*PS1 + DBLE(l+1)*PS2 )/DBLE(2*l+1)
   do
      READ(JPL,*,iostat=istat) RI, ps1
      if( istat /= 0 ) exit
      PSUORG(i,j,it) = ps1
!-->     PSUORG(I,j,it) = PSUORG(I,j,it)*.5d0
!c---                    [Ryd.] -> [Hartree]
      i = i + 1
      IF( i > MESH ) exit
   end do
   PSUORG(i:MESH,j,it) = 0.0d0
   CLOSE(JPL)
end do
end if
end do
!-----------------------------------------------------------------------
!---  read core- and pseudo-charge densities of reference state
!---                      calculated by relativistic atomic calculation.
   fname = 'val'
   call openpp( nfile, myid, nodes,  &
& ifpp, jpl, zatom(it), fname, jgga, aname, ierror )

   do
      READ(JPL,'(a1)') DUMMY
      IF( DUMMY /= '#' ) exit
   end do
   backspace JPL
   i = 1
   do
      READ(JPL,*,iostat=istat) ri, rhov, rhoc(i), rhops(i)
      if( istat /= 0 ) exit
      i = i + 1
      IF( i > MESH ) exit
   end do
   do j = i, MESH
      rhoc(j) = 0.d0
      rhops(j) = 0.d0
   end do
   CLOSE(JPL)
   do j = 1, MESH
      svrhops(j,it) = rhops(j)
      SGCAE(j,it)   = rhoc(j)
   end do
!-----------------------------------------------------------------------

!        ------------ partial core charge ---------------------------
   if( lpcci(it) ) then
       rsgc = rpcc(it)
       CALL PCORE3(rhoc, aaa, bbb, R, rsgc, RMAX, DX, MESH, INDER)
!---      error trap
       call gimax( INDER )
       if( INDER.ne.0 ) go to 998
       call subpcc( nfile, myid, nodes,  &
&                   mx1, ntype, lpcci, rpcc, rhoc, aaa, bbb,  &
&                   ecutdens, lpking, rpking, gpkgmax, gpkgexct,  &
&                   r, pai, it, mesh,  &
&                   vlocli, xitgrd, rr, nvlcl,  &
& betaq, betar, amatrx, vmatrx, ivmatrx, noofdt )
     else
       do i = 1, MESH
          rhoc(i) = 0.d0
       end do
   end if
!        ------------------------------------------------------------

!        ------ potential function from charge density --------------
   CALL POTVAL( rhops, rhoc, VR, QE, work, work2, work3, work4,  &
&               jgga, R, DX, MESH )
!        ------------------------------------------------------------
!   do l=0, lmax(it)
!   if( lchk(l,it) ) then
!       DO I = 1, MESH
!          PSUORG(I,jtab(l+1),it) = PSUORG(I,jtab(l+1),it) - 2.d0*VR(I)/R(I)
!       end do
!   end if
!   end do
!   if(loutfile(1)) write(nfile(1),*) ' check : No. of valence electrons :', QE
!   if(loutfile(2)) write(nfile(2),*) ' check : No. of valence electrons :', QE

!--- calcu. radius in which 99% of electrons are included.
   call calrad( nfile, myid, nodes,  &
& rhops, QE, work, radius, R, DX, MESH )
   rad99(it) = radius
!   if(loutfile(1)) write(nfile(1),*) ' radius for 99% valence electrons :', radius
!   if(loutfile(2)) write(nfile(2),*) ' radius for 99% valence electrons :', radius


!--- set local potential
   if( lrela ) then
       PLOCAL(1:MESH,it) = 0d0
       do j = MJTAB1(lclno(it),it), MJTAB2(lclno(it),it)
          PLOCAL(1:MESH,it) = PLOCAL(1:MESH,it) + PSUORG(1:MESH,j,it)*factl(j,it)
       end do
       !---if lclno(it) /= 0, the nonlocal potentials with lclno(it) are not zero
       if( lclno(it) /= 0 ) lclno(it) = lmax(it) + 1
   else
       do I = 1, MESH
          PLOCAL(I,it) = PSUORG(I,jtab(lclno(it)+1),it)
       end do
   end if

!--- set nonlocal potential
   do l=0, lmax(it)
   if( lchk(l,it) ) then
       do j = MJTAB1(l,it), MJTAB2(l,it)
          PSUORG(1:MESH,j,it) = PSUORG(1:MESH,j,it) - PLOCAL(1:MESH,it)
       end do
   end if
   end do

   DO I = 1, MESH
      PLOCAL(I,it) = PLOCAL(I,it) - 2.d0*VR(I)/R(I)
   end do

   !---DFT+C: an empirical correction to nonlocal pp.
!   if( lplusC ) then
!       do l=0, lmax(it)
!       if( lchk(l,it) ) then
!           do j = MJTAB1(l,it), MJTAB2(l,it)
!              CALL dftC_for_NCPP( nfile, myid, nodes,  &
!& MESH, PSUORG(1,j,it), msh,  &
!& lplusC_at(it), plusC_r(l,it), plusC_e(l,it), R )
!           end do
!       end if
!       end do
!   end if

!=======================================================================
end if typeif
end do typedo
! 2000 FORMAT(4D18.10)


!------deallocate arrays
deallocate( rhoc, rhops, VR, work, work2, work3, work4,  &
& betaq, betar, amatrx, vmatrx, ivmatrx, stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in rdtmpp' )

dealocmem = 8.d0 * msh*7  &
&       + 8.d0 * ( noofdt + 1 ) * ( 2 + noofdt + 1 + noofdt )  &
&       + 4.d0 * noofdt
if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d3),' KB, deallocated (rdtmpp)'
if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d3),' KB, deallocated (rdtmpp)'

call deallocate_unit_number( jpl )


return

998 CONTINUE
if(loutfile(1)) write(nfile(1),*) 'error : INDER in pcore : ', INDER
if(loutfile(2)) write(nfile(2),*) 'error : INDER in pcore : ', INDER
!-----Finalize the parallel environment
call end_parallel(ierr)
stop

999 CONTINUE
if(loutfile(1)) write(nfile(1),*) 'error : file (',FNAME,') is not found.'
if(loutfile(2)) write(nfile(2),*) 'error : file (',FNAME,') is not found.'
!-----Finalize the parallel environment
call end_parallel(ierr)
stop
end




SUBROUTINE POTVAL( SG, SGC, VR, QE, WORK, rho, drhox, d2rhox,  &
&                  jgga, R, DX, MESH )
!-----------------------------------------------------------------------
!                                                          1993/11/15
!      calculating the potential by valence charge density
!      [ a.u.] units are used ( = Ryd * 0.5 )
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION  SG(MESH), SGC(MESH)
DIMENSION  VR(MESH), WORK(MESH), R(MESH), rho(MESH), drhox(MESH),  &
&          d2rhox(MESH)
DIMENSION  DA(3), DB(3)

pi =  3.14159265358979323846d0

S2  =  LOG( SG(2) / SG(1) ) / DX
DO 100 I = 1, 2
   VR(I)   =  SG(I) * R(I) / ( S2 + 1.D0 )
   WORK(I) =  SG(I) / S2
   DB(I) =  DX * SG(I) / 3.D0
   DA(I) =  DB(I) * R(I)
100 CONTINUE
DO 110 I = 3, MESH
   DB(3) = DX * SG(I) / 3.D0
   DA(3) = DB(3) * R(I)
   VR(I)   =  VR(I-2)   + DA(3) + 4.D0 * DA(2) + DA(1)
   WORK(I) =  WORK(I-2) + DB(3) + 4.D0 * DB(2) + DB(1)
   DO 110 L = 1, 2
      DA(L) =  DA(L+1)
      DB(L) =  DB(L+1)
110 CONTINUE
QE  =  VR(MESH)
BM  =  WORK(MESH)
if( jgga.ge.2 ) then
!--- GGA ---------------------------------------------------------------
    do I = 1, MESH
       rho(i) = ( SG(i) + SGC(I) )/( 4.D0 * pi * R(I) * R(I) )
    end do
    call diff12( rho, MESH, DX, drhox, d2rhox )
    do I = 1, MESH
    IF( ( SG(I) + SGC(I) ) .LT. 1.D-30 .or. I.gt.mesh-3 ) THEN
       vxr = 0.d0
      ELSE
        drho   = drhox(i) / R(I)
        d2rho  = ( d2rhox(i) - drhox(i) ) / ( R(I) * R(I) )
!             if( jgga.eq.1 ) then
!                CALL vxcbp(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)
!               else if( jgga.eq.2 ) then
!                CALL pw91(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)
!               else
!       if( jgga.eq.2 ) then
!           CALL pbe(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)
           CALL normal_pbe(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)
!         else
!           CALL rpbe(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)
!       end if
       vxr = R(I) * ( VX + VC )
    end if
      VHARI =  VR(I) + R(I) * ( BM - WORK(I) )
      VR(I) = VHARI + vxr
    end do
!-----------------------------------------------------------------------
else
DO 120 I = 1, MESH
   SGI = SG(I) + SGC(I)
   IF( SGI .LT. 1.D-30 ) THEN
       vxr = 0.d0
    ELSE
       RS  =  ( 3.D0 * R(I) * R(I) / SGI )**( 1.D0 / 3.D0 )
       call vxc( rs, ex, vx, ec, vc )
       vxr = ( vx + vc ) * R(I)
   end if
   VHARI = VR(I) + R(I) * ( BM - WORK(I) )
   VR(I) = VHARI + vxr
120 CONTINUE
end if

!check
!      WRITE(*,*) 'information of valence charge density'
!      WRITE(*,*) ' *** QE',QE
!      WRITE(*,*) ' '


RETURN
END




SUBROUTINE diff12( f, mesh, DX, f1, f2 )
!-----------------------------------------------------------------------
!                                                          1997/03/13
!   calculation of first (f1) and second (f2) derivatives of f(i)
!      f(i) is f(x) at x = x0(arb.) + i*DX
!
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION      f(*), f1(*), f2(*)
dimension c1(5,3),c2(5,3)
save c1, c2
data c1/- 25.0d0,  48.0d0,- 36.0d0,  16.0d0,-  3.0d0,  &
&       -  3.0d0,- 10.0d0,  18.0d0,-  6.0d0,   1.0d0,  &
&          1.0d0,-  8.0d0,   0.0d0,   8.0d0,-  1.0d0/
data c2/  35.0d0,-104.0d0, 114.0d0,- 56.0d0,  11.0d0,  &
&         11.0d0,- 20.0d0,   6.0d0,   4.0d0,-  1.0d0,  &
&       -  1.0d0,  16.0d0,- 30.0d0,  16.0d0,-  1.0d0/


do i = 1, mesh
   if(      i.eq.1 ) then
         f1(i) = ( c1(1,1)*f(1) + c1(2,1)*f(2) + c1(3,1)*f(3)  &
&                + c1(4,1)*f(4) + c1(5,1)*f(5) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( c2(1,1)*f(1) + c2(2,1)*f(2) + c2(3,1)*f(3)  &
&                + c2(4,1)*f(4) + c2(5,1)*f(5) )  &
&                                          / ( 12.d0 * DX * DX )
   else if( i.eq.2 ) then
         f1(i) = ( c1(1,2)*f(1) + c1(2,2)*f(2) + c1(3,2)*f(3)  &
&                + c1(4,2)*f(4) + c1(5,2)*f(5) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( c2(1,2)*f(1) + c2(2,2)*f(2) + c2(3,2)*f(3)  &
&                + c2(4,2)*f(4) + c2(5,2)*f(5) )  &
&                                          / ( 12.d0 * DX * DX )
   else if( i.eq.3 ) then
         f1(i) = ( f(i-2) - 8.D0*f(i-1)  &
&                - f(i+2) + 8.D0*f(i+1) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( -f(i-2) + 16.D0*f(i-1) - 30.d0*f(i)  &
&                  -f(i+2) + 16.D0*f(i+1) )  &
&                                          / ( 12.d0 * DX * DX )
   else if( i.le.mesh-3 ) then
         f1(i) = ( -f(I-3) + 9.d0*f(I-2) - 45.D0*f(I-1)  &
&                +  f(I+3) - 9.d0*f(I+2) + 45.D0*f(I+1) )  &
&                                     / ( 60.D0 * DX )
         f2(i) = ( 2.d0*f(i-3) - 27.d0*f(I-2) + 270.D0*f(I-1)  &
&                   - 490.d0*f(I)  &
&             + 2.d0*f(i+3) - 27.d0*f(I+2) + 270.D0*f(I+1))  &
&                                     / ( 180.d0 * DX * DX )
   else if( i.eq.mesh-2 ) then
         f1(i) = ( f(i-2) - 8.D0*f(i-1)  &
&                - f(i+2) + 8.D0*f(i+1) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( -f(i-2) + 16.D0*f(i-1) - 30.d0*f(i)  &
&                  -f(i+2) + 16.D0*f(i+1) )  &
&                                          / ( 12.d0 * DX * DX )
   else if( i.eq.mesh-1 ) then
         f1(i) = ( -c1(5,2)*f(mesh-4) - c1(4,2)*f(mesh-3)  &
&                 - c1(3,2)*f(mesh-2)  &
&                 - c1(2,2)*f(mesh-1) - c1(1,2)*f(mesh) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( c2(5,2)*f(mesh-4) + c2(4,2)*f(mesh-3)  &
&                + c2(3,2)*f(mesh-2)  &
&                + c2(2,2)*f(mesh-1) + c2(1,2)*f(mesh) )  &
&                                          / ( 12.d0 * DX * DX )
   else
         f1(i) = ( -c1(5,1)*f(mesh-4) - c1(4,1)*f(mesh-3)  &
&                 - c1(3,1)*f(mesh-2)  &
&                 - c1(2,1)*f(mesh-1) - c1(1,1)*f(mesh) )  &
&                                          / ( 12.D0 * DX )
         f2(i) = ( c2(5,1)*f(mesh-4) + c2(4,1)*f(mesh-3)  &
&                + c2(3,1)*f(mesh-2)  &
&                + c2(2,1)*f(mesh-1) + c2(1,1)*f(mesh) )  &
&                                          / ( 12.d0 * DX * DX )
   end if
end do


RETURN
END




SUBROUTINE PCORE3(  &
& SGC, sgc0, sgcd0, R, R0, RMAX, DX, MESH, INDER )
!-----------------------------------------------------------------------
!     partial core correction
!
!        rho_pcc(r) = a * sin(q*r)/(q*r) ) + b*cos(q_0*r)
!
!          a, q, and b are chosen so that the first two derivatives
!          are continuous.
!
! (input)
!      SGC  : AE core density
!      R0   : cutoff length
!
! (output)
!      SGC   : partial core density, 4*pi*r*r * rho_pcc(r)
!      sgc0  : rho_pcc(r=0)
!      sgcd0 : (1/r)*( d rho_pcc/d r ) (r=0)
!
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
real*8, parameter :: pai =  3.14159265358979323846d0
DIMENSION      R(*), SGC(*)


INDER = 0
IF( R0.LT.-1.D-10 ) THEN
    DO 200 I = 1, MESH
       SGC(I) = 0.D0
200     CONTINUE
    RETURN
end if
IF( R0.LT.R(3)+1.D-10 ) RETURN

!--- constnats for numeircal differentiations --------------------------
D    =  EXP( DX )
CR   =  RMAX / D ** MESH
KI  = NINT( LOG(R0/CR)/DX )
IF( KI.LT.3 ) RETURN

!      PAI = ACOS(-1.D0)
!   --- for the first derivative ---
CST1  = 1.0D0/60.0D0/DX
CST11 =  CST1
CST12 = -CST1*9.0D0
CST13 =  CST1*45.0D0
CST15 = -CST1*45.0D0
CST16 =  CST1*9.0D0
CST17 = -CST1

!   --- for the second derivative ---
CST2  = 1.0D0/180.0D0/DX/DX
CST21 =  CST2*2.d0
CST22 = -CST2*27.d0
CST23 =  CST2*270.d0
CST24 = -CST2*490.d0
CST25 =  CST2*270.d0
CST26 = -CST2*27.d0
CST27 =  CST2*2.d0

!-----   numeircal differentiations of SGC   ---------------------------
PNL11 = CST11*SGC(KI+3)/R(KI+3)/R(KI+3)  &
&     + CST12*SGC(KI+2)/R(KI+2)/R(KI+2)  &
&     + CST13*SGC(KI+1)/R(KI+1)/R(KI+1)  &
&     + CST15*SGC(KI-1)/R(KI-1)/R(KI-1)  &
&     + CST16*SGC(KI-2)/R(KI-2)/R(KI-2)  &
&     + CST17*SGC(KI-3)/R(KI-3)/R(KI-3)
PNL1  = PNL11 / R(KI)

PNL21 = CST21*SGC(KI+3)/R(KI+3)/R(KI+3)  &
&     + CST22*SGC(KI+2)/R(KI+2)/R(KI+2)  &
&     + CST23*SGC(KI+1)/R(KI+1)/R(KI+1)  &
&     + CST24*SGC(KI  )/R(KI  )/R(KI  )  &
&     + CST25*SGC(KI-1)/R(KI-1)/R(KI-1)  &
&     + CST26*SGC(KI-2)/R(KI-2)/R(KI-2)  &
&     + CST27*SGC(KI-3)/R(KI-3)/R(KI-3)
PNL2  = ( PNL21 - PNL11 )/ ( R(KI)*R(KI) )

PNL0  = SGC(KI)/( 4.D0*PAI*R(KI)*R(KI) )
PNL1  = PNL1 / ( 4.D0*PAI )
PNL2  = PNL2 / ( 4.D0*PAI )


DVLG = R(KI) * PNL1/PNL0 + 1.D0

IF( DVLG.LE.1.D0 .AND. DVLG.GE.0.D0 ) THEN
    XMN = 0.D0
    XMX = PAI*.5D0
    X1  = 0.D0
    W1  = 1.D0
ELSE IF( DVLG.LT.0.D0 ) THEN
    XMN = PAI*.5D0
    XMX = PAI
    X1  = XMN
    W1  = 0.D0
ELSE
    XMN = PAI
    XMX = PAI*1.5D0
    X1  = XMX
    W1  = 0.D0
end if

ICHK = 0
X2 = ( XMN + XMX ) * .5D0
do
   W2 = X2*COS(X2)/SIN(X2)
!      WRITE(*,*) ICHK,X2,W2,DVLG
   ICHK = ICHK + 1
   IF( ICHK.GT.1000 ) THEN
       INDER = 16
       RETURN
   end if
   IF( ABS( W2 - DVLG ).LT.1.0D-12 ) exit
   CN = X1 + ( DVLG - W1 )*( X2 - X1 )/( W2 - W1 )
   IF( CN.LT.XMN ) CN = ( XMN + X2 )*.5D0
   IF( CN.GT.XMX ) CN = ( XMX + X2 )*.5D0
   X1 = X2
   W1 = W2
   X2 = CN
end do

BBB = X2/R(KI)
AAA = PNL0 * X2/SIN(X2)


q0 = 0.5d0*PAI/R(KI)
cb = ( PNL2 + ( 2.d0*PNL1+x2*x2*PNL0/R(KI) )/R(KI) )/(2.d0*q0*q0)


DO 20 I = 1, KI
   BR = BBB*R(I)
   SGC(I) = 4.D0*PAI*R(I)*R(I)  &
&         * ( AAA*SIN(BR)/BR + cb*cos(q0*R(I))**2 )
20 CONTINUE

sgc0  = AAA + cb
sgcd0 = -AAA*BBB*BBB/3.d0 - 2.d0*cb*q0*q0


!---------------
!      PNL11 = CST11*SGC(KI+3)/R(KI+3)/R(KI+3)
!     &      + CST12*SGC(KI+2)/R(KI+2)/R(KI+2)
!     &      + CST13*SGC(KI+1)/R(KI+1)/R(KI+1)
!     &      + CST15*SGC(KI-1)/R(KI-1)/R(KI-1)
!     &      + CST16*SGC(KI-2)/R(KI-2)/R(KI-2)
!     &      + CST17*SGC(KI-3)/R(KI-3)/R(KI-3)
!
!      PNL21 = CST21*SGC(KI+3)/R(KI+3)/R(KI+3)
!     &      + CST22*SGC(KI+2)/R(KI+2)/R(KI+2)
!     &      + CST23*SGC(KI+1)/R(KI+1)/R(KI+1)
!     &      + CST24*SGC(KI  )/R(KI  )/R(KI  )
!     &      + CST25*SGC(KI-1)/R(KI-1)/R(KI-1)
!     &      + CST26*SGC(KI-2)/R(KI-2)/R(KI-2)
!     &      + CST27*SGC(KI-3)/R(KI-3)/R(KI-3)
!      PNL2c  = ( PNL21 - PNL11 )/ ( 4.d0*PAI*R(KI)*R(KI) )
!      write(*,*) '2nd der.:', PNL2, pnl2c


RETURN
END




SUBROUTINE calrad( nfile, myid, nodes,  &
& SG, QE, WORK, radius, R, DX, MESH )
!-----------------------------------------------------------------------
!                                                          1999/04/12
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
integer :: nfile(*), myid, nodes
DIMENSION  SG(MESH)
DIMENSION  WORK(MESH), R(MESH)
DIMENSION  DA(3), DB(3)

QE99 = 0.99d0*QE

S2  =  LOG( SG(2) / SG(1) ) / DX
DO 100 I = 1, 2
   WORK(I) =  SG(I) * R(I) / ( S2 + 1.D0 )
   DB(I) =  DX * SG(I) / 3.D0
   DA(I) =  DB(I) * R(I)
100 CONTINUE
DO 110 I = 3, MESH
   DB(3) = DX * SG(I) / 3.D0
   DA(3) = DB(3) * R(I)
   WORK(I) =  WORK(I-2) + DA(3) + 4.D0 * DA(2) + DA(1)
   if( WORK(I).ge.QE99 ) then
       radius = R(I)
       return
   end if
   DO 110 L = 1, 2
      DA(L) =  DA(L+1)
      DB(L) =  DB(L+1)
110 CONTINUE

!--- error trap
call fstop( nfile, myid, nodes,  &
& 'something wrong in subroutine calrad.' )


END




subroutine lotmpp( nfile, myid, nodes,  &
& ntype, zv, rctflc, llking,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc, mx1 )
!-----------------------------------------------------------------------
!    cutoff lengths and tables for local pseudopotentials
!-----------------------------------------------------------------------
use outfile
use psvariables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
real*8,  dimension(ntype) :: zv
real*8 :: rctflc
logical, dimension(ntype) :: llking
real*8, dimension(0:mx1,ntype) :: tablc, tablca
real*8, dimension(ntype)       :: dltlc, rmxlc
real*8, dimension(0:mx1,ntype) :: tbflc, tbflca
real*8, dimension(ntype)       :: dltflc, rmxflc

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 1
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip


typedo: do it = 1, ntype
!=======================================================================

MESH =  meshi(it)
RMAX =  rmaxi(it)
DX   =  dxi(it)
DO I = 1, MESH
   R(I)   =  rdial(I,it)
end do

!--- determine cutoff length for local potenital Vloc
vlcmax = 0.d0
do i = 1, mesh
   if( SGCAE(i,it) > vlcmax ) then
       vlcmax = SGCAE(i,it)
       iatmax = i
   end if
end do
if( vlcmax < 1.d-30 ) then
  !--- should be H atom
    therct = 4.d0
else
    i = iatmax
    do
       if( log10(SGCAE(i,it)) < -8.d0 .or. i == mesh ) exit
       i = i + 1
    end do
    therct = r(i) + 2.d0
end if

!if(loutfile(1)) write(nfile(1),*) 'cutoff length for Vloc:', it, therct
!if(loutfile(2)) write(nfile(2),*) 'cutoff length for Vloc:', it, therct

rctflc = max( rctflc, therct )


dltlc(it) = min( therct, RMAX )/dble(mx1)
rmxlc(it) = min( therct, RMAX ) - 5*dltlc(it)
dltflc(it) = dltlc(it)
rmxflc(it) = rmxlc(it)
!--- correction to enforce Vloc=2*z/r in large r
irrec2 = mesh
DO I = 1, MESH
   if( R(I).ge.rmxlc(it) ) then
       irrec2 = i
       exit
   end if
end do
irrec1 = mesh
DO I = 1, MESH
   if( R(I).ge.rmxlc(it)-2.d0 ) then
       irrec1 = i
       exit
   end if
end do
irrec1 = min( irrec1, irrec2-10 )
irrec1 = max( irrec1, 1 )
cfact  = 18.d0/((r(irrec2)-r(irrec1))**2)
!      write(*,*) irrec1,irrec2,cfact
DO I = irrec1, MESH
   if( cfact*(R(I)-r(irrec1))**2 .lt. 40 ) then
       cfilter = exp(-cfact*(R(I)-r(irrec1))**2)
       diffvl  = PLOCAL(i,it) + 2.d0*zv(it)/R(I)
       PLOCAL(i,it) = - 2.d0*zv(it)/R(I) + diffvl*cfilter
     else
       PLOCAL(i,it) = - 2.d0*zv(it)/R(I)
   end if
!      write(*,*) r(i),- 2.d0*zv(it)/R(I) + diffvl,PLOCAL(i,it)
!      write(*,*) r(i), diffvl, diffvl*cfilter
end do


typeif: if( llking(it) ) then
!-----------------------------------------------------------------------
LS     = 5
N0     = 30
NOVERH = LS*2
NOVER  = NOVERH*2

!-----------------------------------------------------------------------
!     interpolation  by  spline
!      tablc(0,it)   = PLOCAL(1,it)

do ir = 1, mx1
   xrr = dltlc(it)*dble(ir)
   if( xrr.gt.r(1) ) then
       irrec = ir
       exit
   end if
end do

LL   = LS
INI  = 1
ILSM = 1
do ir = irrec, mx1
   xrr = dltlc(it)*dble(ir)
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
          YS(Nm) = PLOCAL(LP,it)
670        CONTINUE
       CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,  &
&                   KQM, MM, MM1, NDV )
       INI  = ILS - NOVER
   end if

   CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
   CALL SPLINV( LL, Nm, xrr, ppfsv, QS, AS, WS, KQM, MM, NDV, 1 )
   tablc(ir,it) = ppsv
   tbflc(ir,it) = ppfsv/xrr

end do
!      tbflc(0,it)   = 2.d0*tbflc(1,it) - tbflc(2,it)
dkata = 0.d0
do ir = irrec, irrec+4
   dkata = dkata + tbflc(ir,it)
end do
dkata  = dkata/5.d0

cfa = dkata/2.d0
cfc = tablc(irrec,it) - cfa*(dltlc(it)*irrec)**2
do ir = 0, irrec - 1
   xrr = dltlc(it)*dble(ir)
   tablc(ir,it)   = cfa*xrr*xrr + cfc
   tbflc(ir,it)   = dkata
end do

!-----------------------------------------------------------------------
   do ir = 0, mx1 - 2
      tablca(ir,it) = 2.D0*( tablc(ir,it) + tablc(ir+2,it)  &
&                                    - 2.d0*tablc(ir+1,it) )
      tbflca(ir,it) = 2.D0*( tbflc(ir,it) + tbflc(ir+2,it)  &
&                                    - 2.d0*tbflc(ir+1,it) )
   end do

!=======================================================================
end if typeif
end do typedo


return
end




subroutine set_recnrmex( nfile, myid, nodes,  &
& recnrmex, recnrm, nplw5, ldouble_grid_recip, dgalpha )
!-----------------------------------------------------------------------
!     set recnrmex
!        = 4*pi/recnrm in the usual method
!      or
!        = 4*pi/recnrm * ( 1 - exp(-recnrm/(4*a*a)) )
!          in the double-grid method
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw5
real*8,  dimension(0:nplw5) :: recnrmex
real*8,  dimension(0:nplw5) :: recnrm
logical :: ldouble_grid_recip
real*8  :: dgalpha

!-----declare local variables
real*8  :: pi, pi4, dgal2, c1b4a, eak1
integer :: ig


pi  = acos(-1.d0)
pi4 = 4.d0*pi

if( .not.ldouble_grid_recip ) then

    ig = 0
    recnrmex(ig) = 0.d0
    do ig = 1, nplw5
       recnrmex(ig) = pi4/recnrm(ig)
    end do

  else

    dgal2 = dgalpha*dgalpha
    c1b4a = 0.25d0/dgal2

    ig = 0
    recnrmex(ig) = pi/dgal2
    do ig = 1, nplw5
       eak1 = c1b4a*recnrm(ig)
       if( eak1 < 60.d0 ) then
           recnrmex(ig) = pi4/recnrm(ig)*( 1.d0 - exp(-eak1) )
         else
           recnrmex(ig) = pi4/recnrm(ig)
       end if
    end do

end if


return
end




module for_lotmpp_g
!-----------------------------------------------------------------------
! type declaration of variables in lotmpp_g
!-----------------------------------------------------------------------
implicit none

integer, parameter :: mxqab = 1000
real*8,  allocatable, dimension(:) :: rhocr, rhocra
save

end module




subroutine lotmpp_g( nfile, myid, nodes,  &
& ntype, zv, llking, dltlc, rmxlc,  &
& recnrm, recnrmex, nplw5, vlocli, xitgrd, rr, nvlcl, lstress,  &
& lvshape )
!-----------------------------------------------------------------------
!     Tables for local pseudopotentials in reciprocal space
!-----------------------------------------------------------------------
use outfile
use psvariables
use lopp_variables
use for_lotmpp_g
use planewave_decomp_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
real*8,  dimension(ntype) :: zv
logical, dimension(ntype) :: llking
real*8,  dimension(ntype) :: dltlc, rmxlc
integer :: nplw5
real*8,  dimension(0:nplw5) :: recnrm
real*8,  dimension(0:nplw5) :: recnrmex
integer :: nvlcl
real*8,  dimension(0:nvlcl)  :: vlocli ,xitgrd, rr
logical :: lstress
logical :: lvshape

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 1
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip

!-----declare local variables
integer :: status


!-----if not allocated, allocate arrays
if( .not.allocated(rhocr) ) then
    !------allocate memory
    allocate( rhocr(0:mxqab), rhocra(0:mxqab), stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in lotmpp_g' )

    alocmem = 8.d0 * (mxqab+1) * 2
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d0),' B, allocated (lotmpp_g)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d0),' B, allocated (lotmpp_g)'
end if


LS = 5
N0 = 30
NOVERH = LS*2
NOVER  = NOVERH*2

do it = 1, ntype
if( .not.llking(it) ) then
!=======================================================================

MESH =  meshi(it)
RMAX =  rmaxi(it)
DX   =  dxi(it)
DO I = 1, MESH
   R(I)   =  rdial(I,it)
end do

!-----------------------------------------------------------------------
!     set cutoff length
rhpsmx = rmxlc(it)
xrdel  = rhpsmx/dble(nvlcl-1)
!-----------------------------------------------------------------------
!     interpolation  by  spline
!      xitgrd(0)     = PLOCAL(1,it)
!      xitgrd(nvlcl) = PLOCAL(ix,it)

do ir = 1, nvlcl
   xrr = xrdel*dble(ir)
   if( xrr.gt.r(1) ) then
       irrec = ir
       exit
   end if
end do

LL   = LS
INI  = 1
ILSM = 1
do ir = irrec, nvlcl
   xrr = xrdel*dble(ir)
   IF( xrr.GT.r(ILSM) ) THEN
       do
          ILS  = MIN( INI + N0 - 1, MESH )
          IF( ILS-INI+1.LT.LL+2 ) LL = ILS - INI - 1
          IF( ILS.EQ.MESH ) THEN
              ILSM = ILS
            ELSE
              ILSM = ILS - NOVERH
          end if
          IF( r(ILSM) >= xrr ) exit
          INI = INI + NOVERH
       end do

       Nm = 0
       DO 670 lp = INI, ILS
          Nm = Nm + 1
          XS(Nm) = r(lp)
          YS(Nm) = PLOCAL(lp,it)
670        CONTINUE
       CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,  &
&                   KQM, MM, MM1, NDV )
       INI  = ILS - NOVER
   end if

   CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
   xitgrd(ir) = ppsv

end do

!---for small r
do ir = 0, irrec - 1
   xitgrd(ir) = PLOCAL(1,it)
end do

do ir = 0, nvlcl
   rr(ir) = xrdel*dble(ir)
   vlocli(ir) = ( xitgrd(ir)*rr(ir) + 2.d0*zv(it) )*rr(ir)
end do

   small=1.0d-5
!-----------------------------------------------------------------------
!  create table
   g1maxq = sqrt(recnrm(nplw5))
   dlqab = g1maxq*(1.d0+10.d0/dble(mxqab))/dble(mxqab)
   pi4 = 4.d0*pi

   do ik = 0, mxqab
      rhocr(ik) = 0.d0
   end do

   !-----distribute to nodes
   call divnod( mxqab+1, nodes, myid, ntnod1, ntnod2, ntnod )

!         do ik = 0, mxqab
   do ik = ntnod1-1, ntnod2-1
      ak1 = dlqab * dble(ik)
      if( ak1.le.small ) then
!               -----------------------------------------
          do ir = 0, nvlcl
             xitgrd(ir) = vlocli(ir)
!---                           BSL0( ak1*rr(ir) ) = 1
          end do
          call INTGB3( nvlcl, xrdel, xitgrd, CL )
!               -----------------------------------------
       else

          akdr = ak1*xrdel
          sindr = sin(akdr)
          cosdr = cos(akdr)
           akrr = 0.d0
          singr = 0.d0
          cosgr = 1.d0
              ir = 0
              xitgrd(ir) = vlocli(ir)
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl0 = singr/akrr
                 xitgrd(ir) = vlocli(ir) * qbesl0
              end do
          call INTGB3( nvlcl, xrdel, xitgrd, CL )

      end if

      rhocr(ik) = CL

   end do
   !-----unify rhocr
   call gdsum(rhocr,mxqab+1,rhocra)

   do ik = 0, mxqab - 2
      rhocra(ik) = 2.D0*( rhocr(ik) + rhocr(ik+2)  &
&                     - 2.d0*rhocr(ik+1) )
   end do
   do ik = mxqab - 1, mxqab
      rhocra(ik) = 0.d0
   end do

!--- set svco
!   do ig = 0, nplw5
   do igg = 1, nplw5nod
      ig = igg + nplw5nod1 - 2
      ak1 = sqrt(recnrm(ig))
      m = ak1/(2.d0*dlqab)
      m = 2*m
      d = 0.5d0*( ak1/dlqab - dble(m) )
      rhps = d*( (d-1.d0)*rhocra(m) + rhocr(m+2)  &
&                   - rhocr(m) ) + rhocr(m)

      svco(igg,it) = pi4*rhps - 2.d0*zv(it)*recnrmex(ig)

!            if( ak1.lt.small ) then
!                pslo =  pi4*rhps
!            else
!                pslo =  pi4*( rhps - 2.d0*zv(it)/recnrm(ig) )
!            end if
!            svco(ig,it) = pslo
   end do
!-----------------------------------------------------------------------

   stressif: if( lstress ) then
!-----------------------------------------------------------------------
! for stress calculation
   do ir = 1, nvlcl
      xr = rr(ir)
      vlocli(ir) = xr*vlocli(ir)
   end do

   do ik = 0, mxqab
      rhocr(ik) = 0.d0
   end do

   !-----distribute to nodes
   call divnod( mxqab+1, nodes, myid, ntnod1, ntnod2, ntnod )

!         do ik = 0, mxqab
   do ik = ntnod1-1, ntnod2-1
      ak1 = dlqab * dble(ik)
      if( ak1.le.small ) then
!               -----------------------------------------
          CL = 0.d0
!               -----------------------------------------
       else

          akdr = ak1*xrdel
          sindr = sin(akdr)
          cosdr = cos(akdr)
           akrr = 0.d0
          singr = 0.d0
          cosgr = 1.d0
              ir = 0
              xitgrd(ir) = 0.d0
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 akrrr= 1.d0/akrr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl1 = ( singr*akrrr - cosgr )*akrrr
                 xitgrd(ir) = vlocli(ir) * qbesl1
              end do
          call INTGB3( nvlcl, xrdel, xitgrd, CL )

      end if

      rhocr(ik) = CL

   end do
   !-----unify rhocr
   call gdsum(rhocr,mxqab+1,rhocra)

   do ik = 0, mxqab - 2
      rhocra(ik) = 2.D0*( rhocr(ik) + rhocr(ik+2)  &
&                     - 2.d0*rhocr(ik+1) )
   end do
   do ik = mxqab - 1, mxqab
      rhocra(ik) = 0.d0
   end do

!--- set svcop
   if( myid_pw == 0 ) then
       svcop(1,it) = 0.d0
       ig1 = 2
     else
       ig1 = 1
   end if
!   do ig = 1, nplw5
   do igg = ig1, nplw5nod
      ig = igg + nplw5nod1 - 2
      ak1 = sqrt(recnrm(ig))
      m = ak1/(2.d0*dlqab)
      m = 2*m
      d = 0.5d0*( ak1/dlqab - dble(m) )
      rhps = d*( (d-1.d0)*rhocra(m) + rhocr(m+2)  &
&                   - rhocr(m) ) + rhocr(m)

      if( ak1.lt.small ) then
          pslo =  pi4*rhps
      else
          pslo =  pi4*( rhps - 4.d0*zv(it)/(recnrm(ig)*ak1) )
      end if
      svcop(igg,it) = pslo/ak1
   end do
!-----------------------------------------------------------------------
   end if stressif

!=======================================================================
end if
end do


if( .not.lvshape ) then
    !------deallocate memory
    deallocate( rhocr, rhocra, stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in lotmpp_g' )

    alocmem = 8.d0 * (mxqab+1) * 2
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d0),' B, deallocated (lotmpp_g)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d0),' B, deallocated (lotmpp_g)'
end if


return
end




module for_nlkbpp_g
!-----------------------------------------------------------------------
! type declaration of allocatable variables in nlkbpp_g
!-----------------------------------------------------------------------
implicit none

real*8,  allocatable, dimension(:,:,:) :: ylmr, ylmi
real*8,  allocatable, dimension(:) :: waveli

end module




subroutine nlkbpp_g( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, lkbppi, lking,  &
& ms, aname, mxl, nplw, gx, gy, gz, lplwkmx,  &
& vlocli ,xitgrd, rr, nvlcl, ierror )
!-----------------------------------------------------------------------
!---  Kleinmann and Bylander with Troullier and Martins  ---
!     non-local pseudopotential elements
!
!     The non-local terms are calculated by PP-table and Y_lm.
!-----------------------------------------------------------------------
use outfile
use psvariables
use nlpp_variables
use for_nlkbpp_g
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
integer, dimension(*) :: nplw
real*8,  dimension(0:lplwkmx,*) :: gx, gy, gz
character(1), dimension(7)   :: ms
character(2), dimension(103) :: aname

integer :: nvlcl
real*8, dimension(0:nvlcl)  :: vlocli ,xitgrd, rr

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 1
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip

!-----declare local variables
integer :: status
character(50) :: fname
character(1)  :: dummy
integer :: jpl
complex*16 :: ci
character(1), dimension(0:9) :: num =  &
&                    (/ '0','1','2','3','4','5','6','7','8','9' /)
real*8  :: gvecmax = 0.d0
save gvecmax


!-----if not allocated, allocate arrays
if( .not.allocated(ylmr) ) then

    !------allocate memory for local variables
    allocate( ylmr(3,-mxl:mxl,0:mxl), ylmi(3,-mxl:mxl,0:mxl),  &
& waveli(0:nvlcl), stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in nlkbpp_g' )

    alocmem =  &
&  8.d0 * ( ( 3*(2*mxl+1)*(mxl+1) )*2 + nvlcl + 1 )
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d0),' B, allocated (nlkbpp_g)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d0),' B, allocated (nlkbpp_g)'

end if


small=1.0d-5
sqrtwo = sqrt(2.d0)


ikng1 = 1
ikng2 = 1

ikng = 1
g1max(ikng) = 0.d0
!      do ikng = ikng1, ikng2
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------------------------------
!     set Spherical harmonics

   lxxx(ikng) = 0
   do it = 1, ntype
   if( lkbppi(it) .and. .not.lking(it) ) then

       lmxx = 0
       do l=0, lmax(it)
          if( lchk(l,it) .and. l.ne.lclno(it) ) then
              do j = MJTAB1(l,it), MJTAB2(l,it)
                 lmxx = lmxx + 2*l + 1
              end do
          end if
       end do
       lmx(it) = lmxx

       lxxx(ikng) = max( lxxx(ikng), lmax(it) )

   end if
   end do

do k = 1, nknod
do ig = 0, nplw(k)
!         g1x = bzk(1,k) + gx(ig,k)
!         g1y = bzk(2,k) + gy(ig,k)
!         g1z = bzk(3,k) + gz(ig,k)
   g1x = gx(ig,k)
   g1y = gy(ig,k)
   g1z = gz(ig,k)
   g1xy = g1x*g1x + g1y*g1y
   g1   = g1xy + g1z*g1z
   g1xy = sqrt( g1xy )
   g1   = sqrt( g1 )
   g1max(ikng) = max( g1max(ikng), g1 )
   if( g1.gt.small ) then
       g1r  = 1.d0/g1
       cos1 = g1z/g1
       sin1 = g1xy/g1
       if( g1xy.gt.small ) then
           cos2 = g1x/g1xy
           sin2 = g1y/g1xy
         else
           cos2 = 1.d0
           sin2 = 0.d0
       end if
     else
       g1r  = 0.d0
       cos1 = 1.d0
       sin1 = 0.d0
       cos2 = 1.d0
       sin2 = 0.d0
   end if
!--------------------------------------
!     Spherical harmonics and its derivatives
   call sphrb2( lxxx(ikng), cos1, sin1, cos2, sin2, ylmr, ylmi,  &
&               MXL )

   if( lyreal ) then
       lmxx = 0
       esx = cos1*cos2
       esy = cos1*sin2
       esz = -sin1
       efx = -sin2
       efy =  cos2
       efz = 0.d0
       do l = 0, lxxx(ikng)
          ljdb2 = l/2
          if( mod(ljdb2,2) .eq. 0 ) then
              qlfugo = 1.d0
            else
              qlfugo = -1.d0
          end if
          qrtwo = qlfugo*sqrtwo
!                -------------------------------------------------------
              mj = 1
              lmm = mj - 1
              j   = lmxx + mj
              qlmr(ig,j,k,ikng) = qlfugo*ylmr(1,lmm,l)
              qlmi(ig,j,k,ikng) = 0.d0
!            --- for stress calculation
              consr = qlfugo*ylmr(2,lmm,l)*g1r
              confr = qlfugo*ylmr(3,lmm,l)*g1r
              qlmdgr(ig,j,1,k) = consr*esx + confr*efx
              qlmdgi(ig,j,1,k) = 0.d0
              qlmdgr(ig,j,2,k) = consr*esy + confr*efy
              qlmdgi(ig,j,2,k) = 0.d0
              qlmdgr(ig,j,3,k) = consr*esz + confr*efz
              qlmdgi(ig,j,3,k) = 0.d0
           do mj = 2, l + 1
              lmm = mj - 1
              j   = lmxx + mj
              qlmr(ig,j,k,ikng) = qrtwo*ylmr(1,lmm,l)
              qlmi(ig,j,k,ikng) = qrtwo*ylmi(1,lmm,l)
!            --- for stress calculation
              consr = qrtwo*ylmr(2,lmm,l)*g1r
              consi = qrtwo*ylmi(2,lmm,l)*g1r
              confr = qrtwo*ylmr(3,lmm,l)*g1r
              confi = qrtwo*ylmi(3,lmm,l)*g1r
              qlmdgr(ig,j,1,k) = consr*esx + confr*efx
              qlmdgi(ig,j,1,k) = consi*esx + confi*efx
              qlmdgr(ig,j,2,k) = consr*esy + confr*efy
              qlmdgi(ig,j,2,k) = consi*esy + confi*efy
              qlmdgr(ig,j,3,k) = consr*esz + confr*efz
              qlmdgi(ig,j,3,k) = consi*esz + confi*efz
           end do
!                -------------------------------------------------------
           lmxx = lmxx + l + 1
       end do
     else
       lmxx = 0
       ci   = ( 1.d0, 0.d0 )
       do l = 0, lxxx(ikng)
           cir = dble( ci )
           cii = dimag( ci )
!                ----------------------------------------------------
           do mj = 1, l + 1
              lmm = mj - 1
              j   = lmxx + mj
             qlmr(ig,j,k,ikng)=ylmr(1,lmm,l)*cir+ylmi(1,lmm,l)*cii
             qlmi(ig,j,k,ikng)=ylmr(1,lmm,l)*cii-ylmi(1,lmm,l)*cir
           end do
!                ----------------------------------------------------
           lmxx = lmxx + l + 1
       ci   = ci*( 0.d0, 1.d0 )
       end do
   end if

end do
end do
!-----------------------------------------------------------------------

if( g1max(ikng) <= gvecmax ) return
gvecmax = g1max(ikng)*1.1d0

ifpp = 1
if( lrela ) ifpp = ifpp + 3

call allocate_unit_number( jpl )

!-----zero clear
do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then
    do l  = 1, lnref
    do ik = 0, mxqlm
       tbqlm(ik,l,it) = 0.d0
       tbqld(ik,l,it) = 0.d0
    end do
    end do
end if
end do
!-----set counter for parallel calculation
igcount = 0
!-----------------------------------------------------------------------
LS = 5
N0 = 30
NOVERH = LS*2
NOVER  = NOVERH*2
pi4 = 4.d0*pi

typedo: do it = 1, ntype
typeif: if( lkbppi(it) .and. .not.lking(it) ) then
!=======================================================================

MESH =  meshi(it)
RMAX =  rmaxi(it)
DX   =  dxi(it)
DO I = 1, MESH
   R(I)   =  rdial(I,it)
end do

lmxx = 0
laxx = 0
ldo: do l = 0, lmax(it)
if( lchk(l,it) ) then
jdo: do j = MJTAB1(l,it), MJTAB2(l,it)
   ix = 0
   do I = MESH, 1, -1
      if( abs(PSUORG(I,j,it)).gt.1.0d-10 ) then
          ix = I + 1
          exit
      end if
   end do
   if( ix == 0 ) cycle ldo
!--------------------------------------------
!     set index
   do mj = 1, 2*l + 1
      lmm = mj - ( l + 1 )
      itolm(lmxx+mj,it) = l*(l+1)/2 + abs(lmm) + 1
!cc            itoll(lmxx+mj,it) = l
      litoll(lmxx+mj,it) = mod(l,2) .eq. 0
      if(        lmm.ge.0 ) then
                 itolmf(lmxx+mj,it) = 1
        else if( mod(l+lmm,2).eq.0 ) then
                 itolmf(lmxx+mj,it) = 2
        else
                 itolmf(lmxx+mj,it) = 3
      end if
      itola(lmxx+mj,it) = laxx + 1
   end do
!       --- cf. Vanderbilt pp. ---
   lpmx(it)  = 0
   lptmx(it) = lmx(it)
   do lp = 1, lptmx(it)
      iptv1(lp,it) = l
      iptv2(lp,it) = l
   end do
!--------------------------------------------

!         xrdel = DX*dble(ix-1)/dble(nvlcl)
!         CR   =  RMAX * exp( -DX*( MESH - 1 ) )
!         D    =  EXP( xrdel )
!         do ir = 0, nvlcl
!            rr(ir) = CR
!            CR = CR * D
!         end do
   xrdel = r(ix)/dble(nvlcl)
   do ir = 0, nvlcl
      rr(ir) = xrdel*dble(ir)
   end do
   rdels = rr(nvlcl)/dble(mm-1)
   rdels = rdels*( 1.d0 - 1.d-10 )
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
!     interpolation of pseudo-wave-function  by  spline
!      waveli(0)     = PNL(1)
!      waveli(nvlcl) = PNL(ix)
!
!      LL   = LS
!      INI  = 1
!      ILSM = 1
!      xrr = 0.d0
!      do ir = 1, nvlcl - 1
!         xrr = xrr + xrdel
!         IF( xrr.GT.DX*dble(ILSM-1) ) THEN
!             do
!                ILS  = MIN( INI + N0 - 1, MESH )
!                IF( ILS-INI+1.LT.LL+2 ) LL = ILS - INI - 1
!                IF( ILS.EQ.MESH ) THEN
!                    ILSM = ILS
!                  ELSE
!                    ILSM = ILS - NOVERH
!                end if
!                IF( DX*dble(ILSM-1) >= xrr ) exit
!                INI = INI + NOVERH
!             end do
!
!             Nm = 0
!             do LP = INI, ILS
!                Nm = Nm + 1
!                XS(Nm) = DX*dble(LP-1)
!                YS(Nm) = PNL(LP)
!             end do
!             CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,
!     &                    KQM, MM, MM1, NDV )
!             INI  = ILS - NOVER
!         end if
!
!         CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
!         waveli(ir) = ppsv
!      end do
!-----------------------------------------------------------------------
!     interpolation of pseudo-wave-function  by  spline
   i = 0
   n = 1
   xs(n) = 0.d0
!         ys(n) = PNL(1)
   ys(n) = 0.d0
   rbs  = rdels
   do 
      i = i + 1
      if( r(i).le.rbs ) cycle
      n = n + 1
      rbs = rbs + rdels
      xs(n) = r(i)
      ys(n) = PNL(i)
      if( n.ge.mm ) exit
   end do
!        --------------------------------------
   LL   = LS
   CALL SPLINS( LL, n, XS, YS, QS, BS, AS, WS, IP,  &
&               KQM, MM, MM1, NDV )

   do ir = 0, nvlcl
      xr  = rr(ir)
      call splinv( ll, n, xr, sgc, qs, as, ws, kqm, mm, ndv, 0 )
      waveli(ir) = sgc
   end do
!-----------------------------------------------------------------------
!     interpolation of pseudopotential  by  spline
!      vlocli(0)     = PSUORG(1,l,it)
!      vlocli(nvlcl) = PSUORG(ix,l,it)
!
!      LL   = LS
!      INI  = 1
!      ILSM = 1
!      xrr = 0.d0
!      do ir = 1, nvlcl - 1
!         xrr = xrr + xrdel
!         IF( xrr.GT.DX*dble(ILSM-1) ) THEN
!             do
!             ILS  = MIN( INI + N0 - 1, MESH )
!             IF( ILS-INI+1.LT.LL+2 ) LL = ILS - INI - 1
!             IF( ILS.EQ.MESH ) THEN
!                 ILSM = ILS
!               ELSE
!                 ILSM = ILS - NOVERH
!             end if
!             IF( DX*dble(ILSM-1) >= xrr ) exit
!                 INI = INI + NOVERH
!             end do
!
!             Nm = 0
!             do LP = INI, ILS
!                Nm = Nm + 1
!                XS(Nm) = DX*dble(LP-1)
!                YS(Nm) = PSUORG(LP,l,it)
!             end do
!             CALL SPLINS( LL, Nm, XS, YS, QS, BS, AS, WS, IP,
!     &                    KQM, MM, MM1, NDV )
!             INI  = ILS - NOVER
!         end if
!
!         CALL SPLINV( LL, Nm, xrr, ppsv, QS, AS, WS, KQM, MM, NDV, 0 )
!         vlocli(ir) = ppsv
!      end do
!-----------------------------------------------------------------------
!     interpolation of pseudopotential by  spline
   i = 0
   n = 1
   xs(n) = 0.d0
   ys(n) = PSUORG(1,j,it)
   rbs  = rdels
   do
      i = i + 1
      if( r(i).le.rbs ) cycle
      n = n + 1
      rbs = rbs + rdels
      xs(n) = r(i)
      ys(n) = PSUORG(i,j,it)
      if( n.ge.mm ) exit
   end do
!        --------------------------------------
   LL   = LS
   CALL SPLINS( LL, n, XS, YS, QS, BS, AS, WS, IP,  &
&               KQM, MM, MM1, NDV )

   do ir = 0, nvlcl
      xr  = rr(ir)
      call splinv( ll, n, xr, sgc, qs, as, ws, kqm, mm, ndv, 0 )
      vlocli(ir) = sgc
   end do
!-----------------------------------------------------------------------
   do ir = 0, nvlcl
      xitgrd(ir) = vlocli(ir)*waveli(ir)*waveli(ir)
   end do
   call INTGB3( nvlcl, xrdel, xitgrd, CL )
   clalpt = 1.d0/CL
   do mj = 1, 2*l + 1
      clalp(lmxx+mj,it) = clalpt * factl(j,it)
   end do

   do ir = 0, nvlcl
      vlocli(ir) = vlocli(ir)*waveli(ir)*rr(ir)
   end do
!-----------------------------------------------------------------------
!   create table
   lmm = laxx + 1
   dlqlm(it) = gvecmax*(1.d0+10.d0/dble(mxqlm))/dble(mxqlm)

igcount = igcount + 1

igcountif: if( mod(igcount,nodes) == myid ) then

   ak1 = 0.d0
!cc         ak12 = 0.d0
   do ik = 0, mxqlm
!cc            ak1 = sqrt(ak12)
      if( ak1.le.small ) then
!           -----------------------------------------
      if( l.eq.0 ) then
          do ir = 0, nvlcl
             xitgrd(ir) = vlocli(ir)
!---                             BSL0( ak1*rr(ir) ) = 1
          end do
          call INTGB3( nvlcl, xrdel, xitgrd, CL )
        else
          CL = 0.d0
      end if
!           -----------------------------------------
      else

       akdr = ak1*xrdel
      sindr = sin(akdr)
      cosdr = cos(akdr)
       akrr = 0.d0
      singr = 0.d0
      cosgr = 1.d0
      if( l.eq.0 ) then
          ir = 0
          xitgrd(ir) = vlocli(ir)
          do ir = 1, nvlcl
             akrr = akrr + akdr
             sint = singr*cosdr + cosgr*sindr
             cost = cosgr*cosdr - singr*sindr
             singr = sint
             cosgr = cost
             qbesl0 = singr/akrr
             xitgrd(ir) = vlocli(ir) * qbesl0
          end do
      else if( l.eq.1 ) then
          ir = 0
          xitgrd(ir) = 0.d0
          do ir = 1, nvlcl
             akrr = akrr + akdr
             akrrr= 1.d0/akrr
             sint = singr*cosdr + cosgr*sindr
             cost = cosgr*cosdr - singr*sindr
             singr = sint
             cosgr = cost
             qbesl1 = ( singr*akrrr - cosgr )*akrrr
             xitgrd(ir) = vlocli(ir) * qbesl1
          end do
      else if( l.eq.2 ) then
          ir = 0
          xitgrd(ir) = 0.d0
          do ir = 1, nvlcl
             akrr = akrr + akdr
             akrrr= 1.d0/akrr
             sint = singr*cosdr + cosgr*sindr
             cost = cosgr*cosdr - singr*sindr
             singr = sint
             cosgr = cost
             qbesl0 = singr*akrrr
             qbesl1 = ( qbesl0 - cosgr )*akrrr
             qbesl2 = 3.D0*qbesl1*akrrr - qbesl0
             xitgrd(ir) = vlocli(ir) * qbesl2
          end do
      else if( l.eq.3 ) then
          ir = 0
          xitgrd(ir) = 0.d0
          do ir = 1, nvlcl
             akrr = akrr + akdr
             akrrr= 1.d0/akrr
             sint = singr*cosdr + cosgr*sindr
             cost = cosgr*cosdr - singr*sindr
             singr = sint
             cosgr = cost
             qbesl0 = singr*akrrr
             qbesl1 = ( qbesl0 - cosgr )*akrrr
             qbesl2 = 3.D0*qbesl1*akrrr - qbesl0
             qbesl3 = 5.D0*qbesl2*akrrr - qbesl1
             xitgrd(ir) = vlocli(ir) * qbesl3
          end do
      end if
      call INTGB3( nvlcl, xrdel, xitgrd, CL )

      end if
      CL = CL*pi4

      tbqlm(ik,lmm,it) = CL

!cc            ak12 = ak12 + dlqlm(it)
      ak1 = ak1 + dlqlm(it)
   end do
!         do ik = 0, mxqlm - 2
!            tbqlma(ik,lmm,it) = 2.D0*( tbqlm(ik,lmm,it)
!     &                               + tbqlm(ik+2,lmm,it)
!     &                               - 2.d0*tbqlm(ik+1,lmm,it) )
!         end do

end if igcountif

!-----------------------------------------------------------------------
!   create table for stress calculation
   do ir = 0, nvlcl
      vlocli(ir) = vlocli(ir)*rr(ir)
   end do
   dlqld(it) = gvecmax*(1.d0+10.d0/dble(mxqlm))/dble(mxqlm)

igcount = igcount + 1

igcount2if: if( mod(igcount,nodes) == myid ) then

   ak1 = 0.d0
!cc         ak12 = 0.d0
   do ik = 0, mxqlm
!cc            ak1 = sqrt(ak12)
!cc            if( ikng.eq.1 ) then
      if( ak1.le.small ) then
!               -----------------------------------------
          if( l.eq.1 ) then
              do ir = 0, nvlcl
                 xitgrd(ir) = vlocli(ir) / 3.d0
!---                                 d BSL1( x )/d x = 1/3
              end do
              call INTGB3( nvlcl, xrdel, xitgrd, CL )
            else
              CL = 0.d0
          end if
!               -----------------------------------------
        else

          akdr = ak1*xrdel
          sindr = sin(akdr)
          cosdr = cos(akdr)
           akrr = 0.d0
          singr = 0.d0
          cosgr = 1.d0
          if( l.eq.0 ) then
              ir = 0
              xitgrd(ir) = 0.d0
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 akrrr= 1.d0/akrr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl1 = ( singr*akrrr - cosgr )*akrrr
                 dbessl = - qbesl1
                 xitgrd(ir) = vlocli(ir) * dbessl
              end do
          else if( l.eq.1 ) then
              ir = 0
              xitgrd(ir) = vlocli(ir) / 3.d0
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 akrrr= 1.d0/akrr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl0 = singr*akrrr
                 qbesl1 = ( qbesl0 - cosgr )*akrrr
                 qbesl2 = 3.D0*qbesl1*akrrr - qbesl0
                 dbessl = - qbesl2 + qbesl1*akrrr
                 xitgrd(ir) = vlocli(ir) * dbessl
              end do
          else if( l.eq.2 ) then
              ir = 0
              xitgrd(ir) = 0.d0
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 akrrr= 1.d0/akrr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl0 = singr*akrrr
                 qbesl1 = ( qbesl0 - cosgr )*akrrr
                 qbesl2 = 3.D0*qbesl1*akrrr - qbesl0
                 qbesl3 = 5.D0*qbesl2*akrrr - qbesl1
                 dbessl = - qbesl3 + 2.d0*qbesl2*akrrr
                 xitgrd(ir) = vlocli(ir) * dbessl
              end do
          else if( l.eq.3 ) then
              ir = 0
              xitgrd(ir) = 0.d0
              do ir = 1, nvlcl
                 akrr = akrr + akdr
                 akrrr= 1.d0/akrr
                 sint = singr*cosdr + cosgr*sindr
                 cost = cosgr*cosdr - singr*sindr
                 singr = sint
                 cosgr = cost
                 qbesl0 = singr*akrrr
                 qbesl1 = ( qbesl0 - cosgr )*akrrr
                 qbesl2 = 3.D0*qbesl1*akrrr - qbesl0
                 qbesl3 = 5.D0*qbesl2*akrrr - qbesl1
                 qbesl4 = 7.D0*qbesl3*akrrr - qbesl2
                 dbessl = - qbesl4 + 3.d0*qbesl3*akrrr
                 xitgrd(ir) = vlocli(ir) * dbessl
              end do
          end if
          call INTGB3( nvlcl, xrdel, xitgrd, CL )

      end if
!cc              else
!ccc          --- by spline ---
!cc              call splinv(ll, n, ak1, CL, qs, as, ws, kqm, mm, ndv, 0)
!cc            end if
      CL = CL*pi4

      tbqld(ik,lmm,it) = CL

!cc            ak12 = ak12 + dlqld(it)
      ak1 = ak1 + dlqld(it)
   end do
!         do ik = 0, mxqlm - 2
!            tbqlda(ik,lmm,it) = 2.D0*( tbqld(ik,lmm,it)
!     &                               + tbqld(ik+2,lmm,it)
!     &                               - 2.d0*tbqld(ik+1,lmm,it) )
!         end do

end if igcount2if
!heck
!         do ik = 0, mxqlm
!            write(90+l,'(4e14.6)') dlqld(it)*dble(ik),
!     &             tbqld(ik,lmm,it)
!         end do

   lmxx = lmxx + 2*l + 1
   laxx = laxx + 1
!-----------------------------------------------------------------------
end do jdo
end if
end do ldo
if( lmx(it).ne.lmxx ) then
    ierror = 101
    write(*,*) '*** not match lmx (2) in nlkbpp_g'
end if
lax(it) = laxx
!=======================================================================
end if typeif
end do typedo


do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then
    !-----unify tbqlm by global sum
    call gdsum( tbqlm(0,1,it), (mxqlm+1)*lnref, tbqlma(0,1,it) )

    !-----unify tbqld by global sum
    call gdsum( tbqld(0,1,it), (mxqlm+1)*lnref, tbqlda(0,1,it) )

    do lmm = 1, lnref
       do ik = 0, mxqlm - 2
          tbqlma(ik,lmm,it) = 2.D0*( tbqlm(ik,lmm,it)  &
&                              + tbqlm(ik+2,lmm,it)  &
&                              - 2.d0*tbqlm(ik+1,lmm,it) )
       end do
       do ik = 0, mxqlm - 2
          tbqlda(ik,lmm,it) = 2.D0*( tbqld(ik,lmm,it)  &
&                              + tbqld(ik+2,lmm,it)  &
&                              - 2.d0*tbqld(ik+1,lmm,it) )
       end do
       do ik = mxqlm - 1, mxqlm
          tbqlma(ik,lmm,it) = 0.d0
          tbqlda(ik,lmm,it) = 0.d0
       end do
    end do
end if
end do

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      end do
!heck
!+      ctg1 = cputime()
!+      do it = 1, ntype
!+      do ia = natom1(it), natom2(it)
!+      do l  = 1, 1
!+         do ig=-nplw, nplw
!+            m = dkgnrm(ig)/(2.d0*dlqlm(it))
!+            m = 2*m
!+            d = 0.5d0*( dkgnrm(ig)/dlqlm(it) - dble(m) )
!+            ycosd(ig) = d*( (d-1.d0)*tbqlma(m,l,it) + tbqlm(m+2,l,it)
!+     &                            - tbqlm(m,l,it) ) + tbqlm(m,l,it)
!+c                write(90+l,'(3e20.11)') ak1, CL, value-CL
!+         end do
!+      end do
!+      end do
!+      end do
!+      ctg2 = cputime()
!+      write(*,*) 'cpu-time for cal. qlm :', ctg2 - ctg1
!heck



!      !------deallocate memory for local variables
!      deallocate( ylmr, ylmi, waveli, stat=status )
!
!      !------error trap
!      status = abs(status)
!      call gimax(status)
!      if( status /= 0 ) call fstop( nfile, myid, nodes,
!     & 'memory deallocation error in nlkbpp_g' )
!
!      alocmem =
!     &   8.d0 * ( ( 3*(2*mxl+1)*(mxl+1) )*2 + nvlcl + 1 )
!      if( myid == 0 ) then
!        write(nfile(1),*) nint(alocmem/1d0),' B, deallocated (nlkbpp_g)'
!        write(nfile(2),*) nint(alocmem/1d0),' B, deallocated (nlkbpp_g)'
!      end if

call deallocate_unit_number( jpl )


return
end




subroutine sphrcl( l2, c1, s1, c2, s2, ylmr, ylmi, mxl )
!-----------------------------------------------------------------------
!     Spherical harmonics
!-----------------------------------------------------------------------
implicit real*8(a-h,o-z)
dimension ylmr(-mxl:mxl,0:mxl), ylmi(-mxl:mxl,0:mxl)
parameter( pi = 3.141592653589793d0,  &
&          pi00 = 2.820947917738781d-01,  &
&          pi10 = 4.886025119029199d-01,  &
&          pi11 = 3.454941494713355d-01,  &
&          pi20 = 3.153915652525200d-01,  &
&          pi21 = 7.725484040463791d-01,  &
&          pi22 = 3.862742020231896d-01,  &
&          pi30 = 3.731763325901154d-01,  &
&          pi31 = 3.231801841141507d-01,  &
&          pi32 = 1.021985476433282d0  ,  &
&          pi33 = 4.172238236327841d-01 )

!      ylmr(0,0) = sqrt( 1.d0/(4.d0*pi) )
ylmr(0,0) = pi00
ylmi(0,0) = 0.d0
if( l2.le.0 ) return

!      ylmr(1, 0) = sqrt( 3.d0/(4.d0*pi) )*c1
ylmr(0, 1) = pi10*c1
ylmi(0, 1) = 0.d0
!      a1 = sqrt( 3.d0/(8.d0*pi) )*s1
a1 = pi11*s1
ylmr(1, 1) = a1*c2
ylmi(1, 1) = a1*s2
ylmr(-1,1) = -ylmr(1,1)
ylmi(-1,1) =  ylmi(1,1)
if( l2.le.1 ) return

c12 = c1*c1
s12 = s1*s1
c22p = c2*c2 - s2*s2
s22p = s2*c2 + c2*s2
!      ylmr(0, 2) = sqrt( 5.d0/(16.d0*pi) )*(3.d0*c12 - 1.d0)
ylmr(0, 2) = pi20*(3.d0*c12 - 1.d0)
ylmi(0, 2) = 0.d0
!      a1 = sqrt( 15.d0/(8.d0*pi) )*s1*c1
a1 = pi21*s1*c1
ylmr(1, 2) = a1*c2
ylmi(1, 2) = a1*s2
ylmr(-1,2) = -ylmr(1,2)
ylmi(-1,2) =  ylmi(1,2)
!      a1 = sqrt( 15.d0/(32.d0*pi) )*s12
a1 = pi22*s12
ylmr(2, 2) = a1*c22p
ylmi(2, 2) = a1*s22p
ylmr(-2,2) =  ylmr(2,2)
ylmi(-2,2) = -ylmi(2,2)
if( l2.le.2 ) return

c13 = c12*c1
s13 = s12*s1
c23p = c22p*c2 - s22p*s2
s23p = s22p*c2 + c22p*s2
!      ylmr(0, 3) = sqrt( 7.d0/(16.d0*pi) )*c1*(5.d0*c12 - 3.d0)
ylmr(0, 3) = pi30*c1*(5.d0*c12 - 3.d0)
ylmi(0, 3) = 0.d0
!      a1 = sqrt( 21.d0/(64.d0*pi) )*s1*(5.d0*c12 - 1.d0)
a1 = pi31*s1*(5.d0*c12 - 1.d0)
ylmr(1, 3) = a1*c2
ylmi(1, 3) = a1*s2
ylmr(-1,3) = -ylmr(1,3)
ylmi(-1,3) =  ylmi(1,3)
!      a1 = sqrt( 105.d0/(32.d0*pi) )*s12*c1
a1 = pi32*s12*c1
ylmr(2, 3) = a1*c22p
ylmi(2, 3) = a1*s22p
ylmr(-2,3) =  ylmr(2,3)
ylmi(-2,3) = -ylmi(2,3)
!      a1 = sqrt( 35.d0/(64.d0*pi) )*s13
a1 = pi33*s13
ylmr(3, 3) = a1*c23p
ylmi(3, 3) = a1*s23p
ylmr(-3,3) = -ylmr(3,3)
ylmi(-3,3) =  ylmi(3,3)
!cc      if( l2.le.3 ) return

return
end




DOUBLE PRECISION FUNCTION sbessl( q, l )
!-----------------------------------------------------------------------
!   spherical Bessel functions
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
logical lflag
save lflag
data lflag / .true. /

if( lflag ) call cnstts
lflag = .false.

if( l.eq.0 ) then
    sbessl = BSL0( q )
else if( l.eq.1 ) then
    sbessl = BSL1( q )
else if( l.eq.2 ) then
    sbessl = BSL2( q )
else if( l.eq.3 ) then
    sbessl = BSL3( q )
else if( l.eq.-1 ) then
    if( q.lt.1.d-15 ) then
        sbessl = 0.d0
      else
        sbessl = cos(q)/q
    end if
else if( l.gt.3 ) then
    if( q.lt.1.d-15 ) then
        sbessl = 0.d0
      else
    sbm2 = BSL2( q )
    sbm1 = BSL3( q )
    do i = 4, l
       sbessl = (2.d0*dble(i)-1.d0)*sbm1/q - sbm2
       sbm2 = sbm1
       sbm1 = sbessl
    end do
    end if
else
   write(*,*) ' error in sbessl :',l
   sbessl = 0.d0
!-----Finalize the parallel environment
!         call end_parallel(ierr)
!         stop
end if

return
end




module param_for_Bessel
!-----------------------------------------------------------------------
! type declaration of coefficients for Bessel functions
!-----------------------------------------------------------------------
implicit none

real*8 :: FJ02, FJ04, FJ06, FJ08, FJ010, FJ012
real*8 :: FJ11, FJ13, FJ15, FJ17, FJ19,  FJ111
real*8 :: FJ22, FJ24, FJ26, FJ28, FJ210, FJ212
real*8 :: FJ33, FJ35, FJ37, FJ39,  FJ311
save

end module




SUBROUTINE CNSTTS
!-----------------------------------------------------------------------
!     some constants
!-----------------------------------------------------------------------
use param_for_Bessel
PARAMETER( JISP = 15 )
IMPLICIT REAL*8 (A-H,O-Z)
dimension    CIF(0:JISP)

CN = 0.0
FN = 1.D0
CIF(0) = 1.D0
DO 10 L = 1, JISP
   CN = CN + 1.D0
   FN = FN*CN
   CIF(L) = 1.D0/FN
10 CONTINUE
!-----------------------------------------------------------------------
!             coefficients for Bessel functions

FCC2  = - CIF(2)
FCC4  =   CIF(4)
FCC6  = - CIF(6)
FCC8  =   CIF(8)
FCC10 = - CIF(10)
FCC12 =   CIF(12)

FJ02  = - CIF(3)
FJ04  =   CIF(5)
FJ06  = - CIF(7)
FJ08  =   CIF(9)
FJ010 = - CIF(11)
FJ012 =   CIF(13)

FJ11  =   2.D0*CIF(3)
FJ13  = - 4.D0*CIF(5)
FJ15  =   6.D0*CIF(7)
FJ17  = - 8.D0*CIF(9)
FJ19  =  10.D0*CIF(11)
FJ111 = -12.D0*CIF(13)

FJ22  =   2.D0* 4.D0*CIF(5)
FJ24  = - 4.D0* 6.D0*CIF(7)
FJ26  =   6.D0* 8.D0*CIF(9)
FJ28  = - 8.D0*10.D0*CIF(11)
FJ210 =  10.D0*12.D0*CIF(13)
!     FJ212 = -12.D0*14.D0*CIF(15)

FJ33  =   2.D0* 4.D0* 6.D0*CIF(7)
FJ35  = - 4.D0* 6.D0* 8.D0*CIF(9)
FJ37  =   6.D0* 8.D0*10.D0*CIF(11)
FJ39  = - 8.D0*10.D0*12.D0*CIF(13)
FJ311 =  10.D0*12.D0*14.D0*CIF(15)

RETURN
END


DOUBLE PRECISION FUNCTION BSL0( Q )
use param_for_Bessel
IMPLICIT REAL*8 (A-H,O-Z)

IF( ABS(Q).LT.7.D-01 ) THEN
    XF12 = Q*Q
    BSL0 = 1.0D0 + XF12*( FJ02 + XF12*( FJ04 + XF12*( FJ06  &
&                + XF12*( FJ08 + XF12*FJ010 ) ) ) )
  ELSE
    BSL0 = SIN(Q)/Q
end if

RETURN
END


DOUBLE PRECISION FUNCTION BSL1( Q )
use param_for_Bessel
IMPLICIT REAL*8 (A-H,O-Z)

IF( ABS(Q).LT.7.D-01 ) THEN
    XF12 = Q*Q
    BSL1 = Q*( FJ11 + XF12*( FJ13 + XF12*( FJ15  &
&                      + XF12*( FJ17 + XF12*FJ19 ) ) ) )
  ELSE
    BSL1 = ( SIN(Q)/Q - COS(Q) )/Q
end if

RETURN
END


DOUBLE PRECISION FUNCTION BSL2( Q )
use param_for_Bessel
IMPLICIT REAL*8 (A-H,O-Z)

IF( ABS(Q).LT.7.D-01 ) THEN
    XF12 = Q*Q
    BSL2 = XF12*( FJ22 + XF12*( FJ24 + XF12*( FJ26  &
&                      + XF12*( FJ28 + XF12*FJ210 ) ) ) )
  ELSE
    BS0 = SIN(Q)/Q
    BS1 = ( BS0 - COS(Q) )/Q
    BSL2 = 3.D0*BS1/Q - BS0
end if

RETURN
END


DOUBLE PRECISION FUNCTION BSL3( Q )
use param_for_Bessel
IMPLICIT REAL*8 (A-H,O-Z)

IF( ABS(Q).LT.7.D-01 ) THEN
    XF12 = Q*Q
    BSL3 = Q*XF12*( FJ33 + XF12*( FJ35 + XF12*( FJ37  &
&                           + XF12*( FJ39 + XF12*FJ311 ) ) ) )
  ELSE
    BS0 = SIN(Q)/Q
    BS1 = ( BS0 - COS(Q) )/Q
    BS2 = 3.D0*BS1/Q - BS0
    BSL3 = 5.D0*BS2/Q - BS1
end if

RETURN
END




SUBROUTINE INTGBP( MSH, DX, R, BL, CL )
!-----------------------------------------------------------------------
!                                                          1992/12/08
!  Simpson's quadrature by using log-scale mesh
!
!    output :  CL ... integration of BL
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION      BL(*)
DIMENSION      R(*)
DATA DZERO / 1.0D-30 /

MESH = MSH
IF( MOD(MSH,2).EQ.0 ) MESH = MESH - 1
CL  =  0.D0
DO 10 I = 2, MESH - 1, 2
       CL  =  CL + 2.D0*BL(I)*R(I) + BL(I+1)*R(I+1)
10 CONTINUE
IX = MESH
CL  =  BL(1)*R(1) - BL(IX)*R(IX) + 2.D0*CL
CL  =  CL*DX/3.D0
IF( MOD(MSH,2).EQ.0 ) THEN
    CL  =  CL + ( BL(MSH)*R(MSH) + BL(IX)*R(IX) ) * DX * .5D0
end if
IF( ABS(BL(1)).GT.DZERO .and. BL(1)*BL(2).GT.0.d0 ) THEN
    S2 = LOG( BL(2)/BL(1) )/DX + 1.0D0
    CL = CL + BL(1)*R(1)/S2
end if


RETURN
END




SUBROUTINE INTGB3( MSH, DX, BL, CL )
!-----------------------------------------------------------------------
!                                                          1995/09/18
!  Simpson's quadrature by using normal mesh
!                                ^^^^^^
!    output :  CL ... integration of BL
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION      BL(0:*)
DATA DZERO / 1.0D-30 /

MESH = MSH
IF( MOD(MSH,2).ne.0 ) MESH = MESH - 1
CL  =  0.D0
DO 10 I = 1, MESH - 1, 2
       CL  =  CL + 2.D0*BL(I) + BL(I+1)
10 CONTINUE
IX = MESH
CL  =  BL(0) - BL(IX) + 2.D0*CL
CL  =  CL*DX/3.D0
IF( MOD(MSH,2).ne.0 ) THEN
    CL  =  CL + ( BL(MSH) + BL(IX) ) * DX * .5D0
end if

RETURN
END




subroutine denini( nfile, myid, nodes,  &
& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
& lorthrhmbc, rho, natom, mx1,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!     initial density by atomic charge density
!-----------------------------------------------------------------------
use outfile
use psvariables
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust
integer, dimension(3) :: nd1v
integer :: ntype, natom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3) :: rdel
logical :: lorthrhmbc
real*8  :: rho(mshnod)
integer :: mx1
integer :: mshnod
integer :: mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
integer :: mshx1, mshy1, mshz1, mshx, mshy, mshz

!---declare local variables
real*8,  allocatable, dimension(:) :: tbfnl, tbfnla
real*8  :: b(3,3)
!--- for bulk calculation ----------------------------------------
real*8  :: qlx(8) = (/ 0.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0 /)
real*8  :: qly(8) = (/ 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0 /)
real*8  :: qlz(8) = (/ 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0 /)
!------- dimensions for spline interpolation ---------------------------
integer :: ls = 5
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 0
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip
!-----------------------------------------------------------------------
real*8  :: alocmem
integer :: m1, it, i, n, ll, ir
integer :: ix, iy, iz, m, loop8, im
real*8  :: pai4, rdels, rbs, xsmax, dltchg, xrr, r2, sgc, rvol2
real*8  :: x1, y1, z1, xx, yy, zz, d, cutmin, RK1Di, sgcsum
real*8  :: q1, q2, q3, qq1, qq2, qq3
integer :: status


!------allocate memory
allocate( tbfnl(0:mx1), tbfnla(0:mx1), stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in denini' )

alocmem = 8.d0 * ( mx1 + 1 )*2
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, allocated (denini)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, allocated (denini)'


do m1 = 1, mshnod
   rho( m1 ) = 0.d0
end do

pai4 = 4.d0*acos(-1.d0)
do it = 1, ntype
!=======================================================================

!--- prepare radial meshes ---------------------------------------------
MESH =  meshi(it)
DO I = 1, MESH
   R(I) = rdial(I,it)
end do

   do i = 1, mesh
      pnl(i) = svrhops(i,it)/( pai4*r(i)*r(i) )
   end do
   n = 1
   xs(n) = 0.d0
   ys(n) = pnl(1)
   rdels = rad99(it)/dble(mm-1)
   rbs  = rdels
   do i = 1, mesh
      if( r(i).gt.rbs ) then
          n = n + 1
          rbs = rbs + rdels
          xs(n) = r(i)
          ys(n) = pnl(i)
          if( n.ge.mm ) exit
      end if
   end do
   xsmax = xs(n)*xs(n)
   ll   = ls
   call splins( ll, n, xs, ys, qs, bs, as, ws, ip,  &
&               kqm, mm, mm1, ndv )

   dltchg = xsmax/dble(mx1-2)
   xrr = 0.d0
   do ir = 1, mx1 - 3
      xrr = xrr + dltchg
      r2 = sqrt(xrr)
      call splinv( ll, n, r2, sgc, qs, as, ws, kqm, mm, ndv, 0 )
      tbfnl(ir) = sgc
   end do
   tbfnl(0)     = ys(1)
   tbfnl(mx1-2) = ys(n)
   tbfnl(mx1-1) = 0.d0
   tbfnl(mx1)   = 0.d0
   do ir = 0, mx1 - 2
      tbfnla(ir) = 2.D0*( tbfnl(ir)  &
&                                + tbfnl(ir+2)  &
&                                - 2.d0*tbfnl(ir+1) )
   end do

   if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
   do i = nhk1(it), nhk2(it)
      x1 = ratm(1,i)
      y1 = ratm(2,i)
      z1 = ratm(3,i)
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         xx = dble(ix-1)*rdel(1)  -  x1
         yy = dble(iy-1)*rdel(2)  -  y1
         zz = dble(iz-1)*rdel(3)  -  z1

         r2 = xx*xx + yy*yy + zz*zz
         if( r2.lt.xsmax ) then
!                   r2 = sqrt( r2 )
!              call splinv( ll, n, r2, sgc, qs, as, ws, kqm, mm, ndv, 0 )
             m = r2/(2.d0*dltchg)
             m = 2*m
             d = 0.5d0*( r2/dltchg - dble(m) )
             sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)

             rho( m1 ) = rho( m1 ) + sgc
         end if
      end do
   end do
   else
!--- for bulk calculation ----------------------------------------------
!        maximum distance for half image
   CALL RCIPRL( hcell, b, rvol2 )
   cutmin = 1.d+10
   do i = 1, 3
      RK1Di = B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i)
      rk1di = 0.5d0/sqrt( rk1di )
      cutmin = min( cutmin, rk1di )
   end do
   if( cutmin.ge.sqrt(xsmax) ) then
       loop8 = 1
     else
       loop8 = 8
   end if
   do i = nhk1(it), nhk2(it)
      x1 = ratm(1,i)
      y1 = ratm(2,i)
      z1 = ratm(3,i)

orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
#ifdef VECTOR
      do im = 1, loop8
#endif
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         q1 = dble(ix-1)/dble(nd1v(1))  -  x1
         q2 = dble(iy-1)/dble(nd1v(2))  -  y1
         q3 = dble(iz-1)/dble(nd1v(3))  -  z1
         if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.d0,q1)
         if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.d0,q2)
         if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.d0,q3)

         sgcsum = 0.d0
#ifndef VECTOR
         do im = 1, loop8
#endif
            xx = hcell(1,1)*( q1 - sign(qlx(im),q1) )
            yy = hcell(2,2)*( q2 - sign(qly(im),q2) )
            zz = hcell(3,3)*( q3 - sign(qlz(im),q3) )
            r2 = xx*xx + yy*yy + zz*zz
            if( r2.lt.xsmax ) then
                m = r2/(2.d0*dltchg)
                m = 2*m
                d = 0.5d0*( r2/dltchg - dble(m) )
                sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)
                sgcsum = sgcsum + sgc
            end if
#ifndef VECTOR
         end do
#endif
         rho( m1 ) = rho( m1 ) + sgcsum
      end do
#ifdef VECTOR
      end do
#endif

else orthoif

!------non-orthorhombic super cell -------------------------------
#ifdef VECTOR
      do im = 1, loop8
#endif
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         q1 = dble(ix-1)/dble(nd1v(1))  -  x1
         q2 = dble(iy-1)/dble(nd1v(2))  -  y1
         q3 = dble(iz-1)/dble(nd1v(3))  -  z1
         if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.d0,q1)
         if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.d0,q2)
         if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.d0,q3)

         sgcsum = 0.d0
#ifndef VECTOR
         do im = 1, loop8
#endif
            qq1 = q1 - sign(qlx(im),q1)
            qq2 = q2 - sign(qly(im),q2)
            qq3 = q3 - sign(qlz(im),q3)
            xx = hcell(1,1)*qq1 + hcell(1,2)*qq2 + hcell(1,3)*qq3
            yy = hcell(2,1)*qq1 + hcell(2,2)*qq2 + hcell(2,3)*qq3
            zz = hcell(3,1)*qq1 + hcell(3,2)*qq2 + hcell(3,3)*qq3
            r2 = xx*xx + yy*yy + zz*zz
            if( r2.lt.xsmax ) then
                m = r2/(2.d0*dltchg)
                m = 2*m
                d = 0.5d0*( r2/dltchg - dble(m) )
                sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)
                sgcsum = sgcsum + sgc
            end if
#ifndef VECTOR
         end do
#endif
         rho( m1 ) = rho( m1 ) + sgcsum
      end do
#ifdef VECTOR
      end do
#endif

end if orthoif

   end do
!--- end of setting initial charge density -----------------------------
   end if
end do


!------deallocate memory
deallocate( tbfnl, tbfnla, stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in denini' )

alocmem = 8.d0 * ( mx1 + 1 )*2
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, deallocated (denini)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, deallocated (denini)'


return
end




subroutine spindenini1( nfile, myid, nodes,  &
& rhoud, rho, mshnod, diffud, nel )
!-----------------------------------------------------------------------
!     initial spin density by atomic charge density
!-----------------------------------------------------------------------
use outfile
use param_inispin
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8  :: rho(mshnod), rhoud(mshnod)
real*8  :: diffud
integer :: nel

!---declare local variables
integer :: m1


if( inispin /= 1 ) return

do m1 = 1, mshnod
   rhoud(m1) = rho(m1) * diffud/dble(nel)
end do


return
end




subroutine spindenini( nfile, myid, nodes,  &
& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
& lorthrhmbc, rhoud, natom, mx1,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, lspinini )
!-----------------------------------------------------------------------
!     initial spin density by atomic charge density
!-----------------------------------------------------------------------
use outfile
use param_inispin
use psvariables
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust
integer, dimension(3) :: nd1v
integer :: ntype, natom
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3) :: rdel
logical :: lorthrhmbc
integer :: mx1
integer :: mshnod
integer :: mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
integer :: mshx1, mshy1, mshz1, mshx, mshy, mshz
real*8  :: rhoud(mshnod)
logical :: lspinini

!---declare local variables
real*8,  allocatable, dimension(:) :: tbfnl, tbfnla
logical, allocatable, dimension(:) :: latommagnemom
real*8,  allocatable, dimension(:) :: atommagnemom
real*8  :: b(3,3)
!--- for bulk calculation ----------------------------------------
real*8  :: qlx(8) = (/ 0.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0 /)
real*8  :: qly(8) = (/ 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0 /)
real*8  :: qlz(8) = (/ 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0 /)
!------- dimensions for spline interpolation ---------------------------
integer :: ls = 5
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 0
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip
!-----------------------------------------------------------------------
real*8  :: alocmem
integer :: m1, it, i, n, ll, ir, itt, ia
integer :: ix, iy, iz, m, loop8, im
real*8  :: sum
real*8  :: pai4, rdels, rbs, xsmax, dltchg, xrr, r2, sgc, rvol2
real*8  :: x1, y1, z1, xx, yy, zz, d, cutmin, RK1Di, sgcsum
real*8  :: q1, q2, q3, qq1, qq2, qq3
real*8 :: seed = 25.d0
integer :: status


lspinini = .false.
if( inispin == 1 ) return

lspinini = .true.


!------allocate memory
allocate( tbfnl(0:mx1), tbfnla(0:mx1),  &
& latommagnemom(natom), atommagnemom(natom),  &
& stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in spindenini' )

alocmem = 8.d0 * ( mx1 + 1 )*2 + natom + 8.d0*natom
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, allocated (spindenini)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, allocated (spindenini)'


!   inispin =  1: uniformly polarized (default)
!              2: ferromagnetic alignment
!              3: antiferromagnetic random alignment
!              4: specify spin polarization
if( myid == 0 ) then
    latommagnemom(1:natom) = .false.
     atommagnemom(1:natom) = 0.d0
    if( inispin == 2 ) then
        do itt = 1, nspinat
           it = iatspin(itt)
           if( it <= ntype ) then
               latommagnemom(nhk1(it):nhk2(it)) = .true.
                atommagnemom(nhk1(it):nhk2(it)) = atmmagne(itt)
           else
!               if(loutfile(1)) write(nfile(1),'(a,3i8)') '*** warning in spindenini(1): it > ntype',  &
!&                                                        itt, it, ntype
!               if(loutfile(2)) write(nfile(2),'(a,3i8)') '*** warning in spindenini(1): it > ntype',  &
!&                                                        itt, it, ntype
           end if
        end do
    else if( inispin == 3 ) then
        do itt = 1, nspinat
           it = iatspin(itt)
           if( it <= ntype ) then
               n  = nhk2(it) - nhk1(it) + 1
               do ia = 1, n/2
               do
                 CALL rnd00( xrr, seed )
                 i = int(xrr*n) + 1
                 if( i > n ) i = n
                 i = nhk1(it) + i - 1
                 if( .not.latommagnemom(i) ) then
                     latommagnemom(i) = .true.
                      atommagnemom(i) = atmmagne(itt)
                     exit
                 end if
               end do
               end do
               do i = nhk1(it), nhk2(it)
                  if( .not.latommagnemom(i) ) then
                      latommagnemom(i) = .true.
                       atommagnemom(i) = -atmmagne(itt)
                  end if
               end do
           else
!               if(loutfile(1)) write(nfile(1),'(a,3i8)') '*** warning in spindenini(2): it > ntype',  &
!&                                                        itt, it, ntype
!               if(loutfile(2)) write(nfile(2),'(a,3i8)') '*** warning in spindenini(2): it > ntype',  &
!&                                                        itt, it, ntype
           end if
        end do
    else !if( inispin == 4 ) then
        do ia = 1, nspinat
           i = iatspin(ia)
           if( i <= natom ) then
               latommagnemom(i) = .true.
                atommagnemom(i) = atmmagne(ia)
           else
!               if(loutfile(1)) write(nfile(1),'(a,3i8)') '*** warning in spindenini(3): i > natom',  &
!&                                                        ia, i, natom
!               if(loutfile(2)) write(nfile(2),'(a,3i8)') '*** warning in spindenini(3): i > natom',  &
!&                                                        ia, i, natom
           end if
        end do
    end if
end if
!---unify latommagnemom & atommagnemom by broadcasting
call lbcast(latommagnemom,natom,0)
call dbcast( atommagnemom,natom,0)
call unify_by_lbcast( latommagnemom, natom )
call unify_by_dbcast(  atommagnemom, natom )


do m1 = 1, mshnod
   rhoud( m1 ) = 0.d0
end do

pai4 = 4.d0*acos(-1.d0)
do it = 1, ntype
!=======================================================================

!--- prepare radial meshes ---------------------------------------------
MESH =  meshi(it)
DO I = 1, MESH
   R(I) = rdial(I,it)
end do

   DX   =  dxi(it)
   CALL INTGBP( mesh, DX, R, svrhops(1,it), sum )

   atommagnemom(nhk1(it):nhk2(it)) = atommagnemom(nhk1(it):nhk2(it))/sum

   do i = 1, mesh
      pnl(i) = svrhops(i,it)/( pai4*r(i)*r(i) )
   end do
   n = 1
   xs(n) = 0.d0
   ys(n) = pnl(1)
   rdels = rad99(it)/dble(mm-1)
   rbs  = rdels
   do i = 1, mesh
      if( r(i).gt.rbs ) then
          n = n + 1
          rbs = rbs + rdels
          xs(n) = r(i)
          ys(n) = pnl(i)
          if( n.ge.mm ) exit
      end if
   end do
   xsmax = xs(n)*xs(n)
   ll   = ls
   call splins( ll, n, xs, ys, qs, bs, as, ws, ip,  &
&               kqm, mm, mm1, ndv )

   dltchg = xsmax/dble(mx1-2)
   xrr = 0.d0
   do ir = 1, mx1 - 3
      xrr = xrr + dltchg
      r2 = sqrt(xrr)
      call splinv( ll, n, r2, sgc, qs, as, ws, kqm, mm, ndv, 0 )
      tbfnl(ir) = sgc
   end do
   tbfnl(0)     = ys(1)
   tbfnl(mx1-2) = ys(n)
   tbfnl(mx1-1) = 0.d0
   tbfnl(mx1)   = 0.d0
   do ir = 0, mx1 - 2
      tbfnla(ir) = 2.D0*( tbfnl(ir)  &
&                                + tbfnl(ir+2)  &
&                                - 2.d0*tbfnl(ir+1) )
   end do

   if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
   do i = nhk1(it), nhk2(it)
      if( .not.latommagnemom(i) ) cycle
      x1 = ratm(1,i)
      y1 = ratm(2,i)
      z1 = ratm(3,i)
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         xx = dble(ix-1)*rdel(1)  -  x1
         yy = dble(iy-1)*rdel(2)  -  y1
         zz = dble(iz-1)*rdel(3)  -  z1

         r2 = xx*xx + yy*yy + zz*zz
         if( r2.lt.xsmax ) then
!                   r2 = sqrt( r2 )
!              call splinv( ll, n, r2, sgc, qs, as, ws, kqm, mm, ndv, 0 )
             m = r2/(2.d0*dltchg)
             m = 2*m
             d = 0.5d0*( r2/dltchg - dble(m) )
             sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)

             rhoud( m1 ) = rhoud( m1 ) + sgc * atommagnemom(i)
         end if
      end do
   end do
   else
!--- for bulk calculation ----------------------------------------------
!        maximum distance for half image
   CALL RCIPRL( hcell, b, rvol2 )
   cutmin = 1.d+10
   do i = 1, 3
      RK1Di = B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i)
      rk1di = 0.5d0/sqrt( rk1di )
      cutmin = min( cutmin, rk1di )
   end do
   if( cutmin.ge.sqrt(xsmax) ) then
       loop8 = 1
     else
       loop8 = 8
   end if
   do i = nhk1(it), nhk2(it)
      if( .not.latommagnemom(i) ) cycle
      x1 = ratm(1,i)
      y1 = ratm(2,i)
      z1 = ratm(3,i)

orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
#ifdef VECTOR
      do im = 1, loop8
#endif
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         q1 = dble(ix-1)/dble(nd1v(1))  -  x1
         q2 = dble(iy-1)/dble(nd1v(2))  -  y1
         q3 = dble(iz-1)/dble(nd1v(3))  -  z1
         if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.d0,q1)
         if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.d0,q2)
         if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.d0,q3)

         sgcsum = 0.d0
#ifndef VECTOR
         do im = 1, loop8
#endif
            xx = hcell(1,1)*( q1 - sign(qlx(im),q1) )
            yy = hcell(2,2)*( q2 - sign(qly(im),q2) )
            zz = hcell(3,3)*( q3 - sign(qlz(im),q3) )
            r2 = xx*xx + yy*yy + zz*zz
            if( r2.lt.xsmax ) then
                m = r2/(2.d0*dltchg)
                m = 2*m
                d = 0.5d0*( r2/dltchg - dble(m) )
                sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)
                sgcsum = sgcsum + sgc * atommagnemom(i)
            end if
#ifndef VECTOR
         end do
#endif
         rhoud( m1 ) = rhoud( m1 ) + sgcsum
      end do
#ifdef VECTOR
      end do
#endif

else orthoif

!------non-orthorhombic super cell -------------------------------
#ifdef VECTOR
      do im = 1, loop8
#endif
      do m1 = 1, mshnod
         ix = mshnx(m1)
         ix = mshx1 + ix - 1
         iy = mshny(m1)
         iy = mshy1 + iy - 1
         iz = mshnz(m1)
         iz = mshz1 + iz - 1
         q1 = dble(ix-1)/dble(nd1v(1))  -  x1
         q2 = dble(iy-1)/dble(nd1v(2))  -  y1
         q3 = dble(iz-1)/dble(nd1v(3))  -  z1
         if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.d0,q1)
         if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.d0,q2)
         if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.d0,q3)

         sgcsum = 0.d0
#ifndef VECTOR
         do im = 1, loop8
#endif
            qq1 = q1 - sign(qlx(im),q1)
            qq2 = q2 - sign(qly(im),q2)
            qq3 = q3 - sign(qlz(im),q3)
            xx = hcell(1,1)*qq1 + hcell(1,2)*qq2 + hcell(1,3)*qq3
            yy = hcell(2,1)*qq1 + hcell(2,2)*qq2 + hcell(2,3)*qq3
            zz = hcell(3,1)*qq1 + hcell(3,2)*qq2 + hcell(3,3)*qq3
            r2 = xx*xx + yy*yy + zz*zz
            if( r2.lt.xsmax ) then
                m = r2/(2.d0*dltchg)
                m = 2*m
                d = 0.5d0*( r2/dltchg - dble(m) )
                sgc = d*( (d-1.d0)*tbfnla(m)  &
&                      + tbfnl(m+2)  &
&                      - tbfnl(m)   ) + tbfnl(m)
                sgcsum = sgcsum + sgc * atommagnemom(i)
            end if
#ifndef VECTOR
         end do
#endif
         rhoud( m1 ) = rhoud( m1 ) + sgcsum
      end do
#ifdef VECTOR
      end do
#endif

end if orthoif

   end do
!--- end of setting initial charge density -----------------------------
   end if
end do


!------deallocate memory
deallocate( tbfnl, tbfnla, latommagnemom, atommagnemom, stat=status )

!------error trap
status = abs(status)
call gimax(status)
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in spindenini' )

alocmem = 8.d0 * ( mx1 + 1 )*2 + natom + 8.d0*natom
if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d3),' KB, deallocated (spindenini)'
if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d3),' KB, deallocated (spindenini)'


return
end




