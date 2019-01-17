



module nl_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
implicit none

integer :: lplnplw = 0
integer :: lplnplwex = 0
real*8, allocatable, dimension(:) :: xyzr
real*8, allocatable, dimension(:) :: zzzr, zzzi
real*8, allocatable, dimension(:,:) :: blmr, blmi

integer :: nbnodt
integer :: mxbdat
real*8, allocatable, dimension(:,:,:) :: slmir
real*8, allocatable, dimension(:,:,:) :: slmii

logical :: lbigmm2 = .true.
!-----available, if lbigmm2 = .true.
real*8, allocatable, dimension(:,:,:) :: eqlmr_

real*8, allocatable, dimension(:,:) :: eqlmr, eqlmi

#ifdef VECTOR
real*8, allocatable, dimension(:) :: wrk1r, wrk2r, wrk3
real*8, allocatable, dimension(:) :: wrk1i, wrk2i
real*8, allocatable, dimension(:) :: wrk3_atm
#endif

!-----only for norm-conserving pp.
real*8,  allocatable, dimension(:,:) :: vclalp


!-----for force calculation
real*8,  allocatable, dimension(:) :: fslxr, fslyr, fslzr
real*8,  allocatable, dimension(:) :: fslxi, fslyi, fslzi
#ifdef VECTOR
real*8,  allocatable, dimension(:) :: vclweg
real*8,  allocatable, dimension(:) :: dslmir
real*8,  allocatable, dimension(:) :: dslmii
real*8,  allocatable, dimension(:) :: qslmir, qslmii
#endif

!-----for stress calculation
real*8,  allocatable, dimension(:,:,:) :: blmdr, blmdi
real*8,  allocatable, dimension(:,:,:) :: eqlmdr
real*8,  allocatable, dimension(:,:,:,:) :: slmidr
real*8,  allocatable, dimension(:,:,:,:) :: slmidi

!-----for noncollinear magnetism
real*8, allocatable, dimension(:,:,:) :: slmbr
real*8, allocatable, dimension(:,:,:) :: slmbi
#ifdef VECTOR
real*8, allocatable, dimension(:) :: wrk4r, wrk5r
real*8, allocatable, dimension(:) :: wrk4i, wrk5i
#endif
real*8,  allocatable, dimension(:,:) :: ncfslxr, ncfslyr, ncfslzr
real*8,  allocatable, dimension(:,:) :: ncfslxi, ncfslyi, ncfslzi
#ifdef VECTOR
real*8,  allocatable, dimension(:) :: dslmbr
real*8,  allocatable, dimension(:) :: dslmbi
real*8,  allocatable, dimension(:,:) :: ncvclweg
#endif
real*8,  allocatable, dimension(:,:,:,:) :: slmbdr
real*8,  allocatable, dimension(:,:,:,:) :: slmbdi

save

end module




module nlkbpp_variables
!-----------------------------------------------------------------------
! type declaration of shared variables for nonlocal pp.
!-----------------------------------------------------------------------
implicit none

real*8,  allocatable, dimension(:,:,:,:,:) :: svsumr, svsumi
real*8,  allocatable, dimension(:,:,:,:,:) :: svsubr, svsubi
real*8,  allocatable, dimension(:,:,:) :: slmir_atm
real*8,  allocatable, dimension(:,:,:) :: slmii_atm
real*8,  allocatable, dimension(:,:,:) :: slmir_nod
real*8,  allocatable, dimension(:) :: bufcr2
real*8,  allocatable, dimension(:) :: ncslmir, ncslmii, ncslmbr, ncslmbi
integer :: natnodx

save

end module




subroutine slmir_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, nbnod, ntype, natom, natnod,  &
& lkbpp, lvandi, lvand, lspin, nprvmx, nspnmx, laspc )
!-----------------------------------------------------------------------
!     allocate memory for shared variables in nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nband, nbnod
integer :: ntype, natom, natnod
logical :: lkbpp, lvandi(ntype), lvand
logical :: lspin
integer :: nprvmx, nspnmx
logical :: laspc

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: natomx = -1
integer :: natomxx
save natomx


nbnodt = nbnod
!-----global maximum
call gimax(nbnodt)

if( natom /= natomx ) then

    the_mem = 0.d0
    if( allocated(slmir) ) then
        the_mem = the_mem + 8.d0 * ( size(slmir) + size(slmii) )
        deallocate( slmir, slmii, stat=status )
    end if

    if( allocated(slmbr) ) then
        the_mem = the_mem + 8.d0 * ( size(slmbr) + size(slmbi) )
        deallocate( slmbr, slmbi, stat=status )
    end if

    if( the_mem > 0.5d0 ) then
        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'slmir_alloc', .true. )
    end if


    natomx = natom
    natomxx = max( natomx, 1)
    !------allocate memory
    allocate( slmir(lmmax,natomxx,nbnodt), stat=status )

    the_mem = 8.d0 * ( size(slmir) )

    if( status == 0 ) then
        if( .not.lgamma .or. lnoncollinear ) then
            allocate( slmii(lmmax,natomxx,nbnodt), stat=status )
          else
            allocate( slmii(1,1,1), stat=status )
        end if

        the_mem = the_mem + 8.d0 * ( size(slmii) )
    end if

    if( status == 0 ) then
        !-----noncollinear magnetism
        if( lnoncollinear ) then
            allocate( slmbr(lmmax,natomxx,nbnodt), slmbi(lmmax,natomxx,nbnodt), stat=status )
            the_mem = the_mem + 8.d0 * ( size(slmbr) + size(slmbi) )
        end if
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'slmir_alloc', .true. )

end if


mxbdat = lmmax * max( natom*nbnodt, max(natnod,1)*nband )


slmir(:,:,:) = 0.d0
slmii(:,:,:) = 0.d0
if( lnoncollinear ) then
    slmbr(:,:,:) = 0.d0
    slmbi(:,:,:) = 0.d0
end if


!call slmir_k_alloc( nfile, myid, nodes,  &
!& alloc_mem, nband, nbnod, lmmax, natom, natnod, lspin )

!if( lvand ) then
!    call uspp_slmir_alloc( nfile, myid, nodes,  &
!& alloc_mem, nband, nbnod, ntype, natom, natnod,  &
!& lvandi, lspin, nprvmx, nspnmx, laspc )
!end if

if( lkbpp ) then
    call slmir_atm_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, natom, natnod )
end if


return
end subroutine




subroutine slmir_atm_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, natom, natnod )
!-----------------------------------------------------------------------
!     allocate memory for shared variables in nlkbpp.f
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
use nlkbpp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nband
integer :: natom, natnod

!-----declare local variables
integer :: nbnodt2
integer :: status
real*8  :: the_mem, the_deallocmem
integer :: natomx = -1, natomxx
integer :: natnodxx = 0
integer :: mxbdatx = 0
save natomx, natnodxx, mxbdatx


natnodx = max( natnod, 1 )

nbnodt2 = nbnodt
if( natom > 0 ) then
    do
       if( lmmax*natom*nbnodt2 > lmmax*natnodx*nband ) exit
       nbnodt2 = nbnodt2 + 1
    end do
end if

the_mem = 0d0
the_deallocmem = 0d0
status = 0

if( natom /= natomx ) then
    if( allocated(slmir_nod) ) then
        the_deallocmem = the_deallocmem + 8.d0 * ( size(slmir_nod) )
        deallocate( slmir_nod, stat=status )
    end if

    natomx = natom
    natomxx = max( natomx, 1)
    if( status == 0 ) then
        !------allocate memory
        allocate( slmir_nod(lmmax,natomxx,nbnodt2), stat=status )
        the_mem = the_mem + 8.d0 * ( size(slmir_nod) )
    end if

    if( status == 0 ) then
        slmir_nod(1:lmmax,1:natom,1:nbnodt) = 0.d0
    end if
end if


if( status == 0 ) then
if( natnodx /= natnodxx ) then
    if( allocated(slmir_atm) ) then
        the_deallocmem = the_deallocmem + 8.d0 * ( size(slmir_atm) )
        deallocate( slmir_atm, stat=status )

        if( status == 0 ) then
        if( .not.lgamma .or. ltddft .or. lnoncollinear ) then
            the_deallocmem = the_deallocmem + 8.d0 * ( size(slmii_atm) )
            deallocate( slmii_atm, stat=status )
        end if
        end if
    end if

    natnodxx = natnodx
    if( status == 0 ) then
        !------allocate memory
        allocate( slmir_atm(lmmax,natnodxx,nband), stat=status )
        the_mem = 8.d0 * ( size(slmir_atm) )
    end if

    if( status == 0 ) then
    if( .not.lgamma .or. ltddft .or. lnoncollinear ) then
        allocate( slmii_atm(lmmax,natnodxx,nband), stat=status )
        the_mem = the_mem + 8.d0 * size(slmii_atm)
    end if
    end if

    if( status == 0 ) then
        slmir_atm(1:lmmax,1:natnodx,1:nband) = 0.d0
    end if
end if
end if


if( status == 0 ) then
if( mxbdat > mxbdatx ) then
    if( allocated(bufcr2) ) then
        the_deallocmem = the_deallocmem + 8.d0 * ( size(bufcr2) )
        deallocate( bufcr2, stat=status )
    end if

    mxbdatx = mxbdat
    if( status == 0 ) then
        !------allocate memory
        allocate( bufcr2(mxbdatx), stat=status )
        the_mem = 8.d0 * ( size(bufcr2) )
    end if

    if( status == 0 ) then
        bufcr2(1:mxbdatx) = 0.d0
    end if
end if
end if


!------error trap
if( the_deallocmem > 0.5d0 ) then
    call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_deallocmem, 'slmir_atm_alloc', .true. )
end if
if( the_mem > 0.5d0 ) then
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'slmir_atm_alloc', .true. )
end if


if( .not.allocated(slmii_atm) ) allocate( slmii_atm(1,1,1), stat=status )


!-----first-order SO effects by purturbation calculation
!call slmso_atm_alloc( nfile, myid, nodes,  &
!& alloc_mem, natom, nband, natnodx )


return
end subroutine




subroutine nlkbpp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, mxl, nylmmx, nbnod, nspnmx )
!-----------------------------------------------------------------------
!     allocate memory for variables for nonlocal pp.
!-----------------------------------------------------------------------
use nlkbpp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: natom, mxl, nylmmx, nbnod
integer :: nspnmx

!-----declare local variables
integer :: nknod
logical :: lgamma
integer :: status
real*8  :: the_mem
logical :: lnoncollinear


!---get No. of k points
call get_nknod( nknod )
call get_lgamma( lgamma )

call get_lnoncollinear( lnoncollinear )


!if( .not.lgamma .or. nspnmx == 2 ) then
    !------allocate memory
    allocate( svsumr(nylmmx,natom,nbnod,nspnmx,nknod), stat=status )

    the_mem = 8.d0 * ( size(svsumr) )

    if( status == 0 ) then
        if( .not.lgamma .or. lnoncollinear ) then
            !------allocate memory
            allocate( svsumi(nylmmx,natom,nbnod,nspnmx,nknod), stat=status )
        else
            allocate( svsumi(1,1,1,1,1), stat=status )
        end if

        the_mem = the_mem + 8.d0 * ( size(svsumi) )
    end if

    if( status == 0 ) then
        if( lnoncollinear ) then
            !------allocate memory
            allocate( svsubr(nylmmx,natom,nbnod,nspnmx,nknod),  &
&                     svsubi(nylmmx,natom,nbnod,nspnmx,nknod),  &
&                     ncslmir(nylmmx), ncslmii(nylmmx),  &
&                     ncslmbr(nylmmx), ncslmbi(nylmmx),  &
& stat=status )

            the_mem = the_mem + 8.d0 * ( size(svsubr) + size(svsubi)  &
& + size(ncslmir) + size(ncslmii) + size(ncslmbr) + size(ncslmbi) )
        end if
    end if
!end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nlkbpp_variables_alloc', .true. )

svsumr(:,:,:,:,:) = 0d0
svsumi(:,:,:,:,:) = 0d0


return
end subroutine




subroutine nl_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, nplw, nplwex, nband, nbnod, nbncnt, nbndsp,  &
& ntype, nhk1, nhk2, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& nhk1_nat, nhk2_nat, lvandi, lstress, lspin, natom_alloc, natnod_alloc )
!-----------------------------------------------------------------------
!     allocate memory for shared variables in nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nplw, nplwex
integer :: nband, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: ntype
integer :: natom, natnod1, natnod2, natnod
integer, dimension(nodes) :: natcnt, natdsp
integer, dimension(ntype) :: nhk1, nhk2, nhk1_nat, nhk2_nat
logical, dimension(ntype) :: lvandi
logical :: lstress
logical :: lspin
integer :: natom_alloc, natnod_alloc

!-----declare local variables
integer :: i, j, k, l, ia, iab, ib
integer :: status
real*8  :: the_mem
logical :: licall = .true.
#ifdef VECTOR
integer :: ncnt, it, la, lb
#endif
save licall


if( licall ) then
    !------allocate memory
#ifdef VECTOR
    allocate(  &
& wrk1r(lptmax*natom_alloc), wrk2r(lptmax*natom_alloc), wrk3(lptmax*natom_alloc),  &
& wrk3_atm(lptmax*natnod_alloc),  &
& stat=status )
    the_mem =  &
& + 8.d0 * ( size(wrk1r) + size(wrk2r) + size(wrk3) + size(wrk3_atm) )
#else
    status = 0
    the_mem = 0.d0
#endif

    if( status == 0 ) then
        if( lstress ) then
            allocate( slmidr(lmmax,natom,6,nbnod), stat=status )
        else
            allocate( slmidr(1,1,1,1), stat=status )
        end if

        the_mem = the_mem + 8.d0 * ( size(slmidr) )
    end if

    if( status == 0 ) then
#ifdef VECTOR
        if( .not.lgamma .or. ltddft .or. lnoncollinear ) then
            allocate( wrk1i(lptmax*natom_alloc), wrk2i(lptmax*natom_alloc), stat=status )
            the_mem = the_mem + 8.d0 * ( size(wrk1i) + size(wrk2i) )
        end if
#endif
        if( .not.lgamma .or. lnoncollinear ) then
            if( status == 0 .and. lstress ) then
                allocate( slmidi(lmmax,natom,6,nbnod), stat=status )
                if( lnoncollinear ) then
                    allocate( slmbdr(lmmax,natom,6,nbnod), slmbdi(lmmax,natom,6,nbnod), stat=status )
                    the_mem = the_mem + 8.d0 * ( size(slmbdr) + size(slmbdi) )
                end if
            else if( status == 0 .and. .not.lstress ) then
                allocate( slmidi(1,1,1,1), stat=status )
            end if
        else
            allocate( slmidi(1,1,1,1), stat=status )
        end if
        the_mem = the_mem + 8.d0 * ( size(slmidi) )
    end if


#ifdef VECTOR
    if( lnoncollinear ) then
        allocate(  &
& wrk4r(lptmax*natom_alloc), wrk5r(lptmax*natom_alloc),  &
& wrk4i(lptmax*natom_alloc), wrk5i(lptmax*natom_alloc),  &
& stat=status )
        the_mem = the_mem + 8.d0 * ( size(wrk4r) + size(wrk5r)  &
&                                  + size(wrk4i) + size(wrk5i) )
    end if
#endif

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_variables_alloc', .true. )

    licall = .false.
end if


!-----for TDDFT
!call nl_variables_tddft_alloc( nfile, myid, nodes,  &
!& alloc_mem, lmmax, natom, nband, lspin )


#ifdef VECTOR
ncnt = 0
do it = 1, ntype
   if( lvandi(it) ) then
       do l = 1, lptmx(it)
          la = iptv1(l,it)
          lb = iptv2(l,it)
          do ia = nhk1(it), nhk2(it)
             ncnt = ncnt + 1
             wrk3(ncnt) = qablm(l,it)
          end do
       end do
   end if
end do
ncnt = 0
do it = 1, ntype
   if( lvandi(it) ) then
       do l = 1, lptmx(it)
          la = iptv1(l,it)
          lb = iptv2(l,it)
          do ia = nhk1_nat(it), nhk2_nat(it)
             ncnt = ncnt + 1
             wrk3_atm(ncnt) = qablm(l,it)
          end do
       end do
   end if
end do
#endif


if( lstress ) then
    slmidr = 0.d0
    if( lnoncollinear ) then
        slmbdr = 0.d0
        slmbdi = 0.d0
    end if
end if

if( .not.lgamma .or. lnoncollinear ) then
    if( lstress ) then
        slmidi = 0.d0
    end if
end if


return
end subroutine




subroutine nlg_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, nplw, nplwex, nband, nbnod, nbncnt, nbndsp,  &
& ntype, lvandi, lstress, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for shared variables in nlkbpp.f & nlvand.f
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nplw, nplwex
integer :: nband, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: ntype
logical, dimension(ntype) :: lvandi
logical :: lstress
real*8  :: pwscale

!-----declare local variables
integer :: i, j, k, l, ia, iab, ib
real*8  :: pkscale
integer :: status
real*8  :: the_mem
logical :: licall = .true.
save licall


if( licall ) then

    !------allocate memory
    allocate( vclalp(nylmmx,ntype),  &
#ifdef VECTOR
& dslmir(lptmax*nbnod),  &
#endif
& stat=status )

    the_mem =  &
& + 8.d0 * ( size(vclalp) )
#ifdef VECTOR
    if( allocated(dslmir) ) the_mem = the_mem + 8.d0 * size(dslmir)
#endif

    if( status == 0 ) then
    if( .not.lnoncollinear ) then
        allocate(  &
#ifdef VECTOR
& vclweg(lptmax*nbnod),  &
& fslxr(lptmax*nbnod), fslyr(lptmax*nbnod), fslzr(lptmax*nbnod),  &
#else
& fslxr(lmmax), fslyr(lmmax), fslzr(lmmax),  &
#endif
& stat=status )

#ifdef VECTOR
        if( allocated(vclweg) ) the_mem = the_mem + 8.d0 * size(vclweg)
#endif
        the_mem = the_mem  &
& + 8.d0 * ( size(fslxr) + size(fslyr) + size(fslzr) )
    else
        !-----noncollinear magnetism
        allocate(  &
#ifdef VECTOR
& ncvclweg(lptmax*nbnod,4), dslmbr(lptmax*nbnod),  &
& ncfslxr(lptmax*nbnod,2), ncfslyr(lptmax*nbnod,2), ncfslzr(lptmax*nbnod,2),  &
#else
& ncfslxr(lmmax,2), ncfslyr(lmmax,2), ncfslzr(lmmax,2),  &
#endif
& stat=status )

#ifdef VECTOR
        if( allocated(dslmbr) ) the_mem = the_mem + 8.d0 * ( size(dslmbr) + size(ncvclweg) )
#endif
        the_mem = the_mem  &
& + 8.d0 * ( size(ncfslxr) + size(ncfslyr) + size(ncfslzr) )
    end if
    end if

    if( status == 0 ) then
    if( .not.lgamma .or. ltddft .or. lnoncollinear ) then
#ifdef VECTOR
        allocate( dslmii(lptmax*nbnod), stat=status )
        the_mem = the_mem + 8.d0 * ( size(dslmii) )
#endif

        if( status == 0 ) then
        if( .not.lnoncollinear ) then
            allocate(  &
#ifdef VECTOR
& fslxi(lptmax*nbnod), fslyi(lptmax*nbnod), fslzi(lptmax*nbnod),  &
#else
& fslxi(lmmax), fslyi(lmmax), fslzi(lmmax),  &
#endif
& stat=status )

            the_mem = the_mem  &
& + 8.d0 * ( size(fslxi) + size(fslyi) + size(fslzi) )
        else
            !-----noncollinear magnetism
            allocate(  &
#ifdef VECTOR
& dslmbi(lptmax*nbnod),  &
& ncfslxi(lptmax*nbnod,2), ncfslyi(lptmax*nbnod,2), ncfslzi(lptmax*nbnod,2),  &
#else
& ncfslxi(lmmax,2), ncfslyi(lmmax,2), ncfslzi(lmmax,2),  &
#endif
& stat=status )

#ifdef VECTOR
            if( allocated(dslmbi) ) the_mem = the_mem + 8.d0 * size(dslmbi)
#endif
            the_mem = the_mem  &
& + 8.d0 * ( size(ncfslxi) + size(ncfslyi) + size(ncfslzi) )
        end if
        end if

    end if
    end if

    if( status == 0 .and. ltddft ) then
#ifdef VECTOR
        allocate(  &
& qslmir(lptmax*nbnod), qslmii(lptmax*nbnod),  &
& stat=status )
        the_mem = the_mem + 8.d0*( size(qslmir) + size(qslmii) )
#endif
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_variables_alloc', .true. )

end if


!------allocate memory
call get_pkscale( pkscale )

!-----check array sizes
if( nplw   * pkscale > lplnplw   .or.  &
&   nplwex * pkscale > lplnplwex  ) then

    !-----if already allocated, deallocate arrays
    if( allocated(zzzr) ) then
        deallocate( zzzr, zzzi, xyzr, blmr, blmi,  &
& stat=status )

        the_mem =  &
&  8.d0 * ( (lplnplw+1)*2 + lplnplwex )  &
& + 8.d0 * ( (lplnplw+1)*lmmax*2 )

        if( allocated(eqlmr) .and. status == 0 ) then
            deallocate( eqlmr, eqlmi, stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( (lplnplw+1)*lmmax*2 )
        end if

        if( allocated(blmdr) .and. status == 0 ) then
            deallocate( blmdr, blmdi, eqlmdr, stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( (lplnplw+1)*6*lmmax*4 )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_variables_alloc', .true. )
    end if

    lplnplw   = nplw   * pwscale * pkscale
    lplnplwex = nplwex * pwscale * pkscale
    allocate( zzzr(0:lplnplw), zzzi(0:lplnplw), xyzr(lplnplwex),  &
& eqlmr(0:lplnplw,lmmax), eqlmi(0:lplnplw,lmmax),  &
& blmr(0:lplnplw,lmmax), blmi(0:lplnplw,lmmax),  &
& stat=status )

    if( status == 0 .and. lstress ) allocate(  &
& blmdr(0:lplnplw,6,lmmax),  blmdi(0:lplnplw,6,lmmax),  &
& eqlmdr(2*(lplnplw+1),6,lmmax),  &
& stat=status )

    the_mem =  &
&  8.d0 * ( (lplnplw+1)*2 + lplnplwex )  &
& + 8.d0 * ( (lplnplw+1)*lmmax*4 )

    if( lstress ) the_mem = the_mem  &
& + 8.d0 * ( (lplnplw+1)*6*lmmax*4 )

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_variables_alloc', .true. )

end if


licall = .false.

return
end subroutine




subroutine nl_bigmem_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nplw, pwscale, natom_alloc )
!-----------------------------------------------------------------------
!     allocate big memory for shared variables in nlkbpp.f & nlvand.f, if possible
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: natom
integer :: nplw
real*8  :: pwscale
integer :: natom_alloc

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: lplnplwbig = 0
save lplnplwbig


!      call get_lgamma( lgamma )
if( .not.lgamma .or. lnoncollinear ) lbigmm2 = .false.

!-----check array sizes
if( nplw > lplnplwbig .and. lbigmm2 ) then

    !-----if already allocated, deallocate arrays
    if( allocated(eqlmr_) ) then

        the_mem = 8.d0 * size(eqlmr_)
        deallocate( eqlmr_, stat=status )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_bigmem_alloc', .false. )
    end if

    !------allocate memory
    lplnplwbig = nplw * pwscale
    allocate( eqlmr_(2*(lplnplwbig+1),lmmax,natom_alloc),  &
& stat=status )

    the_mem = 8.d0 * size(eqlmr_)

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_bigmem_alloc', .false. )

    lbigmm2 = status == 0


    if( lbigmm2 ) then

    if( allocated(eqlmr) ) then
        the_mem = 8.d0 * ( size(eqlmr) + size(eqlmi) )
        !------deallocate memory
        deallocate( eqlmr, eqlmi, stat=status )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nl_bigmem_alloc', .true. )
    end if

    end if

end if


return
end subroutine




subroutine vlocal_g( nfile, myid, nodes,  &
& vext, ntype, nhk, nhk1, nhk2, iatoit, iatmpt, natom, nion,  &
& llking, nplw5, nplw5ex, nplw, kfft1d, kfft2d, kfft3d, kfft0d,  &
& nga, ngb, ngc, ijkgd, thrhgr, rvol,  &
& mshnod, mftnod, mftdsp, mfd2ft, ntotfd,  &
& fft3x, fft3y, fftwork, lstress,  & !& ycos, ysin, nplwcs,
& DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3, kmax1, kmax2, kmax3 )
!-----------------------------------------------------------------------
!    local pseudopotential in reciprocal space : vext
!-----------------------------------------------------------------------
use lopp_variables
use planewave_decomp_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
integer, dimension(ntype) :: nhk
integer, dimension(ntype) :: nhk1, nhk2
integer :: nion, natom
integer, dimension(natom) :: iatoit
integer, dimension(nion)  :: iatmpt
logical, dimension(ntype) :: llking
real*8  :: rvol
integer :: nplw5, nplw
integer :: kfft1d, kfft2d, kfft3d, kfft0d
integer, dimension(0:nplw5) :: nga, ngb, ngc
integer, dimension(-nplw5:nplw5) :: ijkgd
real*8,  dimension(nplw5ex) :: thrhgr
integer :: mshnod, ntotfd
integer, dimension(*) :: mftnod, mftdsp, mfd2ft
real*8,  dimension(mshnod) :: vext
real*8,  dimension(*) :: fft3x, fft3y
complex*16,  dimension(*) :: fftwork
logical :: lstress
!real*8,  dimension(0:nplwcs,natom) :: ycos, ysin
integer :: kmax1, kmax2, kmax3
real*8,  dimension(-kmax1:kmax1,natom) :: DCQ1, DSQ1
real*8,  dimension(-kmax2:kmax2,natom) :: DCQ2, DSQ2
real*8,  dimension(-kmax3:kmax3,natom) :: DCQ3, DSQ3

!-----declare local variables
integer :: myid_, nodes_, nkd_

do ir = 1, kfft0d
   fft3x(ir) = 0.d0
   fft3y(ir) = 0.d0
end do

rvolr = 1.d0/rvol

if( myid_pw == 0 ) then
    ig  = 0
    ijk = ijkgd(ig)
    do it = 1, ntype
    if( .not.llking(it) ) then
       DCSN = dble( nhk(it) )*svco(1,it)
       fft3x(ijk) = fft3x(ijk) + DCSN*rvolr
    end if
    end do
    ig1 = 2
  else
    ig1 = 1
end if

thrhgr(1:nplw5ex) = 0.d0
do it = 1, ntype
if( .not.llking(it) ) then
   do ig = 1, 2*nplw5nod
      rhmur(ig) = 0.d0
   end do
   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)

!      do ig = 1, nplwcs
!         rhmur(2*ig-1) = rhmur(2*ig-1) + ycos(ig,ia)
!         rhmur(2*ig  ) = rhmur(2*ig  ) + ysin(ig,ia)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
      do igg = ig1, nplw5nod
         ig = igg + nplw5nod1 - 2
         K1M = nga(ig)
         K2M = ngb(ig)
         K3M = ngc(ig)
         DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                DSQ1(K1M,ia)*DSQ2(K2M,ia)
         DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                DCQ1(K1M,ia)*DSQ2(K2M,ia)

         ycosiK = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
         ysineK = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)

         rhmur(2*igg-1) = rhmur(2*igg-1) + ycosiK
         rhmur(2*igg  ) = rhmur(2*igg  ) + ysineK
      end do
   end do

!   do ig = 1, nplw5
   do igg = ig1, nplw5nod
      DGGG =  svco(igg,it)*rvolr
      phar =  rhmur(2*igg-1)*DGGG
      phai = -rhmur(2*igg  )*DGGG
      thrhgr(2*igg-1) = thrhgr(2*igg-1) + phar
      thrhgr(2*igg  ) = thrhgr(2*igg  ) + phai
   end do

end if
end do
if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid_, nodes_ )
!    !-----global sum
!    call gdsum(thrhgr,2*nplw5,rhmur)
    do igg = 2*ig1-1, 2*nplw5nod
       rhmur(igg) = thrhgr(igg)
    end do
    call alldgatherv(rhmur,2*nplw5nod,thrhgr,nplw5excnt,nplw5exdsp)

    ig = 0
    ijk = ijkgd(ig)
    call ibcast(ijk,1,0)
    call dbcast(fft3x(ijk),1,0)

    !-----set communicator
    call get_worldkd( myid_, nodes_, nkd_ )
end if
!-----global sum
call gdsum(thrhgr(3),2*nplw5,rhmur)
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  = fft3x(ijk)  + thrhgr(2*ig+1)
   fft3y(ijk)  = fft3y(ijk)  + thrhgr(2*ig+2)
end do
do ig = 1, nplw5
   ijkm = ijkgd(-ig)
   fft3x(ijkm) = fft3x(ijkm) + thrhgr(2*ig+1)
   fft3y(ijkm) = fft3y(ijkm) - thrhgr(2*ig+2)
end do

!        ..............................
!        ... vlocal(g) to vlocal(r) ...
!        ..............................
!          .................
!          ... constants ...
!          .................
nml = 1
inv = 2

call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1d, kfft2d, kfft3d, kfft0d, ierrft )

!      --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3y(ifd) = fft3x(ift)
end do

call distlc( nfile, myid, nodes,  &
& fft3x, mshnod, fft3y, mftdsp )

do i = 1, mshnod
   vext(i) = vext(i) + fft3x(i)
end do

!-----for stress calculation
if( lstress ) then
    call save_vext( nfile, myid, nodes, fft3x, mshnod )
end if


return
end




subroutine flocal_g( nfile, myid, nodes,  &
& floc, strloc, lcstress, ntype, nhk1, nhk2, iatoit, iatmpt,  &
& natom, nion, llking, nplw5, nplw5ex, nplw,  &
& nga, ngb, ngc, gx, gy, gz, rhgr, elclv, & !ycos, ysin, nplwcs,  &
& kmax1, kmax2, kmax3, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& ltimecnt )
!-----------------------------------------------------------------------
!    floc : force by local pseudopotential in reciprocal space
!-----------------------------------------------------------------------
use outfile
use lopp_variables
use planewave_decomp_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ntype
logical :: lcstress
integer, dimension(ntype) :: nhk1, nhk2
integer :: nion, natom
real*8,  dimension(3,natom) :: floc
real*8,  dimension(3,3)     :: strloc
integer, dimension(natom) :: iatoit
integer, dimension(nion)  :: iatmpt
logical, dimension(ntype) :: llking
integer :: nplw5, nplw
integer, dimension(0:nplw5) :: nga, ngb, ngc
real*8,  dimension(0:nplw5) :: gx, gy, gz
real*8,  dimension(nplw5ex) :: rhgr
real*8  :: elclv
!real*8,  dimension(0:nplwcs,natom) :: ycos, ysin
integer :: kmax1, kmax2, kmax3
real*8,  dimension(-kmax1:kmax1,natom) :: DCQ1, DSQ1
real*8,  dimension(-kmax2:kmax2,natom) :: DCQ2, DSQ2
real*8,  dimension(-kmax3:kmax3,natom) :: DCQ3, DSQ3
logical :: ltimecnt

!-----declare local variables
real*8,  dimension(6) :: bufst, bufstr
integer :: myid_, nodes_, nkd_


if( ltimecnt ) then
    ct0 = timecnt()
end if

if( myid_pw == 0 ) then
    ig1 = 2
  else
    ig1 = 1
end if

!-----set communicator
call get_worldpw( myid_, nodes_ )

do it = 1, ntype
if( .not.llking(it) ) then

   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)
!      do ig = 1, nplwcs
!         ftmp = svco(ig,it)*( rhgr(2*ig+1)*ysin(ig,ia)  &
!&                           + rhgr(2*ig+2)*ycos(ig,ia) )
!         floc(1,ia) = floc(1,ia) + ftmp*gx(ig)
!         floc(2,ia) = floc(2,ia) + ftmp*gy(ig)
!         floc(3,ia) = floc(3,ia) + ftmp*gz(ig)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
      bufst(1:3) = 0.d0
      do igg = ig1, nplw5nod
         ig = igg + nplw5nod1 - 2
         K1M = nga(ig)
         K2M = ngb(ig)
         K3M = ngc(ig)
         DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                DSQ1(K1M,ia)*DSQ2(K2M,ia)
         DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                DCQ1(K1M,ia)*DSQ2(K2M,ia)
         zcosk = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
         zsink = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)
         ftmp = svco(igg,it)*( rhgr(2*ig+1)*zsink  &
&                            + rhgr(2*ig+2)*zcosk )
!         floc(1,ia) = floc(1,ia) + ftmp*gx(ig)
!         floc(2,ia) = floc(2,ia) + ftmp*gy(ig)
!         floc(3,ia) = floc(3,ia) + ftmp*gz(ig)
         bufst(1) = bufst(1) + ftmp*gx(ig)
         bufst(2) = bufst(2) + ftmp*gy(ig)
         bufst(3) = bufst(3) + ftmp*gz(ig)
      end do
      if( nodes_pw > 1 ) then
          !-----global sum
          call gdsum(bufst,3,bufstr)
      end if
      floc(1,ia) = floc(1,ia) + bufst(1)
      floc(2,ia) = floc(2,ia) + bufst(2)
      floc(3,ia) = floc(3,ia) + bufst(3)
   end do

   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)
      floc(1,ia) = floc(1,ia)*2.d0
      floc(2,ia) = floc(2,ia)*2.d0
      floc(3,ia) = floc(3,ia)*2.d0
   end do

end if
end do

!-----set communicator
call get_worldkd( myid_, nodes_, nkd_ )


!--- stress calculation
if( lcstress ) then

    !-----set communicator
    call get_worldpw( myid_, nodes_ )

    bufst(1:6) = 0.d0

   do it = 1, ntype
   if( .not.llking(it) ) then

   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)

!      do ig = 1, nplwcs
!         stmp = svcop(ig,it)*( rhgr(2*ig+1)*ycos(ig,ia)  &
!&                            - rhgr(2*ig+2)*ysin(ig,ia) )
!         bufst(1) = bufst(1) + stmp*gx(ig)*gx(ig)
!         bufst(2) = bufst(2) + stmp*gy(ig)*gy(ig)
!         bufst(3) = bufst(3) + stmp*gz(ig)*gz(ig)
!         bufst(4) = bufst(4) + stmp*gy(ig)*gz(ig)
!         bufst(5) = bufst(5) + stmp*gz(ig)*gx(ig)
!         bufst(6) = bufst(6) + stmp*gx(ig)*gy(ig)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
      do igg = ig1, nplw5nod
         ig = igg + nplw5nod1 - 2
         K1M = nga(ig)
         K2M = ngb(ig)
         K3M = ngc(ig)
         DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                DSQ1(K1M,ia)*DSQ2(K2M,ia)
         DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                DCQ1(K1M,ia)*DSQ2(K2M,ia)
         zcosk = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
         zsink = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)
         stmp = svcop(igg,it)*( rhgr(2*ig+1)*zcosk  &
&                             - rhgr(2*ig+2)*zsink )
         bufst(1) = bufst(1) + stmp*gx(ig)*gx(ig)
         bufst(2) = bufst(2) + stmp*gy(ig)*gy(ig)
         bufst(3) = bufst(3) + stmp*gz(ig)*gz(ig)
         bufst(4) = bufst(4) + stmp*gy(ig)*gz(ig)
         bufst(5) = bufst(5) + stmp*gz(ig)*gx(ig)
         bufst(6) = bufst(6) + stmp*gx(ig)*gy(ig)
      end do

   end do

   end if
   end do
   if( nodes_pw > 1 ) then
       !-----global sum
       call gdsum(bufst,6,bufstr)
   end if

   !-----set communicator
   call get_worldkd( myid_, nodes_, nkd_ )

   call gdsum(bufst,6,bufstr)
   strloc(1,1) = strloc(1,1) + bufst(1)*2.d0 - elclv
   strloc(2,2) = strloc(2,2) + bufst(2)*2.d0 - elclv
   strloc(3,3) = strloc(3,3) + bufst(3)*2.d0 - elclv
   strloc(2,3) = strloc(2,3) + bufst(4)*2.d0
   strloc(3,1) = strloc(3,1) + bufst(5)*2.d0
   strloc(1,2) = strloc(1,2) + bufst(6)*2.d0
   strloc(3,2) = strloc(2,3)
   strloc(1,3) = strloc(3,1)
   strloc(2,1) = strloc(1,2)

end if

if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '    local pp in reciprocal space ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '    local pp in reciprocal space ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if


return
end




subroutine setnlc_g( nfile, myid, nodes,  &
& nplw, ntype, nhk1, nhk2, natom, lkbppi, lking,  &
& ycos, ysin, nplwcs, rvol )
!-----------------------------------------------------------------------
!     initialize reciprocal-space calculation
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw
integer :: ntype
integer :: natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
integer :: nplwcs
real(8) :: ycos(0:nplwcs,natom), ysin(0:nplwcs,natom)
real(8) :: rvol

!------declare local variables
integer :: it, l


do it = 1, ntype
   if( lkbppi(it) .and. .not.lking(it) ) then
       do l = 1, lmx(it)
          vclalp(l,it) = clalp(l,it)/rvol
       end do
   end if
end do

!---fully relativistic SO term 
!call fullrela_setnlc_g( nfile, myid, nodes,  &
!& ntype, lkbppi, lking, rvol )


!-----set eqlmr_
call setbgq( nfile, myid, nodes,  &
& nplw, ntype, nhk1, nhk2, natom, lkbppi, lking,  &
& ycos, ysin, nplwcs )


return
end subroutine




subroutine hckbpp_g( nfile, myid, nodes,  &
& cgjr, rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!     Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw
integer :: nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(nplwex,nbnod), rhcr(nplwex,nbnod)
integer :: ntype, natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
integer :: nplwcs
real*8  :: ycos(0:nplwcs,natom), ysin(0:nplwcs,natom)

!------declare local variables
integer :: nk, it, ikng, l, ia, ig, ib
real*8  :: tslmir, tslmii


call get_kvec( nk )

do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then

    if( .not.lbigmm2 ) then
        ikng = 1
        do l  = 1, lmx(it)
           call setqlm( nfile, myid, nodes,  &
& blmr(0,l), blmi(0,l), it, l, 0, nplw, nk, ikng )
        end do
    end if


do ia = nhk1(it), nhk2(it)

   if( .not.lbigmm2 ) then
       do l  = 1, lmx(it)
       if( litoll(l,it) ) then
           do ig = 0, nplw
              eqlmr(ig,l) = ycos(ig,ia)*blmr(ig,l)
              eqlmi(ig,l) = ysin(ig,ia)*blmr(ig,l)
           end do
         else
           do ig = 0, nplw
              eqlmr(ig,l) = - ysin(ig,ia)*blmi(ig,l)
              eqlmi(ig,l) =   ycos(ig,ia)*blmi(ig,l)
           end do
       end if
       end do
   end if


   do l = 1, lmx(it)
      if( lbigmm2 ) then
          do ig = 0, nplw
             zzzr(ig) = eqlmr_(2*ig+1,l,ia) * vclalp(l,it)
             zzzi(ig) = eqlmr_(2*ig+2,l,ia) * vclalp(l,it)
          end do
      else
          do ig = 0, nplw
             zzzr(ig) =  eqlmr(ig,l) * vclalp(l,it)
             zzzi(ig) = -eqlmi(ig,l) * vclalp(l,it)
          end do
      end if

      do ib = 1, nbnod
         tslmir = 0.d0
         tslmii = 0.d0

         if( lbigmm2 ) then
             do ig = 3, nplwex
                tslmir = tslmir + eqlmr_(ig,l,ia)*cgjr(ig,ib)
             end do
             ig = 1
             tslmir = 2.d0*tslmir  &
&                   + eqlmr_(ig,l,ia)*cgjr(ig,ib)
         else
             do ig = 1, nplw
                tslmir = tslmir + eqlmr(ig,l)*cgjr(2*ig+1,ib)  &
&                               - eqlmi(ig,l)*cgjr(2*ig+2,ib)
             end do
             ig = 0
             tslmir = 2.d0*tslmir  &
&                   + eqlmr(ig,l)*cgjr(2*ig+1,ib)
         end if
         slmir(l,ia,ib) = tslmir

         do ig = 0, nplw
            rhcr(2*ig+1,ib) = rhcr(2*ig+1,ib) + zzzr(ig)*tslmir
            rhcr(2*ig+2,ib) = rhcr(2*ig+2,ib) + zzzi(ig)*tslmir
         end do
      end do
   end do

end do
end if
end do


!      !-----unify slmir
!      call alldgatherv(slmir_nod,lmmax*natom*nbnod,slmir,nsmcnt,nsmdsp)


return
end




subroutine savesumr_g( nfile, myid, nodes,  &
& nbnod, ntype, nhk1, nhk2, natom, lkbppi, lking, lspin, nspin )
!-----------------------------------------------------------------------
!     store sumr
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
use nlkbpp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nbnod
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lkbppi, lking
logical :: lspin
integer :: nspin
!------declare local variables
integer :: it, nk


!if( lgamma .and. .not.lspin ) return

call get_kvec( nk )

do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then
    svsumr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&  = slmir(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
    if( .not.lgamma .or. lnoncollinear ) then
        svsumi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&      = slmii(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
    end if

    if( lnoncollinear ) then
        svsubr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&      = slmbr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
        svsubi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&      = slmbi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
    end if
end if
end do


return
end




subroutine loadsumr_g( nfile, myid, nodes,  &
& ntype, nhk1, nhk2, lkbppi, lking, lspin, nspin,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, iod, iodg, ioag, idstnd )
!-----------------------------------------------------------------------
!     load sumr
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
use nlkbpp_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ntype, nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
logical :: lspin
integer :: nspin
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(nodes), nbndsp(nodes)
integer :: natom, natnod1, natnod2, natnod, natcnt(nodes), natdsp(nodes)
integer :: iods(nbnod), iodsg(nband)
integer :: iod(nbnod), iodg(nband), ioag(natom), idstnd(nodes)

!------declare local variables
integer :: it, nk
real(8) :: ct0


!if( lgamma .and. .not.lspin ) return

call get_kvec( nk )

do it = 1, ntype
if( lkbppi(it) ) then
!     slmir(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)  &
!& = svsumr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)
     slmir(1:lmmax,nhk1(it):nhk2(it),1:nbnod)  &
& = svsumr(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)
end if
end do

!--- to convert band decomposition to atom decomposition
call slm_bdtovgd_org( nfile, myid, nodes, ct0,  &
& slmir, slmir_atm, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp, natnodx,  &
& iodg, ioag, idstnd, .true., .false. )

!--- to convert atom decomposition to band decomposition
call slm_vgdtobd_org( nfile, myid, nodes, ct0,  &
& slmir_atm, slmir, slmir_nod, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, ioag, idstnd, .false. )

do it = 1, ntype
if( lkbppi(it) ) then
!     svsumr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
!&   = slmir(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
     svsumr(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&   = slmir(1:lmmax,nhk1(it):nhk2(it),1:nbnod)
end if
end do


if( .not.lgamma .or. lnoncollinear ) then

    do it = 1, ntype
    if( lkbppi(it) ) then
!        slmii(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)  &
!&    = svsumi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)
        slmii(1:lmmax,nhk1(it):nhk2(it),1:nbnod)  &
&    = svsumi(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)
    end if
    end do

    !--- to convert band decomposition to atom decomposition
    call slm_bdtovgd_org( nfile, myid, nodes, ct0,  &
& slmii, slmir_atm, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp, natnodx,  &
& iodg, ioag, idstnd, .true., .false. )

    !--- to convert atom decomposition to band decomposition
    call slm_vgdtobd_org( nfile, myid, nodes, ct0,  &
& slmir_atm, slmii, slmir_nod, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, ioag, idstnd, .false. )

    do it = 1, ntype
    if( lkbppi(it) ) then
!     svsumi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
!&   = slmii(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
       svsumi(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&     = slmii(1:lmmax,nhk1(it):nhk2(it),1:nbnod)
    end if
    end do

end if


if( lnoncollinear ) then
    do it = 1, ntype
    if( lkbppi(it) ) then
!        slmbr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)  &
!&    = svsubr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)
        slmbr(1:lmmax,nhk1(it):nhk2(it),1:nbnod)  &
&    = svsubr(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)
    end if
    end do

    !--- to convert band decomposition to atom decomposition
    call slm_bdtovgd_org( nfile, myid, nodes, ct0,  &
& slmbr, slmir_atm, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp, natnodx,  &
& iodg, ioag, idstnd, .true., .false. )

    !--- to convert atom decomposition to band decomposition
    call slm_vgdtobd_org( nfile, myid, nodes, ct0,  &
& slmir_atm, slmbr, slmir_nod, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, ioag, idstnd, .false. )

    do it = 1, ntype
    if( lkbppi(it) ) then
!        svsubr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
!&      = slmbr(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
        svsubr(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&      = slmbr(1:lmmax,nhk1(it):nhk2(it),1:nbnod)
    end if
    end do


    do it = 1, ntype
    if( lkbppi(it) ) then
!        slmbi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)  &
!&    = svsubi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)
        slmbi(1:lmmax,nhk1(it):nhk2(it),1:nbnod)  &
&    = svsubi(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)
    end if
    end do

    !--- to convert band decomposition to atom decomposition
    call slm_bdtovgd_org( nfile, myid, nodes, ct0,  &
& slmbi, slmir_atm, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp, natnodx,  &
& iodg, ioag, idstnd, .true., .false. )

    !--- to convert atom decomposition to band decomposition
    call slm_vgdtobd_org( nfile, myid, nodes, ct0,  &
& slmir_atm, slmbi, slmir_nod, bufcr2,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lmmax, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, ioag, idstnd, .false. )

    do it = 1, ntype
    if( lkbppi(it) ) then
!     svsubi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
!&   = slmbi(1:lmx(it),nhk1(it):nhk2(it),1:nbnod)
       svsubi(1:lmmax,nhk1(it):nhk2(it),1:nbnod,nspin,nk)  &
&     = slmbi(1:lmmax,nhk1(it):nhk2(it),1:nbnod)
    end if
    end do
end if


return
end




subroutine hckbpp_ib_g( nfile, myid, nodes,  &
& cgjr, rhcr, nplwex, nplw,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!     Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw
real*8  :: cgjr(nplwex), rhcr(nplwex)
integer :: ntype, nhk1(ntype), nhk2(ntype), natom
logical :: lkbppi(ntype), lking(ntype)
integer :: nplwcs
real*8  :: ycos(0:nplwcs,natom), ysin(0:nplwcs,natom)

!------declare local variables
integer :: nk, it, ikng, l, ia, ig
real*8  :: tslmir, tslmii


call get_kvec( nk )

do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then

    if( .not.lbigmm2 ) then
        ikng = 1
        do l  = 1, lmx(it)
           call setqlm( nfile, myid, nodes,  &
& blmr(0,l), blmi(0,l), it, l, 0, nplw, nk, ikng )
        end do
    end if


do ia = nhk1(it), nhk2(it)

   if( .not.lbigmm2 ) then
       do l  = 1, lmx(it)
       if( litoll(l,it) ) then
           do ig = 0, nplw
              eqlmr(ig,l) = ycos(ig,ia)*blmr(ig,l)
              eqlmi(ig,l) = ysin(ig,ia)*blmr(ig,l)
           end do
         else
           do ig = 0, nplw
              eqlmr(ig,l) = - ysin(ig,ia)*blmi(ig,l)
              eqlmi(ig,l) =   ycos(ig,ia)*blmi(ig,l)
           end do
       end if
       end do
   end if


   do l = 1, lmx(it)
      if( lbigmm2 ) then
          do ig = 0, nplw
             zzzr(ig) = eqlmr_(2*ig+1,l,ia) * vclalp(l,it)
             zzzi(ig) = eqlmr_(2*ig+2,l,ia) * vclalp(l,it)
          end do
      else
          do ig = 0, nplw
             zzzr(ig) =  eqlmr(ig,l) * vclalp(l,it)
             zzzi(ig) = -eqlmi(ig,l) * vclalp(l,it)
          end do
      end if

         tslmir = 0.d0
         tslmii = 0.d0

         if( lbigmm2 ) then
             do ig = 3, nplwex
                tslmir = tslmir + eqlmr_(ig,l,ia)*cgjr(ig)
             end do
             ig = 1
             tslmir = 2.d0*tslmir + eqlmr_(ig,l,ia)*cgjr(ig)
         else
             do ig = 1, nplw
                tslmir = tslmir + eqlmr(ig,l)*cgjr(2*ig+1)  &
&                               - eqlmi(ig,l)*cgjr(2*ig+2)
             end do
             ig = 0
             tslmir = 2.d0*tslmir + eqlmr(ig,l)*cgjr(2*ig+1)
         end if

         do ig = 0, nplw
            rhcr(2*ig+1) = rhcr(2*ig+1) + zzzr(ig)*tslmir
            rhcr(2*ig+2) = rhcr(2*ig+2) + zzzi(ig)*tslmir
         end do
   end do

end do
end if
end do


return
end




subroutine enkbpp_g( nfile, myid, nodes,  &
& ecenl, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, rvol )
!-----------------------------------------------------------------------
!     energy from Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ecenl
integer :: nband, nbnod1, nbnod2, nbnod
real*8  :: occ(nband)
integer :: iod(nbnod)
integer :: ntype, natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
real*8  :: rvol


!if( .not.lnoncollinear ) then
    call enkbpp_g2( nfile, myid, nodes,  &
& ecenl, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, rvol )
!else
!    !-----noncollinear magnetism
!    call ncenkbpp_g2( nfile, myid, nodes,  &
!& ecenl, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
!& ntype, nhk1, nhk2, natom, lkbppi, lking, rvol )
!end if


return
end




subroutine enkbpp_g2( nfile, myid, nodes,  &
& ecenl, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, rvol )
!-----------------------------------------------------------------------
!     energy from Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ecenl
integer :: nband, nbnod1, nbnod2, nbnod
real*8  :: occ(nband)
integer :: iod(nbnod)
integer :: ntype, natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
real*8  :: rvol

!------declare local variables
logical :: lgammak
integer :: nk, i, ib, it, ia, l
real*8  :: ecenib, ecenia, dbuf1r


call get_kvec( nk )
call get_lgammak( lgammak )


do i = 1, nbnod
   ib = iod(i)
   if( abs(occ(ib)).gt.1.d-15 ) then

       ecenib = 0.0d0
       do it = 1, ntype
       if( lkbppi(it) .and. .not. lking(it) ) then
           do l = 1, lmx(it)
              ecenia = 0.d0
              if( lgamma .or. lgammak ) then
                  do ia = nhk1(it), nhk2(it)
                     ecenia = ecenia + slmir(l,ia,i)*slmir(l,ia,i)
                  end do
              else
                  do ia = nhk1(it), nhk2(it)
                     ecenia = ecenia + slmir(l,ia,i)*slmir(l,ia,i)  &
&                                    + slmii(l,ia,i)*slmii(l,ia,i)
                  end do
              end if
              ecenib = ecenib + ecenia * clalp(l,it)
           end do
       end if
       end do
       ecenl = ecenl + occ(ib)*ecenib/rvol

    end if
end do

call gdsum(ecenl,1,dbuf1r)


return
end




subroutine kbnlpp_frc( nfile, myid, nodes,  &
& fnlc, str, lcalstr, cgjr, nplwex, nplw, gx, gy, gz,  &
& nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!     Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
real*8  :: fnlc(3,natom)
real*8  :: str(6)
logical :: lcalstr
integer :: nplwex, nplw
real*8  :: gx(0:nplw),  gy(0:nplw),  gz(0:nplw)
integer :: nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(nplwex,nbnod)
real*8  :: occ(nband)
integer :: iod(nbnod)
integer :: ntype
integer :: nhk1(ntype), nhk2(ntype)
logical :: lkbppi(ntype), lking(ntype)
integer :: nplwcs
real*8  :: ycos(0:nplwcs,natom), ysin(0:nplwcs,natom)


!if( .not.lnoncollinear ) then
    call kbnlpp_frc2( nfile, myid, nodes,  &
& fnlc, str, lcalstr, cgjr, nplwex, nplw, gx, gy, gz,  &
& nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!else
!    !-----noncollinear magnetism
!    call nckbnlpp_frc_k22( nfile, myid, nodes,  &
!& fnlc, str, lcalstr, cgjr, nplwex, nplw, dkgx, dkgy, dkgz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iod,  &
!& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!end if


return
end




subroutine kbnlpp_frc2( nfile, myid, nodes,  &
& fnlc, str, lcalstr, cgjr, nplwex, nplw, gx, gy, gz,  &
& nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!     Kleinmann and Bylander type nonlocal pseudopotential
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
real*8,  dimension(3,natom) :: fnlc
real*8,  dimension(6) :: str
logical :: lcalstr
integer :: nplwex, nplw
real*8,  dimension(0:nplw) :: gx,  gy,  gz
integer :: nband, nbnod1, nbnod2, nbnod
real*8,  dimension(nplwex,nbnod) :: cgjr
real*8,  dimension(nband) :: occ
integer, dimension(nbnod) :: iod
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lkbppi, lking
real*8,  dimension(0:nplwcs,natom) :: ycos, ysin


call get_kvec( nk )

nplws1 = 0

do it = 1, ntype
if( lkbppi(it) .and. .not.lking(it) ) then
    ikng = 1
    if( .not.lbigmm2 ) then
        do l  = 1, lmx(it)
           call setqlm( nfile, myid, nodes,  &
& blmr(0,l), blmi(0,l), it, l, 0, nplw, nk, ikng )
        end do
    end if

    !--- stress calculation
    if( lcalstr ) then
        do l  = 1, lmx(it)
           call setqldm( nfile, myid, nodes,  &
& blmdr(0,1,l), blmdi(0,1,l), it, l, 0, nplw, nk, ikng, gx,gy,gz,  &
& lplnplw )
        end do
    end if


do ia = nhk1(it), nhk2(it)

   fnlc1 = 0.d0
   fnlc2 = 0.d0
   fnlc3 = 0.d0
   if( .not.lbigmm2 ) then
       do l  = 1, lmx(it)
       if( litoll(l,it) ) then
           do ig = 0, nplw
              eqlmr(ig,l) = ycos(ig,ia)*blmr(ig,l)
              eqlmi(ig,l) = ysin(ig,ia)*blmr(ig,l)
           end do
         else
           do ig = 0, nplw
              eqlmr(ig,l) = - ysin(ig,ia)*blmi(ig,l)
              eqlmi(ig,l) =   ycos(ig,ia)*blmi(ig,l)
           end do
       end if
       end do
   end if

   !-----h-f forces in non-local pseudopotential
#ifndef VECTOR
   do ibb = 1, nbnod
      ib = iod(ibb)
      if( abs(occ(ib)).gt.1.d-15 ) then
      do l = 1, lmx(it)
         fslxrl = 0.d0
         fslyrl = 0.d0
         fslzrl = 0.d0
         wegib2 = occ(ib)*2.d0 * vclalp(l,it)
         if( lbigmm2 ) then
             do ig = 1, nplw
                fzzr = eqlmr_(2*ig+1,l,ia)*cgjr(2*ig+2,ibb)  &
&                    - eqlmr_(2*ig+2,l,ia)*cgjr(2*ig+1,ibb)
                fslxrl = fslxrl + fzzr*gx(ig)
                fslyrl = fslyrl + fzzr*gy(ig)
                fslzrl = fslzrl + fzzr*gz(ig)
             end do
           else
             do ig = 1, nplw
                fzzr = eqlmr(ig,l)*cgjr(2*ig+2,ibb)  &
&                    + eqlmi(ig,l)*cgjr(2*ig+1,ibb)
                fslxrl = fslxrl + fzzr*gx(ig)
                fslyrl = fslyrl + fzzr*gy(ig)
                fslzrl = fslzrl + fzzr*gz(ig)
             end do
         end if
         fnltx = fslxrl*slmir(l,ia,ibb)
         fnlty = fslyrl*slmir(l,ia,ibb)
         fnltz = fslzrl*slmir(l,ia,ibb)
         fnlc1 = fnlc1 + fnltx*wegib2
         fnlc2 = fnlc2 + fnlty*wegib2
         fnlc3 = fnlc3 + fnltz*wegib2
!         fnlc(1,ia) = fnlc(1,ia) + fnltx*wegib2
!         fnlc(2,ia) = fnlc(2,ia) + fnlty*wegib2
!         fnlc(3,ia) = fnlc(3,ia) + fnltz*wegib2
      end do
      end if
   end do
#else
   do i = 1, lmx(it)*nbnod
      fslxr(i) = 0.d0
      fslyr(i) = 0.d0
      fslzr(i) = 0.d0
   end do
   ibl = 0
   do l =1,lmx(it)
      do ib = 1, nbnod
         ibl = ibl + 1
         if( lbigmm2 ) then
             do ig = 1, nplw
                fzzr = eqlmr_(2*ig+1,l,ia)*cgjr(2*ig+2,ib)  &
&                    - eqlmr_(2*ig+2,l,ia)*cgjr(2*ig+1,ib)
                fslxr(ibl) = fslxr(ibl) + fzzr*gx(ig)
                fslyr(ibl) = fslyr(ibl) + fzzr*gy(ig)
                fslzr(ibl) = fslzr(ibl) + fzzr*gz(ig)
             end do
           else
             do ig = 1, nplw
                fzzr = eqlmr(ig,l)*cgjr(2*ig+2,ib)  &
&                    + eqlmi(ig,l)*cgjr(2*ig+1,ib)
                fslxr(ibl) = fslxr(ibl) + fzzr*gx(ig)
                fslyr(ibl) = fslyr(ibl) + fzzr*gy(ig)
                fslzr(ibl) = fslzr(ibl) + fzzr*gz(ig)
             end do
         end if
      end do
   end do

   ibl = 0
   do l =1,lmx(it)
      do ibb = 1, nbnod
         ibl = ibl + 1
         ib = iod(ibb)
         dslmir(ibl) = slmir(l,ia,ibb)
         vclweg(ibl) = vclalp(l,it)*occ(ib)*2.d0
      end do
   end do

   do ibl = 1, lmx(it)*nbnod
      fnlc1 = fnlc1 + fslxr(ibl)*dslmir(ibl)*vclweg(ibl)
      fnlc2 = fnlc2 + fslyr(ibl)*dslmir(ibl)*vclweg(ibl)
      fnlc3 = fnlc3 + fslzr(ibl)*dslmir(ibl)*vclweg(ibl)
!      fnlc(1,ia) = fnlc(1,ia) + fslxr(ibl)*dslmir(ibl)*vclweg(ibl)
!      fnlc(2,ia) = fnlc(2,ia) + fslyr(ibl)*dslmir(ibl)*vclweg(ibl)
!      fnlc(3,ia) = fnlc(3,ia) + fslzr(ibl)*dslmir(ibl)*vclweg(ibl)
   end do
#endif
   fnlc(1,ia) = fnlc(1,ia) + 2.d0*fnlc1
   fnlc(2,ia) = fnlc(2,ia) + 2.d0*fnlc2
   fnlc(3,ia) = fnlc(3,ia) + 2.d0*fnlc3

   !--- stress calculation
   if( lcalstr ) then
       do l = 1, lmx(it)
       if( litoll(l,it) ) then
!OCL UNROLL(6)
           do iab = 1, 6
           do ig = 0, nplw
              eqlmdr(2*ig+1,iab,l) =   ycos(ig,ia)*blmdr(ig,iab,l)
              eqlmdr(2*ig+2,iab,l) = - ysin(ig,ia)*blmdr(ig,iab,l)
           end do
           end do
         else
!OCL UNROLL(6)
           do iab = 1, 6
           do ig = 0, nplw
              eqlmdr(2*ig+1,iab,l) = - ysin(ig,ia)*blmdi(ig,iab,l)
              eqlmdr(2*ig+2,iab,l) = - ycos(ig,ia)*blmdi(ig,iab,l)
           end do
           end do
       end if
       end do

       do ib = 1, nbnod
          do iab = 1, 6
          do l = 1, lmx(it)
             tslmir = 0.d0
             do ig = 3, nplwex
                tslmir = tslmir + eqlmdr(ig,iab,l)*cgjr(ig,ib)
             end do
             ig = 1
             tslmir = 2.d0*tslmir + eqlmdr(ig,iab,l)*cgjr(ig,ib)
             slmidr(l,ia,iab,ib) = tslmir
          end do
          end do
!                do iab = 1, 3
!                do l = 1, lmx(it)
!                   slmidr(l,ia,ib,iab) = slmidr(l,ia,ib,iab)
!     &                                 + 0.5d0*slmir(l,ia,ib)
!                   slmidi(l,ia,ib,iab) = slmidi(l,ia,ib,iab)
!     &                                 + 0.5d0*slmii(l,ia,ib)
!                end do
!                end do
       end do
   end if

end do
end if
end do


!do it = 1, ntype
!if( lkbppi(it) .and. .not.lking(it) ) then
!    do ia = nhk1(it), nhk2(it)
!       fnlc(1,ia) = 2.d0*fnlc(1,ia)
!       fnlc(2,ia) = 2.d0*fnlc(2,ia)
!       fnlc(3,ia) = 2.d0*fnlc(3,ia)
!    end do
!end if
!end do


!--- stress calculation
if( lcalstr ) then
    do ibb = 1, nbnod
       ib = iod(ibb)
       if( abs(occ(ib)).gt.1.d-15 ) then
       do iab = 1, 6
          ecenib = 0.0d0
          do it = 1,ntype
          if( lkbppi(it) .and. .not.lking(it) ) then
              do l = 1, lmx(it)
                 ecenia = 0.d0
                 do ia = nhk1(it), nhk2(it)
                    ecenia = ecenia  &
&                          + slmidr(l,ia,iab,ibb)*slmir(l,ia,ibb)
                 end do
                 ecenib = ecenib + 2.d0*ecenia * clalp(l,it)
              end do
          end if
          end do
          str(iab) = str(iab) + occ(ib)*ecenib
       end do
       end if
    end do
end if


return
end




subroutine setbgq( nfile, myid, nodes,  &
& nplw, ntype, nhk1, nhk2, natom, lnlcni, lking,  &
& ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!     initial set for nonlocal pp. for machines with BIG memory
!             set eqlmr_
!-----------------------------------------------------------------------
use nlpp_variables
use nl_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: nplw
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lnlcni, lking
real*8,  dimension(0:nplwcs,natom) :: ycos, ysin


if( .not.lbigmm2 ) return

do it = 1, ntype
if( lnlcni(it) .and. .not.lking(it) ) then

    nk   = 1
    ikng = 1
    do l  = 1, lmx(it)

       call setqlm( nfile, myid, nodes,  &
& zzzr, zzzi, it, l, 0, nplw, nk, ikng )

       do ia = nhk1(it), nhk2(it)
          if( litoll(l,it) ) then
              do ig = 0, nplw
                 eqlmr_(2*ig+1,l,ia) =   ycos(ig,ia)*zzzr(ig)
                 eqlmr_(2*ig+2,l,ia) = - ysin(ig,ia)*zzzr(ig)
              end do
            else
              do ig = 0, nplw
                 eqlmr_(2*ig+1,l,ia) = - ysin(ig,ia)*zzzi(ig)
                 eqlmr_(2*ig+2,l,ia) = - ycos(ig,ia)*zzzi(ig)
              end do
          end if
       end do

    end do

end if
end do


return
end




subroutine setqlm( nfile, myid, nodes,  &
& zzzr, zzzi, it, la, nplws1, noofpw, nk, ikng )
!-----------------------------------------------------------------------
!     calculate non-local pseudopotential term by table
!-----------------------------------------------------------------------
use nlpp_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: noofpw
real*8, dimension(0:noofpw) :: zzzr, zzzi
integer :: it, la, nk, ikng


l   = itola(la,it)
lm  = itolm(la,it)
lmf = itolmf(la,it)
do ig = nplws1, noofpw
       m = abgvec(ig,nk)/(2.d0*dlqlm(it))
       m = 2*m
       d = 0.5d0*( abgvec(ig,nk)/dlqlm(it) - dble(m) )
       zzzr(ig) = d*( (d-1.d0)*tbqlma(m,l,it) + tbqlm(m+2,l,it)  &
&                           - tbqlm(m,l,it) ) + tbqlm(m,l,it)
end do
if( lmf.eq.1 ) then
        if( litoll(la,it) ) then
            do ig = nplws1, noofpw
               buffr = zzzr(ig)
               zzzr(ig) = buffr*qlmr(ig,lm,nk,ikng)
               zzzi(ig) = 0.d0
            end do
          else
            do ig = nplws1, noofpw
               buffr = zzzr(ig)
               zzzr(ig) = 0.d0
               zzzi(ig) = buffr*qlmr(ig,lm,nk,ikng)
            end do
        end if
  else
        if( litoll(la,it) ) then
            do ig = nplws1, noofpw
               buffr = zzzr(ig)
               zzzr(ig) = buffr*qlmi(ig,lm,nk,ikng)
               zzzi(ig) = 0.d0
            end do
          else
            do ig = nplws1, noofpw
               buffr = zzzr(ig)
               zzzr(ig) = 0.d0
               zzzi(ig) = buffr*qlmi(ig,lm,nk,ikng)
            end do
        end if
end if


return
end




subroutine setqldm( nfile, myid, nodes,  &
& zzzr, zzzi, it, la, nplws1, noofpw, nk, ikng, gx,  gy,  gz,  &
& ndimpw )
!-----------------------------------------------------------------------
!     --- stress calculation ---
!     calculate non-local pseudopotential term by table
!-----------------------------------------------------------------------
use nlpp_variables
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: nplws1, noofpw, ndimpw
real*8, dimension(0:ndimpw,6) :: zzzr, zzzi
real*8,  dimension(0:noofpw) :: gx,  gy,  gz
integer :: it, la, nk, ikng


l   = itola(la,it)
lm  = itolm(la,it)
lmf = itolmf(la,it)
do ig = nplws1, noofpw
       m = abgvec(ig,nk)/(2.d0*dlqlm(it))
       m = 2*m
       d = 0.5d0*( abgvec(ig,nk)/dlqlm(it) - dble(m) )
       zzzr(ig,1) = d*( (d-1.d0)*tbqlma(m,l,it) + tbqlm(m+2,l,it)  &
&                           - tbqlm(m,l,it) ) + tbqlm(m,l,it)
!             m = abgvec(ig,nk)/(2.d0*dlqld(it))
!             m = 2*m
!             d = 0.5d0*( abgvec(ig,nk)/dlqld(it) - dble(m) )
       zzzr(ig,2) = d*( (d-1.d0)*tbqlda(m,l,it) + tbqld(m+2,l,it)  &
&                           - tbqld(m,l,it) ) + tbqld(m,l,it)
end do
do ig = nplws1, noofpw
!             if( ig.ne.0 ) then
       if( abgvec(ig,nk) > 1.d-12 ) then
           zzzr(ig,3) = 1.d0/abgvec(ig,nk)
         else
           zzzr(ig,3) = 0.d0
       end if
end do
if( lmf.eq.1 ) then
    if( litoll(la,it) ) then
            do ig = nplws1, noofpw
               buffr  = zzzr(ig,1)
               buffrh = buffr*0.5d0
               buffs  = zzzr(ig,2)
               stmp   = zzzr(ig,3)
               zzzr(ig,1) = buffr*qlmdgr(ig,lm,1,nk)*gx(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gx(ig)*gx(ig)
               zzzr(ig,2) = buffr*qlmdgr(ig,lm,2,nk)*gy(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gy(ig)*gy(ig)
               zzzr(ig,3) = buffr*qlmdgr(ig,lm,3,nk)*gz(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gz(ig)*gz(ig)
               zzzr(ig,4) = buffrh*( qlmdgr(ig,lm,2,nk)*gz(ig)  &
&                                  + qlmdgr(ig,lm,3,nk)*gy(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gy(ig)*gz(ig)
               zzzr(ig,5) = buffrh*( qlmdgr(ig,lm,3,nk)*gx(ig)  &
&                                  + qlmdgr(ig,lm,1,nk)*gz(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gz(ig)*gx(ig)
               zzzr(ig,6) = buffrh*( qlmdgr(ig,lm,1,nk)*gy(ig)  &
&                                  + qlmdgr(ig,lm,2,nk)*gx(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gx(ig)*gy(ig)
            end do
            do iab = 1, 6
            do ig = nplws1, noofpw
               zzzi(ig,iab) = 0.d0
            end do
            end do
      else
            do ig = nplws1, noofpw
               buffr  = zzzr(ig,1)
               buffrh = buffr*0.5d0
               buffs  = zzzr(ig,2)
               stmp   = zzzr(ig,3)
               zzzi(ig,1) = buffr*qlmdgr(ig,lm,1,nk)*gx(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gx(ig)*gx(ig)
               zzzi(ig,2) = buffr*qlmdgr(ig,lm,2,nk)*gy(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gy(ig)*gy(ig)
               zzzi(ig,3) = buffr*qlmdgr(ig,lm,3,nk)*gz(ig)  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gz(ig)*gz(ig)
               zzzi(ig,4) = buffrh*( qlmdgr(ig,lm,2,nk)*gz(ig)  &
&                                  + qlmdgr(ig,lm,3,nk)*gy(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gy(ig)*gz(ig)
               zzzi(ig,5) = buffrh*( qlmdgr(ig,lm,3,nk)*gx(ig)  &
&                                  + qlmdgr(ig,lm,1,nk)*gz(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gz(ig)*gx(ig)
               zzzi(ig,6) = buffrh*( qlmdgr(ig,lm,1,nk)*gy(ig)  &
&                                  + qlmdgr(ig,lm,2,nk)*gx(ig) )  &
&                   + buffs*qlmr(ig,lm,nk,ikng)*stmp*gx(ig)*gy(ig)
            end do
            do iab = 1, 6
            do ig = nplws1, noofpw
               zzzr(ig,iab) = 0.d0
            end do
            end do
    end if
  else
    if( litoll(la,it) ) then
            do ig = nplws1, noofpw
               buffr  = zzzr(ig,1)
               buffrh = buffr*0.5d0
               buffs  = zzzr(ig,2)
               stmp   = zzzr(ig,3)
               zzzr(ig,1) = buffr*qlmdgi(ig,lm,1,nk)*gx(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gx(ig)*gx(ig)
               zzzr(ig,2) = buffr*qlmdgi(ig,lm,2,nk)*gy(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gy(ig)*gy(ig)
               zzzr(ig,3) = buffr*qlmdgi(ig,lm,3,nk)*gz(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gz(ig)*gz(ig)
               zzzr(ig,4) = buffrh*( qlmdgi(ig,lm,2,nk)*gz(ig)  &
&                                  + qlmdgi(ig,lm,3,nk)*gy(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gy(ig)*gz(ig)
               zzzr(ig,5) = buffrh*( qlmdgi(ig,lm,3,nk)*gx(ig)  &
&                                  + qlmdgi(ig,lm,1,nk)*gz(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gz(ig)*gx(ig)
               zzzr(ig,6) = buffrh*( qlmdgi(ig,lm,1,nk)*gy(ig)  &
&                                  + qlmdgi(ig,lm,2,nk)*gx(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gx(ig)*gy(ig)
            end do
            do iab = 1, 6
            do ig = nplws1, noofpw
               zzzi(ig,iab) = 0.d0
            end do
            end do
      else
            do ig = nplws1, noofpw
               buffr  = zzzr(ig,1)
               buffrh = buffr*0.5d0
               buffs  = zzzr(ig,2)
               stmp   = zzzr(ig,3)
               zzzi(ig,1) = buffr*qlmdgi(ig,lm,1,nk)*gx(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gx(ig)*gx(ig)
               zzzi(ig,2) = buffr*qlmdgi(ig,lm,2,nk)*gy(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gy(ig)*gy(ig)
               zzzi(ig,3) = buffr*qlmdgi(ig,lm,3,nk)*gz(ig)  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gz(ig)*gz(ig)
               zzzi(ig,4) = buffrh*( qlmdgi(ig,lm,2,nk)*gz(ig)  &
&                                  + qlmdgi(ig,lm,3,nk)*gy(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gy(ig)*gz(ig)
               zzzi(ig,5) = buffrh*( qlmdgi(ig,lm,3,nk)*gx(ig)  &
&                                  + qlmdgi(ig,lm,1,nk)*gz(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gz(ig)*gx(ig)
               zzzi(ig,6) = buffrh*( qlmdgi(ig,lm,1,nk)*gy(ig)  &
&                                  + qlmdgi(ig,lm,2,nk)*gx(ig) )  &
&                   + buffs*qlmi(ig,lm,nk,ikng)*stmp*gx(ig)*gy(ig)
            end do
            do iab = 1, 6
            do ig = nplws1, noofpw
               zzzr(ig,iab) = 0.d0
            end do
            end do
    end if
end if


return
end




