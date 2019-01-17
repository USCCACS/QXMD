



module fermi_variables
!-----------------------------------------------------------------------
! type declaration of variables in fermi.f
!-----------------------------------------------------------------------
implicit none

!-----variables for Gaussian broadening
integer, parameter :: ntbfi = 10000
real*8,  allocatable, dimension(:) :: tbfi ,tbfia
real*8 :: tbxe
real*8 :: tbmx
real*8 :: piru

!-----work arrays
real*8,  allocatable, dimension(:) :: awork
real*8,  allocatable, dimension(:) :: bwork

!-----for Methfessel & Paxton smearing, PRB 40, 3616 (1989) (lfermi > 3)
real*8,  allocatable, dimension(:) :: coefdelta  ! = An
real*8,  allocatable, dimension(:) :: hpoly      ! = Hn
integer :: npolyorder

!--- for sorting ------------------------------
integer, allocatable, dimension(:) :: key, ix, dstkey, x1, x2
real*8,  parameter :: dkeymx = 1000000000.d0

real*8  :: small = 1.d-13

logical :: lnoncollinear

logical :: lcal_sop = .false.
save

end module




subroutine set_lnoncollinear_in_fermi( nfile, myid, nodes, lnoncollinear_ )
!-----------------------------------------------------------------------
! set lnoncollinear
!-----------------------------------------------------------------------
use fermi_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lnoncollinear_

lnoncollinear = lnoncollinear_

return
end subroutine




subroutine set_lcal_sop_in_fermi( nfile, myid, nodes, lcal_sop_ )
!-----------------------------------------------------------------------
! set lcal_sop
!-----------------------------------------------------------------------
use fermi_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lcal_sop_

lcal_sop = lcal_sop_

return
end subroutine




subroutine fermi_alloc( nfile, myid, nodes, alloc_mem, nband )
!-----------------------------------------------------------------------
!     allocate memory for work arrays in fermi.f
!-----------------------------------------------------------------------
use fermi_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: nband

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( tbfi(0:ntbfi) ,tbfia(0:ntbfi),  &
& awork(nband), bwork(nband),  &
& key(nband), ix(nband), dstkey(nband),  &
& x1(nband), x2(nband), stat=status )

the_mem = 8.d0 * nband * 2 + 4.d0 * nband * 5

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'fermi_alloc', .true. )


return
end subroutine




subroutine iniocc( nfile, myid, nodes,  &
& nel, nocc, nband, occ, norder, nkpnt, wbzk )
!-----------------------------------------------------------------------
!     set initial occupancies
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: nel
integer :: nocc
integer :: nband
integer :: nkpnt
real*8,  dimension(nband,nkpnt) :: occ
integer, dimension(nband,nkpnt) :: norder
real*8,  dimension(nkpnt) :: wbzk

!-----declare local variables
integer :: ib, k


do k = 1, nkpnt

nocc  = nel/2
do ib = 1, nocc
   occ(ib,k) = 2.d0*wbzk(k)
end do
if( mod(nel,2) /= 0 ) then
    nocc = nocc + 1
    occ(nocc,k) = 1.d0*wbzk(k)
end if
do ib = nocc + 1, nband
   occ(ib,k) = 0.d0
end do

do ib = 1, nband
   norder(ib,k) = ib
end do

end do
!-----------------------------------------------------------------------
!      nocc = 4
!      occ(1) = 2.d0
!      occ(2) = 2.d0/3.d0
!      occ(3) = 2.d0/3.d0
!      occ(4) = 2.d0/3.d0


return
end subroutine




subroutine iniocc_spin( nfile, myid, nodes,  &
& nel, nocc, nband, occ, norder, wegud, nsordr, nkpnt, wbzk )
!-----------------------------------------------------------------------
!     set initial occupancies
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: nel
integer :: nocc
integer :: nband
integer :: nkpnt
real*8,  dimension(nband,nkpnt) :: occ
integer, dimension(nband,nkpnt) :: norder
real*8,  dimension(nband,2,nkpnt) :: wegud
integer, dimension(nband,2,nkpnt) :: nsordr
real*8,  dimension(nkpnt) :: wbzk

!-----declare local variables
integer :: ib, k


do k = 1, nkpnt

nocc  = nel/2
do ib = 1, nocc
   occ(ib,k) = 1.d0*wbzk(k)
   wegud(ib,1,k) = 1.d0*wbzk(k)
   wegud(ib,2,k) = 1.d0*wbzk(k)
end do
if( mod(nel,2) /= 0 ) then
    nocc = nocc + 1
    occ(nocc,k) = 0.5d0*wbzk(k)
   wegud(nocc,1,k) = 1.d0*wbzk(k)
   wegud(nocc,2,k) = 0.d0
end if
do ib = nocc + 1, nband
   occ(ib,k) = 0.d0
   wegud(ib,1,k) = 0.d0
   wegud(ib,2,k) = 0.d0
end do

do ib = 1, nband
   norder(ib,k) = ib
   nsordr(ib,1,k) = ib
   nsordr(ib,2,k) = ib
end do

end do
!-----------------------------------------------------------------------
!      nocc = 4
!      occ(1) = 2.d0
!      occ(2) = 2.d0/3.d0
!      occ(3) = 2.d0/3.d0
!      occ(4) = 2.d0/3.d0


return
end subroutine




subroutine fermi_k( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, eig, norder, occ, wbzk, nocc,  &
& entrpy, fermie, lfermi, tfermi )
!-----------------------------------------------------------------------
!     occupancies
!
! (input)
!     nband  ...... No. of electronic bands
!     nel    ...... No. of electrons
!     eig    ...... eigenvalues
!     lfermi ...... = 1  non-metallic
!                     2  metallic with Fermi distribution
!                     3  metallic with Gaussian distribution
!            lfermi > 3  Methfessel & Paxton smearing, PRB 40, 3616 (1989)
!                        The order of approximant = lfermi - 3
!     tfermi ...... Fermi temperature, if lfermi = 2 or 3.
!
! (output)
!     norder ...... order of each state
!     occ    ...... occupancies
!     nocc   ...... No. of occupied states
!     entrpy ...... entropy by fractional occupations
!-----------------------------------------------------------------------
use outfile
use fermi_variables
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8,  dimension(nband*nkpnt) :: eig, occ
integer, dimension(nband,nkpnt) :: norder
real*8,  dimension(nkpnt) :: wbzk

!----- for Gaussian broadening -----
INTEGER :: status
logical :: licall = .true.
save       licall


if( .not.lnoncollinear .and. .not.lcal_sop ) then
    two  = 2.d0
    half = 0.5d0
else
    two  = 1.d0
    half = 1.d0
end if

nbandk = nband*nkpnt
!--- sorting by eigenvalues : eig --------------------------------------
!     ix     : data pointer of record

ddmx = 0.d0
do ib = 1, nbandk
   ddmx = max( ddmx, abs(eig(ib)) )
end do

fctmax = dkeymx/ddmx
do ib = 1, nbandk
   key(ib) = fctmax*eig(ib)
   ix(ib)  = ib
end do

if( nbandk.gt.1 ) then
    call vsgar( nbandk, key, ix, nbandk, dstkey, x1, x2, x2 )
    call vsoe(  nbandk, key, ix, nbandk )
end if
!-----------------------------------------------------------------------

do k = 1, nkpnt
   key(k) = 0
end do
do ibk = 1, nbandk
   k  = ( ix(ibk) - 1 )/nband + 1
   ib = ix(ibk) - (k-1)*nband
   key(k) = key(k) + 1
   norder(key(k),k) = ib
end do
ii = 0
do k = 1, nkpnt
do ib = 1, nband
   ii = ii + 1
   bwork(ii) = wbzk(k)
end do
end do


entrpy = 0.d0

!nocc = nel/2
!if( mod(nel,2).ne.0 ) nocc = nocc + 1
nocc = nel/nint(two)
if( mod(nel,nint(two)).ne.0 ) nocc = nocc + 1
if( nband.lt.nocc ) then
    call fstop( nfile, myid, nodes, 'nband.lt.nocc in fermi_k' )
end if

!--- non metallic ------------------------------------------------------
    occt = 0.d0
    noccbk = nbandk
    do ibk = 1, nbandk
       ib = ix(ibk)
!       occ(ib) = 2.d0*bwork(ib)
       occ(ib) = two*bwork(ib)
       occt = occt + occ(ib)
       if( occt+small >= dble(nel) ) then
           noccbk = ibk
           occ(ib) = occ(ib) - ( occt - dble(nel) )
           exit
       end if
    end do
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       occ(ib) = 0.d0
    end do

if( lfermi.eq.1 .or. nband.eq.nocc ) then
!--- check degeneracy
    toldeg = 1.d-02
    refeig = eig(ix(noccbk))
    ib1 = 1
    do ibk = noccbk - 1, 1, -1
       ib = ix(ibk)
       if( abs(refeig-eig(ib)) > toldeg ) then
           ib1 = ibk + 1
           exit
       end if
    end do
    ib2 = nbandk
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       if( abs(refeig-eig(ib)) > toldeg ) then
           ib2 = ibk - 1
           exit
       end if
    end do
    del = 0.d0
    wei = 0.d0
    do ibk = ib1, ib2
       ib = ix(ibk)
       del = del + occ(ib)
       wei = wei + bwork(ib)
    end do
    del = del/wei
    do ibk = ib1, ib2
       ib = ix(ibk)
       occ(ib) = del*bwork(ib)
    end do
!-----------------------------------------------------------------------
    return
end if

!--- metallic ----------------------------------------------------------
eps   = 1.0d-08
if( lfermi.eq.2 ) then
    eps2   = 1.0d-12
  else
    eps2   = 1.0d-12
    if( licall ) then
        !----- initial set for Gaussian broadening -----
        call setgbr( tbfi ,tbfia, ntbfi, tbxe, tbmx, piru )

        !---the oder of approximant
        npolyorder = lfermi - 3

        !---alocate memory
        allocate( coefdelta(0:npolyorder), hpoly(0:2*npolyorder),  &
& stat=status )

        !---set coefficients An
        coefdelta(0) = 1.d0/sqrt(acos(-1.d0))
        do i = 1, npolyorder
           coefdelta(i) = coefdelta(i-1)*(-1.d0)/4.d0/dble(i)
        end do
    end if
    licall = .false.
end if


enel = nel
zero = 0.d0
one  = 1.d0
!half = 0.5d0
!--- estimate Fermi energy by False position method
nocclw = noccbk
nocchg = noccbk + 1
ib = ix(nocclw)
e1 = eig(ib)
ib = ix(nocchg)
e2 = eig(ib)
!---      fermie = 0.5d0*(e1+e2)
itt = 0
do
   itt = itt + 1
   if( itt == 1 ) then
      fermie = e1
   else if( itt == 2 ) then
      fermie = e2
   else if( itt == 3 ) then
      fermie = (e1 + e2)*0.5d0
   else
      fermie = e1 - (e2-e1)/(d2-d1)*d1
   end if

      sum = zero
      if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
          do ib = 1, nbandk
              a = (eig(ib)-fermie)/tfermi
              if(a.gt.3.5d+01) then
                  occ(ib) = zero
              else if(a.lt.-3.3d+01) then
                  occ(ib) = one
              else
                  occ(ib) = one/(one+dexp(a))
              end if
              occ(ib) = occ(ib)*bwork(ib)
              sum = sum + occ(ib)
          end do
        else
!     --- Gaussian broadening ---
          do ib = 1, nbandk
              a = (eig(ib)-fermie)/tfermi
              if(a.gt.tbmx) then
                  occ(ib) = zero
              else if(a.lt.-tbmx) then
                  occ(ib) = one
              else
                  dm = (a+tbmx)/tbxe
                  m  = 0.5d0*dm
                  m  = 2*m
                  d  = 0.5d0*( dm - dble(m) )
                  occ(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
     &                                            - tbfi(m) ) + tbfi(m)
              end if
              occ(ib) = occ(ib)*bwork(ib)
              sum = sum + occ(ib)
          end do

          !---Methfessel & Paxton smearing
          if( npolyorder > 0 ) then
              do ib = 1, nbandk
                 a = (eig(ib)-fermie)/tfermi
                 if( abs(a) > tbmx ) cycle
                 hpoly(0) = 1.d0
                 hpoly(1) = 2.d0*a
                 do j = 2, 2*npolyorder - 1
                    hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
                 end do
                 gfi = 0.d0
                 do j = 1, npolyorder
                    gfi = gfi + coefdelta(j)*hpoly(2*j-1)
                 end do
                 gfi = gfi * exp(-a*a) * bwork(ib)
                 occ(ib) = occ(ib) + gfi
                 sum     = sum     + gfi
              end do
          end if
      end if

      d = sum - enel*half
!      if( myid == 0 ) then
!          write(*,*) itt, d, sum
!          write(*,*) fermie, enel*half
!      end if
      if( dabs(d) < 1.d-03 ) exit

   if( itt == 1 ) then
       d1 = d
       if( d > 0.d0) then
           if( nocclw == 1 ) then
               ib = ix(nocclw)
               e1 = e1 - abs(eig(ib))*0.5d0
           else
               nocclw = nocclw - 1
               ib = ix(nocclw)
               e1 = eig(ib)
           end if
           itt = 0
       end if
   else if( itt == 2 ) then
       d2 = d
       if( d < 0.d0) then
           if( nocchg == nbandk ) then
               ib = ix(nocchg)
               e2 = e2 + abs(eig(ib))*0.5d0
           else
               nocchg = nocchg + 1
               ib = ix(nocchg)
               e2 = eig(ib)
           end if
           itt = 1
       end if
   else
       if( d1*d < 0.d0 ) then
           e2 = fermie
           d2 = d
         else
           e1 = fermie
           d1 = d
       end if
   end if
end do


if(.not.(dabs(d).lt.eps)) then
!--- calculate fermi energy by newton iteration
    itr =0
    sumv=zero
    nchk=0
    do ib = 1, nbandk
       a = (eig(ib)-fermie)/tfermi
       if( a.lt.-3.3d+01 ) then
           sumv = sumv + one*bwork(ib)
       else if( a.le.3.5d+01 ) then
            nchk=nchk+1
            awork(nchk)=eig(ib)
            bwork(nchk)=bwork(ib)
       end if
    end do
    c = enel*half - sumv

    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                px =px+one/(one+y)   * bwork(ib)
                pdx=pdx+y/(one+y)**2 * bwork(ib)
            end if
        end do

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

      else
!         --- Gaussian broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.tbmx)) then
                if(a.lt.-tbmx) then
                    px  = px + one * bwork(ib)
                  else
                    dm = (a+tbmx)/tbxe
                    m  = 0.5d0*dm
                    m  = 2*m
                    d  = 0.5d0*( dm - dble(m) )
                    gfi = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
                    px  = px + gfi            * bwork(ib)
                    pdx = pdx+ exp(-a*a)*piru * bwork(ib)
                end if
            end if
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nchk
               a = ( awork(ib) - fermie )/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
                gfi = 0.d0
               dgfi = 0.d0
               do j = 1, npolyorder
                   gfi =  gfi + coefdelta(j)*hpoly(2*j-1)
                  dgfi = dgfi + coefdelta(j)*hpoly(2*j)
               end do
                gfi =  gfi * exp(-a*a) * bwork(ib)
               dgfi = dgfi * exp(-a*a) * bwork(ib)
               px  = px  +  gfi
               pdx = pdx + dgfi
            end do
        end if

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

    end if

!--- check the sum of weighting function
ii = 0
do k = 1, nkpnt
do ib = 1, nband
   ii = ii + 1
   bwork(ii) = wbzk(k)
end do
end do

    sum=zero
    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---
        do ib = 1, nbandk
           a = (eig(ib)-fermie)/tfermi
           if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                occ(ib)= one/(one+y) * bwork(ib)
                sum    = sum + occ(ib)
           else
                occ(ib) = 0.d0
           end if
        end do
      else
!         --- Gaussian broadening ---
        do ib = 1, nbandk
            a = (eig(ib)-fermie)/tfermi
            if(a.gt.tbmx) then
                occ(ib) = zero
            else if(a.lt.-tbmx) then
                occ(ib) = one
            else
                dm = (a+tbmx)/tbxe
                m  = 0.5d0*dm
                m  = 2*m
                d  = 0.5d0*( dm - dble(m) )
                occ(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
            end if
            occ(ib) = occ(ib) * bwork(ib)
            sum = sum + occ(ib)
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nbandk
               a = (eig(ib)-fermie)/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder - 1
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
               gfi = 0.d0
               do j = 1, npolyorder
                  gfi = gfi + coefdelta(j)*hpoly(2*j-1)
               end do
               gfi = gfi * exp(-a*a) * bwork(ib)
               occ(ib) = occ(ib) + gfi
               sum     = sum     + gfi
            end do
        end if
    end if

    d=sum-enel*half
end if

!--- small correction
weizzz = 0.d0
do ib = 1, nbandk
   weizzz = weizzz + ( bwork(ib) - occ(ib) )*( occ(ib) - zero )
end do
if( abs(weizzz).gt.1.d-30 ) then
    do ib = 1, nbandk
       wzzz = ( bwork(ib) - occ(ib) )*( occ(ib) - zero )/weizzz
       occ(ib) = occ(ib) - wzzz*d
    end do
end if

if( dabs(d).gt.eps ) then
    if(loutfile(1)) write(nfile(1),1020) sum
    if(loutfile(2)) write(nfile(2),1020) sum
end if
1020     format('*** sub. fermi, error in sum of weighting'/  &
&                  '*** sum=',es20.12)

!--- calculate entropy term
entrpy=zero
if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
    do ib = 1, nbandk
       occib = occ(ib)/bwork(ib)
       if(.not.( occib.lt.eps .or. occib.gt.one-eps)) then
            entrpy = entrpy  &
&                  + bwork(ib)*( occib*dlog(occib)  &
&                              + (one-occib)*dlog(one-occib) )
       end if
    end do
!    entrpy = -2.0d0*tfermi*entrpy
    entrpy = -two*tfermi*entrpy
else
!     --- Gaussian broadening ---
    do ib = 1, nbandk
       a = (eig(ib)-fermie)/tfermi
       a = a*a
       if( a.lt.3.3d+01 ) then
           entrpy = entrpy + exp(-a) * bwork(ib)
        end if
    end do
!    entrpy = tfermi*entrpy*piru
    entrpy = tfermi*entrpy*piru*two*0.5d0

    !---Methfessel & Paxton smearing
    if( npolyorder > 0 ) then
        entrpymp = 0.d0
        do ib = 1, nbandk
           a = (eig(ib)-fermie)/tfermi
           if( abs(a) > tbmx ) cycle
           hpoly(0) = 1.d0
           hpoly(1) = 2.d0*a
           do j = 2, 2*npolyorder
              hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
           end do
           j = npolyorder
           entrpymp = entrpymp + hpoly(2*j) * exp(-a*a) * bwork(ib)
        end do
        j = npolyorder
!        entrpy = tfermi*coefdelta(j)*entrpymp
        entrpy = tfermi*coefdelta(j)*entrpymp*two*0.5d0
    end if
end if


!--- spin up & down double occupancy
do ib = 1, nbandk
!   occ(ib) = occ(ib)*2.d0
   occ(ib) = occ(ib)*two
end do

nocc = nband


return
end




subroutine fermi_k_spin( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, egvlud, nsordr, wegud, wbzk, nocc,  &
& eig, norder, occ,  &
& entrpy, feneud, fermie, lfermi, tfermi, difupd, lfixud )
!-----------------------------------------------------------------------
!     occupancies
!
! (input)
!     nband  ...... No. of electronic bands
!     nel    ...... No. of electrons
!     eig    ...... eigenvalues
!     lfermi ...... = 1  non-metallic
!                     2  metallic with Fermi distribution
!                     3  metallic with Gaussian distribution
!            lfermi > 3  Methfessel & Paxton smearing, PRB 40, 3616 (1989)
!                        The order of approximant = lfermi - 3
!     tfermi ...... Fermi temperature, if lfermi = 2 or 3.
!
! (output)
!     norder ...... order of each state
!     occ    ...... occupancies
!     nocc   ...... No. of occupied states
!     entrpy ...... entropy by fractional occupations
!-----------------------------------------------------------------------
use outfile
use fermi_variables
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8,  dimension(nband,2,nkpnt) :: egvlud
real*8,  dimension(nband*2*nkpnt) :: wegud
integer, dimension(nband,2,nkpnt) :: nsordr
real*8,  dimension(nband*nkpnt) :: eig, occ
integer, dimension(nband,nkpnt) :: norder
real*8,  dimension(nkpnt) :: wbzk
real*8,  dimension(2) :: deleud, feneud
logical :: lfixud

!----- for Gaussian broadening -----
integer :: status
logical :: licall = .true.
save       licall


nocc = nel/2
if( mod(nel,2).ne.0 ) nocc = nocc + 1
if( nband.lt.nocc ) then
    call fstop( nfile, myid, nodes, 'nband.lt.nocc in fermi_k_s' )
end if


!     No. of electrons
enel = nel
deleud(1) = 0.5d0*( enel + difupd )
deleud(2) = 0.5d0*( enel - difupd )

nbandk = nband*nkpnt

entrpy = 0.d0

spindo: do nspin = 1, 2
   ii = 0
   do k = 1, nkpnt
   do ib = 1, nband
      ii = ii + 1
      eig(ii) = egvlud(ib,nspin,k)
   end do
   end do
!--- sorting by eigenvalues : eig --------------------------------------
!     ix     : data pointer of record

ddmx = 0.d0
do ib = 1, nbandk
   ddmx = max( ddmx, abs(eig(ib)) )
end do

fctmax = dkeymx/ddmx
do ib = 1, nbandk
   key(ib) = fctmax*eig(ib)
   ix(ib)  = ib
end do

if( nbandk.gt.1 ) then
    call vsgar( nbandk, key, ix, nbandk, dstkey, x1, x2, x2 )
    call vsoe(  nbandk, key, ix, nbandk )
end if
!-----------------------------------------------------------------------

do k = 1, nkpnt
   key(k) = 0
end do
do ibk = 1, nbandk
   k  = ( ix(ibk) - 1 )/nband + 1
   ib = ix(ibk) - (k-1)*nband
   key(k) = key(k) + 1
   nsordr(key(k),nspin,k) = ib
end do


fixif: if( lfixud ) then

   ii = 0
   do k = 1, nkpnt
   do ib = 1, nband
      ii = ii + 1
      bwork(ii) = wbzk(k)
   end do
   end do
!--- non metallic ------------------------------------------------------
    occt = 0.d0
    noccbk = nbandk
    do ibk = 1, nbandk
       ib = ix(ibk)
       occ(ib) = bwork(ib)
       occt = occt + occ(ib)
       if( occt+small >= deleud(nspin) ) then
           noccbk = ibk
           occ(ib) = occ(ib) - ( occt - deleud(nspin) )
           exit
       end if
    end do
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       occ(ib) = 0.d0
    end do

if( lfermi.eq.1 .or. nband.eq.nocc ) then
!--- check degeneracy
    toldeg = 1.d-02
    refeig = eig(ix(noccbk))
    ib1 = 1
    do ibk = noccbk - 1, 1, -1
       ib = ix(ibk)
       if( abs(refeig-eig(ib)) > toldeg ) then
           ib1 = ibk + 1
           exit
       end if
    end do
    ib2 = nbandk
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       if( abs(refeig-eig(ib)) > toldeg ) then
           ib2 = ibk - 1
           exit
       end if
    end do
    del = 0.d0
    wei = 0.d0
    do ibk = ib1, ib2
       ib = ix(ibk)
       del = del + occ(ib)
       wei = wei + bwork(ib)
    end do
    del = del/wei
    do ibk = ib1, ib2
       ib = ix(ibk)
       occ(ib) = del*bwork(ib)
    end do
!-----------------------------------------------------------------------
end if
!--- metallic ----------------------------------------------------------
eps   = 1.0d-08
if( lfermi.eq.2 ) then
    eps2   = 1.0d-12
  else
    eps2   = 1.0d-12
    if( licall ) then
        !----- initial set for Gaussian broadening -----
        call setgbr( tbfi ,tbfia, ntbfi, tbxe, tbmx, piru )

        !---the oder of approximant
        npolyorder = lfermi - 3

        !---alocate memory
        allocate( coefdelta(0:npolyorder), hpoly(0:2*npolyorder),  &
& stat=status )

        !---set coefficients An
        coefdelta(0) = 1.d0/sqrt(acos(-1.d0))
        do i = 1, npolyorder
           coefdelta(i) = coefdelta(i-1)*(-1.d0)/4.d0/dble(i)
        end do
    end if
    licall = .false.
end if


enel = deleud(nspin)
zero = 0.d0
one  = 1.d0
half = 0.5d0
!--- estimate Fermi energy by False position method
nocclw = noccbk
nocchg = noccbk + 1
ib = ix(nocclw)
e1 = eig(ib)
ib = ix(nocchg)
e2 = eig(ib)
!---      fermie = 0.5d0*(e1+e2)
itt = 0
do
   itt = itt + 1
   if( itt == 1 ) then
      fermie = e1
   else if( itt == 2 ) then
      fermie = e2
   else if( itt == 3 ) then
      fermie = (e1 + e2)*0.5d0
   else
      fermie = e1 - (e2-e1)/(d2-d1)*d1
   end if

      sum = zero
      if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
          do ib = 1, nbandk
              a = (eig(ib)-fermie)/tfermi
              if(a.gt.3.5d+01) then
                  occ(ib) = zero
              else if(a.lt.-3.3d+01) then
                  occ(ib) = one
              else
                  occ(ib) = one/(one+dexp(a))
              end if
              occ(ib) = occ(ib)*bwork(ib)
              sum = sum + occ(ib)
          end do
        else
!     --- Gaussian broadening ---
          do ib = 1, nbandk
              a = (eig(ib)-fermie)/tfermi
              if(a.gt.tbmx) then
                  occ(ib) = zero
              else if(a.lt.-tbmx) then
                  occ(ib) = one
              else
                  dm = (a+tbmx)/tbxe
                  m  = 0.5d0*dm
                  m  = 2*m
                  d  = 0.5d0*( dm - dble(m) )
                  occ(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
     &                                            - tbfi(m) ) + tbfi(m)
              end if
              occ(ib) = occ(ib)*bwork(ib)
              sum = sum + occ(ib)
          end do

          !---Methfessel & Paxton smearing
          if( npolyorder > 0 ) then
              do ib = 1, nbandk
                 a = (eig(ib)-fermie)/tfermi
                 if( abs(a) > tbmx ) cycle
                 hpoly(0) = 1.d0
                 hpoly(1) = 2.d0*a
                 do j = 2, 2*npolyorder - 1
                    hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
                 end do
                 gfi = 0.d0
                 do j = 1, npolyorder
                    gfi = gfi + coefdelta(j)*hpoly(2*j-1)
                 end do
                 gfi = gfi * exp(-a*a) * bwork(ib)
                 occ(ib) = occ(ib) + gfi
                 sum     = sum     + gfi
              end do
          end if
      end if

      d = sum - enel
!      if( myid == 0 ) then
!          write(*,*) itt, d, sum
!          write(*,*) fermie, enel*half
!      end if
      if( dabs(d) < 1.d-03 ) exit

   if( itt == 1 ) then
       d1 = d
       if( d > 0.d0) then
           if( nocclw == 1 ) then
               ib = ix(nocclw)
               e1 = e1 - abs(eig(ib))*0.5d0
           else
               nocclw = nocclw - 1
               ib = ix(nocclw)
               e1 = eig(ib)
           end if
           itt = 0
       end if
   else if( itt == 2 ) then
       d2 = d
       if( d < 0.d0) then
           if( nocchg == nbandk ) then
               ib = ix(nocchg)
               e2 = e2 + abs(eig(ib))*0.5d0
           else
               nocchg = nocchg + 1
               ib = ix(nocchg)
               e2 = eig(ib)
           end if
           itt = 1
       end if
   else
       if( d1*d < 0.d0 ) then
           e2 = fermie
           d2 = d
         else
           e1 = fermie
           d1 = d
       end if
   end if
end do


if(.not.(dabs(d).lt.eps)) then
!--- calculate fermi energy by newton iteration
    itr =0
    sumv=zero
    nchk=0
    do ib = 1, nbandk
       a = (eig(ib)-fermie)/tfermi
       if( a.lt.-3.3d+01 ) then
           sumv = sumv + one*bwork(ib)
       else if( a.le.3.5d+01 ) then
            nchk=nchk+1
            awork(nchk)=eig(ib)
            bwork(nchk)=bwork(ib)
       end if
    end do
    c = enel - sumv

    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                px =px+one/(one+y)   * bwork(ib)
                pdx=pdx+y/(one+y)**2 * bwork(ib)
            end if
        end do

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

      else
!         --- Gaussian broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.tbmx)) then
                if(a.lt.-tbmx) then
                    px  = px + one * bwork(ib)
                  else
                    dm = (a+tbmx)/tbxe
                    m  = 0.5d0*dm
                    m  = 2*m
                    d  = 0.5d0*( dm - dble(m) )
                    gfi = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
                    px  = px + gfi            * bwork(ib)
                    pdx = pdx+ exp(-a*a)*piru * bwork(ib)
                end if
            end if
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nchk
               a = ( awork(ib) - fermie )/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
                gfi = 0.d0
               dgfi = 0.d0
               do j = 1, npolyorder
                   gfi =  gfi + coefdelta(j)*hpoly(2*j-1)
                  dgfi = dgfi + coefdelta(j)*hpoly(2*j)
               end do
                gfi =  gfi * exp(-a*a) * bwork(ib)
               dgfi = dgfi * exp(-a*a) * bwork(ib)
               px  = px  +  gfi
               pdx = pdx + dgfi
            end do
        end if

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

    end if

!--- check the sum of weighting function
ii = 0
do k = 1, nkpnt
do ib = 1, nband
   ii = ii + 1
   bwork(ii) = wbzk(k)
end do
end do

    sum=zero
    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---
        do ib = 1, nbandk
           a = (eig(ib)-fermie)/tfermi
           if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                occ(ib)= one/(one+y) * bwork(ib)
                sum    = sum + occ(ib)
           else
                occ(ib) = 0.d0
           end if
        end do
      else
!         --- Gaussian broadening ---
        do ib = 1, nbandk
            a = (eig(ib)-fermie)/tfermi
            if(a.gt.tbmx) then
                occ(ib) = zero
            else if(a.lt.-tbmx) then
                occ(ib) = one
            else
                dm = (a+tbmx)/tbxe
                m  = 0.5d0*dm
                m  = 2*m
                d  = 0.5d0*( dm - dble(m) )
                occ(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
            end if
            occ(ib) = occ(ib) * bwork(ib)
            sum = sum + occ(ib)
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nbandk
               a = (eig(ib)-fermie)/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder - 1
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
               gfi = 0.d0
               do j = 1, npolyorder
                  gfi = gfi + coefdelta(j)*hpoly(2*j-1)
               end do
               gfi = gfi * exp(-a*a) * bwork(ib)
               occ(ib) = occ(ib) + gfi
               sum     = sum     + gfi
            end do
        end if
    end if

    d=sum-enel
end if

feneud(nspin) = fermie

!--- small correction
weizzz = 0.d0
do ib = 1, nbandk
   weizzz = weizzz + ( bwork(ib) - occ(ib) )*( occ(ib) - zero )
end do
if( abs(weizzz).gt.1.d-30 ) then
    do ib = 1, nbandk
       wzzz = ( bwork(ib) - occ(ib) )*( occ(ib) - zero )/weizzz
       occ(ib) = occ(ib) - wzzz*d
    end do
end if

if( dabs(d).gt.eps ) then
    if(loutfile(1)) write(nfile(1),1020) sum
    if(loutfile(2)) write(nfile(2),1020) sum
end if

!--- calculate entropy term
tentrpy=zero
if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
    do ib = 1, nbandk
       occib = occ(ib)/bwork(ib)
       if(.not.( occib.lt.eps .or. occib.gt.one-eps)) then
            tentrpy = tentrpy  &
&                  + bwork(ib)*( occib*dlog(occib)  &
&                              + (one-occib)*dlog(one-occib) )
       end if
    end do
    entrpy = entrpy - tfermi*tentrpy
else
!     --- Gaussian broadening ---
    do ib = 1, nbandk
       a = (eig(ib)-fermie)/tfermi
       a = a*a
       if( a.lt.3.3d+01 ) then
           tentrpy = tentrpy + exp(-a) * bwork(ib)
        end if
    end do
    if( npolyorder == 0 ) then
        entrpy = entrpy + tfermi*tentrpy*piru*0.5d0

    else
        !---Methfessel & Paxton smearing
        entrpymp = 0.d0
        do ib = 1, nbandk
           a = (eig(ib)-fermie)/tfermi
           if( abs(a) > tbmx ) cycle
           hpoly(0) = 1.d0
           hpoly(1) = 2.d0*a
           do j = 2, 2*npolyorder
              hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
           end do
           j = npolyorder
           entrpymp = entrpymp + hpoly(2*j) * exp(-a*a) * bwork(ib)
        end do
        j = npolyorder
        entrpy = entrpy + tfermi*coefdelta(j)*entrpymp*0.5d0
    end if
end if

call set_wegud( wegud, occ, nspin, nband, nkpnt )

end if fixif

end do spindo


if( lfixud ) then
    fermie = 0.5d0*( feneud(1) + feneud(2) )
    !--- if( lfermi.eq.1 .or. nband.eq.nocc ) return
    return
end if



!================
nbandk = nbandk*2
!================

ii = 0
do k = 1, nkpnt
do nspin = 1, 2
   do ib = 1, nband
      ii = ii + 1
      awork(ii) = egvlud(ib,nspin,k)
   end do
end do
end do
!--- sorting by eigenvalues : awork --------------------------------------
!     ix     : data pointer of record

ddmx = 0.d0
do ib = 1, nbandk
   ddmx = max( ddmx, abs(awork(ib)) )
end do

fctmax = dkeymx/ddmx
do ib = 1, nbandk
   key(ib) = fctmax*awork(ib)
   ix(ib)  = ib
end do

if( nbandk.gt.1 ) then
    call vsgar( nbandk, key, ix, nbandk, dstkey, x1, x2, x2 )
    call vsoe(  nbandk, key, ix, nbandk )
end if
!-----------------------------------------------------------------------


ii = 0
do k = 1, nkpnt
do nspin = 1, 2
   do ib = 1, nband
      ii = ii + 1
      bwork(ii) = wbzk(k)
   end do
end do
end do


!--- non metallic ------------------------------------------------------
    occt = 0.d0
    noccbk = nbandk
    do ibk = 1, nbandk
       ib = ix(ibk)
       wegud(ib) = bwork(ib)
       occt = occt + wegud(ib)
       if( occt >= dble(nel) ) then
           noccbk = ibk
           wegud(ib) = wegud(ib) - ( occt - dble(nel) )
           exit
       end if
    end do
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       wegud(ib) = 0.d0
    end do

if( lfermi.eq.1 .or. nband.eq.nocc ) then
!--- check degeneracy
    toldeg = 1.d-02
    refeig = awork(ix(noccbk))
    ib1 = 1
    do ibk = noccbk - 1, 1, -1
       ib = ix(ibk)
       if( abs(refeig-awork(ib)) > toldeg ) then
           ib1 = ibk + 1
           exit
       end if
    end do
    ib2 = nbandk
    do ibk = noccbk + 1, nbandk
       ib = ix(ibk)
       if( abs(refeig-awork(ib)) > toldeg ) then
           ib2 = ibk - 1
           exit
       end if
    end do
    del = 0.d0
    wei = 0.d0
    do ibk = ib1, ib2
       ib = ix(ibk)
       del = del + wegud(ib)
       wei = wei + bwork(ib)
    end do
    del = del/wei
    do ibk = ib1, ib2
       ib = ix(ibk)
       wegud(ib) = del*bwork(ib)
    end do
!-----------------------------------------------------------------------
    return
end if

!--- metallic ----------------------------------------------------------
eps   = 1.0d-08
if( lfermi.eq.2 ) then
    eps2   = 1.0d-12
  else
    eps2   = 1.0d-12
    if( licall ) then
        !----- initial set for Gaussian broadening -----
        call setgbr( tbfi ,tbfia, ntbfi, tbxe, tbmx, piru )

        !---the oder of approximant
        npolyorder = lfermi - 3

        !---alocate memory
        allocate( coefdelta(0:npolyorder), hpoly(0:2*npolyorder),  &
& stat=status )

        !---set coefficients An
        coefdelta(0) = 1.d0/sqrt(acos(-1.d0))
        do i = 1, npolyorder
           coefdelta(i) = coefdelta(i-1)*(-1.d0)/4.d0/dble(i)
        end do
    end if
    licall = .false.
end if

enel = nel
zero = 0.d0
one  = 1.d0
half = 0.5d0
!--- estimate Fermi energy by False position method
nocclw = noccbk
nocchg = noccbk + 1
ib = ix(nocclw)
e1 = awork(ib)
ib = ix(nocchg)
e2 = awork(ib)
!---      fermie = 0.5d0*(e1+e2)
itt = 0
do
   itt = itt + 1
   if( itt == 1 ) then
      fermie = e1
   else if( itt == 2 ) then
      fermie = e2
   else
      fermie = e1 - (e2-e1)/(d2-d1)*d1
   end if

      sum = zero
      if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
          do ib = 1, nbandk
              a = (awork(ib)-fermie)/tfermi
              if(a.gt.3.5d+01) then
                  wegud(ib) = zero
              else if(a.lt.-3.3d+01) then
                  wegud(ib) = one
              else
                  wegud(ib) = one/(one+dexp(a))
              end if
              wegud(ib) = wegud(ib)*bwork(ib)
              sum = sum + wegud(ib)
          end do
        else
!     --- Gaussian broadening ---
          do ib = 1, nbandk
              a = (awork(ib)-fermie)/tfermi
              if(a.gt.tbmx) then
                  wegud(ib) = zero
              else if(a.lt.-tbmx) then
                  wegud(ib) = one
              else
                  dm = (a+tbmx)/tbxe
                  m  = 0.5d0*dm
                  m  = 2*m
                  d  = 0.5d0*( dm - dble(m) )
                  wegud(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
     &                                            - tbfi(m) ) + tbfi(m)
              end if
              wegud(ib) = wegud(ib)*bwork(ib)
              sum = sum + wegud(ib)
          end do

          !---Methfessel & Paxton smearing
          if( npolyorder > 0 ) then
              do ib = 1, nbandk
                 a = (awork(ib)-fermie)/tfermi
                 if( abs(a) > tbmx ) cycle
                 hpoly(0) = 1.d0
                 hpoly(1) = 2.d0*a
                 do j = 2, 2*npolyorder - 1
                    hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
                 end do
                 gfi = 0.d0
                 do j = 1, npolyorder
                    gfi = gfi + coefdelta(j)*hpoly(2*j-1)
                 end do
                 gfi = gfi * exp(-a*a) * bwork(ib)
                 wegud(ib) = wegud(ib) + gfi
                 sum       = sum       + gfi
              end do
          end if
      end if

      d = sum - enel
!      if( myid == 0 ) then
!          write(*,*) itt, d, sum
!          write(*,*) fermie, enel*half
!      end if
      if( dabs(d) < 1.d-03 ) exit

   if( itt == 1 ) then
       d1 = d
       if( d > 0.d0) then
           if( nocclw == 1 ) then
               ib = ix(nocclw)
               e1 = e1 - abs(awork(ib))*0.5d0
           else
               nocclw = nocclw - 1
               ib = ix(nocclw)
               e1 = awork(ib)
           end if
           itt = 0
       end if
   else if( itt == 2 ) then
       d2 = d
       if( d < 0.d0) then
           if( nocchg == nbandk ) then
               ib = ix(nocchg)
               e2 = e2 + abs(awork(ib))*0.5d0
           else
               nocchg = nocchg + 1
               ib = ix(nocchg)
               e2 = awork(ib)
           end if
           itt = 1
       end if
   else
       if( d1*d < 0.d0 ) then
           e2 = fermie
           d2 = d
         else
           e1 = fermie
           d1 = d
       end if
   end if
end do


if(.not.(dabs(d).lt.eps)) then
!--- calculate fermi energy by newton iteration
    itr =0
    sumv=zero
    nchk=0
    do ib = 1, nbandk
       a = (awork(ib)-fermie)/tfermi
       if( a.lt.-3.3d+01 ) then
           sumv = sumv + one*bwork(ib)
       else if( a.le.3.5d+01 ) then
            nchk=nchk+1
            awork(nchk)=awork(ib)
            bwork(nchk)=bwork(ib)
       end if
    end do
    c = enel - sumv

    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                px =px+one/(one+y)   * bwork(ib)
                pdx=pdx+y/(one+y)**2 * bwork(ib)
            end if
        end do

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

      else
!         --- Gaussian broadening ---

     do
        itr  =itr+1
        efold=fermie
        px   =zero
        pdx  =zero

        do ib = 1, nchk
            a = ( awork(ib) - fermie )/tfermi
            if(.not.(a.gt.tbmx)) then
                if(a.lt.-tbmx) then
                    px  = px + one * bwork(ib)
                  else
                    dm = (a+tbmx)/tbxe
                    m  = 0.5d0*dm
                    m  = 2*m
                    d  = 0.5d0*( dm - dble(m) )
                    gfi = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
                    px  = px + gfi            * bwork(ib)
                    pdx = pdx+ exp(-a*a)*piru * bwork(ib)
                end if
            end if
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nchk
               a = ( awork(ib) - fermie )/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
                gfi = 0.d0
               dgfi = 0.d0
               do j = 1, npolyorder
                   gfi =  gfi + coefdelta(j)*hpoly(2*j-1)
                  dgfi = dgfi + coefdelta(j)*hpoly(2*j)
               end do
                gfi =  gfi * exp(-a*a) * bwork(ib)
               dgfi = dgfi * exp(-a*a) * bwork(ib)
               px  = px  +  gfi
               pdx = pdx + dgfi
            end do
        end if

        px    =px-c
        pdx   =pdx/tfermi
        !---If gradient is small, the energy gap should be large, and then exit.
        if( dabs(pdx) <= eps2 ) exit
        d     =px/pdx
        fermie=fermie-d
        if( dabs(d) <= eps2 .or. itr > 50 ) exit
     end do

    end if

!--- check the sum of weighting function
ii = 0
do k = 1, nkpnt
do nspin = 1, 2
   do ib = 1, nband
      ii = ii + 1
      bwork(ii) = wbzk(k)
      awork(ii) = egvlud(ib,nspin,k)
   end do
end do
end do

    sum=zero
    if( lfermi.eq.2 ) then
!         --- Fermi broadening ---
        do ib = 1, nbandk
           a = (awork(ib)-fermie)/tfermi
           if(.not.(a.gt.3.5d+01)) then
                if( a.lt.-3.3d+01 ) then
                    y = 0.d0
                  else
                    y = dexp(a)
                end if
                wegud(ib)= one/(one+y) * bwork(ib)
                sum    = sum + wegud(ib)
           else
                wegud(ib) = 0.d0
           end if
        end do
      else
!         --- Gaussian broadening ---
        do ib = 1, nbandk
            a = (awork(ib)-fermie)/tfermi
            if(a.gt.tbmx) then
                wegud(ib) = zero
            else if(a.lt.-tbmx) then
                wegud(ib) = one
            else
                dm = (a+tbmx)/tbxe
                m  = 0.5d0*dm
                m  = 2*m
                d  = 0.5d0*( dm - dble(m) )
                wegud(ib) = d*( (d-1.d0)*tbfia(m) + tbfi(m+2)  &
&                                           - tbfi(m) ) + tbfi(m)
            end if
            wegud(ib) = wegud(ib) * bwork(ib)
            sum = sum + wegud(ib)
        end do

        !---Methfessel & Paxton smearing
        if( npolyorder > 0 ) then
            do ib = 1, nbandk
               a = (awork(ib)-fermie)/tfermi
               if( abs(a) > tbmx ) cycle
               hpoly(0) = 1.d0
               hpoly(1) = 2.d0*a
               do j = 2, 2*npolyorder - 1
                  hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
               end do
               gfi = 0.d0
               do j = 1, npolyorder
                  gfi = gfi + coefdelta(j)*hpoly(2*j-1)
               end do
               gfi = gfi * exp(-a*a) * bwork(ib)
               wegud(ib) = wegud(ib) + gfi
               sum       = sum       + gfi
            end do
        end if
    end if

    d=sum-enel
end if
feneud(1) = fermie
feneud(2) = fermie

!--- small correction
weizzz = 0.d0
do ib = 1, nbandk
  weizzz = weizzz + ( bwork(ib) - wegud(ib) )*( wegud(ib) - zero )
end do
if( abs(weizzz).gt.1.d-30 ) then
    do ib = 1, nbandk
      wzzz = ( bwork(ib) - wegud(ib) )*( wegud(ib) - zero )/weizzz
      wegud(ib) = wegud(ib) - wzzz*d
    end do
end if

if( dabs(d).gt.eps ) then
    if(loutfile(1)) write(nfile(1),1020) sum
    if(loutfile(2)) write(nfile(2),1020) sum
end if
1020     format('*** sub. fermi, error in sum of weighting'/  &
&                  '*** sum=',es20.12)

!--- calculate entropy term
entrpy=zero
if( lfermi.eq.2 ) then
!     --- Fermi broadening ---
    do ib = 1, nbandk
       occib = wegud(ib)/bwork(ib)
       if(.not.( occib.lt.eps .or. occib.gt.one-eps)) then
            entrpy = entrpy  &
&                  + bwork(ib)*( occib*dlog(occib)  &
&                              + (one-occib)*dlog(one-occib) )
       end if
    end do
    entrpy = -tfermi*entrpy
else
!     --- Gaussian broadening ---
    do ib = 1, nbandk
       a = (awork(ib)-fermie)/tfermi
       a = a*a
       if( a.lt.3.3d+01 ) then
           entrpy = entrpy + exp(-a) * bwork(ib)
        end if
    end do
    entrpy = tfermi*entrpy*piru*0.5d0

    !---Methfessel & Paxton smearing
    if( npolyorder > 0 ) then
        entrpymp = 0.d0
        do ib = 1, nbandk
           a = (awork(ib)-fermie)/tfermi
           if( abs(a) > tbmx ) cycle
           hpoly(0) = 1.d0
           hpoly(1) = 2.d0*a
           do j = 2, 2*npolyorder
              hpoly(j) = 2.d0*( a*hpoly(j-1) - dble(j-1)*hpoly(j-2) )
           end do
           j = npolyorder
           entrpymp = entrpymp + hpoly(2*j) * exp(-a*a) * bwork(ib)
        end do
        j = npolyorder
        entrpy = tfermi*coefdelta(j)*entrpymp*0.5d0
    end if
end if


nocc = nband


return
end




subroutine out_eigocc( nfile, myid, nodes,  &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, iunit)
!-----------------------------------------------------------------------
!     output eigenvalues and occupancies
!-----------------------------------------------------------------------
use outfile
use fermi_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nband, nkpnt
real*8  :: eig(nband,nkpnt), occ(nband,nkpnt)
integer :: norder(nband,nkpnt)
real*8  :: bzk(3,nkpnt)
real*8  :: wbzk(nkpnt)
integer :: nstep, jgcycl
real*8  :: fermie
integer :: iunit  ! = 3 or 20


!if( .not.lnoncollinear ) then
    call out_eigocc2( nfile, myid, nodes,  &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, iunit)
!else
!    !-----noncollinear magnetism
!    call ncout_eigocc2( nfile, myid, nodes,  &
!& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, iunit)
!end if


return
end




subroutine out_eigocc2( nfile, myid, nodes,  &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, iunit)
!-----------------------------------------------------------------------
!     output eigenvalues and occupancies
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: nband, nkpnt
real*8  :: eig(nband,nkpnt), occ(nband,nkpnt)
integer :: norder(nband,nkpnt)
real*8  :: bzk(3,nkpnt)
real*8  :: wbzk(nkpnt)
integer :: nstep, jgcycl
real*8  :: fermie
integer :: iunit  ! = 3 or 20

!-----declare local variables
integer :: k, j, ib, jb
character(20) :: form


if( loutfile(1) ) then

    if( nkpnt > 1 ) then
        write(nfile(iunit),'(10i7)') nstep, jgcycl, nband, nkpnt
      else
        write(nfile(iunit),'(10i7)') nstep, jgcycl, nband
    end if
    do k = 1, nkpnt
    if( nkpnt > 1 ) then
        write(nfile(iunit),'(i5,3f10.6,es13.5)')  &
& k, ( bzk(j,k), j = 1, 3 ), wbzk(k)
    end if
    do ib = 1, nband
       jb = norder(ib,k)
       if( occ(jb,k) >= 0.d0 ) then
           form = '(i6,es13.5,f6.3)'
       else
           form = '(i6,es13.5,f7.3)'
       end if
       write(nfile(iunit),form) jb, eig(jb,k), occ(jb,k)/wbzk(k)
    end do
    end do

    if( iunit == 3 ) then
        !--- output Fermi energy
        write(nfile(16),'(2i7,2es13.5)') nstep, jgcycl, fermie
    end if

end if


return
end




subroutine out_eigocc_spin( nfile, myid, nodes,  &
& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk,  &
& nstep, jgcycl, feneud, iunit )
!-----------------------------------------------------------------------
!     output eigenvalues and occupancies
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: nband
integer :: nkpnt
real*8  :: egvlud(nband,2,nkpnt)
real*8  :: wegud(nband,2,nkpnt)
integer :: nsordr(nband,2,nkpnt)
real*8  :: bzk(3,nkpnt)
real*8  :: wbzk(nkpnt)
integer :: nstep, jgcycl
real*8  :: feneud(2)
integer :: iunit

!-----declare local variables
integer :: k, j, ib, jb1, jb2
character(50) :: form
real*8  :: ecorr
integer :: nunocc


if( loutfile(1) ) then

    if( nkpnt > 1 ) then
        write(nfile(iunit),'(10i7)') nstep, jgcycl, nband, nkpnt
      else
        write(nfile(iunit),'(10i7)') nstep, jgcycl, nband
    end if
    do k = 1, nkpnt
    if( nkpnt > 1 ) then
        write(nfile(iunit),'(i5,3f10.6,es13.5)')  &
& k, ( bzk(j,k), j = 1, 3 ), wbzk(k)
    end if

    ecorr  = 0.d0
    nunocc = nband + 1
!    if( iunit == 20 ) call lr_tddft_ecorr( nfile, myid, nodes, ecorr, nunocc )

    do ib = 1, nband
       jb1 = nsordr(ib,1,k)
       jb2 = nsordr(ib,2,k)
       if( wegud(jb1,1,k) >= 0.d0 ) then
           form = '(i5,es13.5,f6.3'
       else
           form = '(i5,es13.5,f7.3'
       end if
       if( wegud(jb2,2,k) >= 0.d0 ) then
           form = form(1:len_trim(form)) //',i5,es13.5,f6.3)'
       else
           form = form(1:len_trim(form)) //',i5,es13.5,f7.3)'
       end if
       if( ib < nunocc ) then
           write(nfile(iunit),form)  &
&                    jb1, egvlud(jb1,1,k), wegud(jb1,1,k)/wbzk(k),  &
&                    jb2, egvlud(jb2,2,k), wegud(jb2,2,k)/wbzk(k)
       else
           write(nfile(iunit),form)  &
&                    jb1, egvlud(jb1,1,k)+ecorr, wegud(jb1,1,k)/wbzk(k),  &
&                    jb2, egvlud(jb2,2,k)+ecorr, wegud(jb2,2,k)/wbzk(k)
       end if
    end do
    end do

    if( iunit == 3 ) then
        !--- output Fermi energy
        write(nfile(16),'(2i7,2es13.5)')  &
&                              nstep, jgcycl, feneud(1), feneud(2)
    end if

end if


return
end




!subroutine out_eigocc_fssh_gsscf( nfile, myid, nodes, &
!& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie)
!!-----------------------------------------------------------------------
!!     output eigenvalues and occupancies
!!-----------------------------------------------------------------------
!implicit none
!
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!integer :: nband, nkpnt
!real*8,  dimension(nband,nkpnt) :: eig, occ
!integer, dimension(nband,nkpnt) :: norder
!real*8,  dimension(3,nkpnt) :: bzk
!real*8,  dimension(nkpnt)   :: wbzk
!integer :: nstep, jgcycl
!real*8  :: fermie
!
!!-----declare local variables
!integer :: myid_, nodes_, nkd_, k, j, ib, jb
!
!
!call get_worldqm( myid_, nodes_ )
!
!if( myid_ == 0 ) then
!
!    if( nkpnt > 1 ) then
!        write(nfile(20),'(10i7)') nstep, jgcycl, nband, nkpnt
!      else
!        write(nfile(20),'(10i7)') nstep, jgcycl, nband
!    end if
!    do k = 1, nkpnt
!    if( nkpnt > 1 ) then
!        write(nfile(20),'(i5,3f10.6,es13.5)')  &
!& k, ( bzk(j,k), j = 1, 3 ), wbzk(k)
!    end if
!    do ib = 1, nband
!       jb = norder(ib,k)
!       write(nfile(20),'(i6,es13.5,f6.3)')  &
!& jb, eig(jb,k), occ(jb,k)/wbzk(k)
!    end do
!    end do
!
!!    !--- output Fermi energy
!!    write(nfile(16),'(2i7,2es13.5)') nstep, jgcycl, fermie
!
!end if
!
!call get_worldkd( myid_, nodes_, nkd_ )
!
!
!return
!end
!
!
!
!
!subroutine out_eigocc_spin_fssh_gsscf( nfile, myid, nodes, &
!& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk, &
!& nstep, jgcycl, feneud )
!!-----------------------------------------------------------------------
!!     output eigenvalues and occupancies
!!-----------------------------------------------------------------------
!implicit none
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!
!integer :: nband
!integer :: nkpnt
!real*8,  dimension(nband,2,nkpnt) :: egvlud
!real*8,  dimension(nband,2,nkpnt) :: wegud
!integer, dimension(nband,2,nkpnt) :: nsordr
!real*8,  dimension(3,nkpnt) :: bzk
!real*8,  dimension(nkpnt) :: wbzk
!integer :: nstep, jgcycl
!real*8,  dimension(2) :: feneud
!
!!-----declare local variables
!integer :: myid_, nodes_, nkd_, k, j, ib, jb1, jb2
!
!
!call get_worldqm( myid_, nodes_ )
!
!if( myid_ == 0 ) then
!
!    if( nkpnt > 1 ) then
!        write(nfile(20),'(10i7)') nstep, jgcycl, nband, nkpnt
!      else
!        write(nfile(20),'(10i7)') nstep, jgcycl, nband
!    end if
!    do k = 1, nkpnt
!    if( nkpnt > 1 ) then
!        write(nfile(20),'(i5,3f10.6,es13.5)')  &
!& k, ( bzk(j,k), j = 1, 3 ), wbzk(k)
!    end if
!    do ib = 1, nband
!       jb1 = nsordr(ib,1,k)
!       jb2 = nsordr(ib,2,k)
!       write(nfile(20),1005)  &
!&                    jb1, egvlud(jb1,1,k), wegud(jb1,1,k)/wbzk(k),  &
!&                    jb2, egvlud(jb2,2,k), wegud(jb2,2,k)/wbzk(k)
!1005            format( 2(i5,es13.5,f6.3) )
!    end do
!    end do
!
!!    !--- output Fermi energy
!!     write(nfile(16),'(2i7,2es13.5)')  &
!!&                              nstep, jgcycl, feneud(1), feneud(2)
!
!end if
!
!call get_worldkd( myid_, nodes_, nkd_ )
!
!
!return
!end




subroutine set_wegud( wegud, occ, nspin, nband, nkpnt )
!-----------------------------------------------------------------------
implicit none
integer :: nspin, nband, nkpnt
real*8,  dimension(nband,2,nkpnt) :: wegud
real*8,  dimension(nband,nkpnt) :: occ

!------declare local variables
integer :: k, i


do k = 1, nkpnt
   do i = 1, nband
      wegud(i,nspin,k) = occ(i,k)
   end do
end do


return
end




subroutine setgbr( tbfi ,tbfia, ntbfi, tbxe, tbmx, piru )
!-----------------------------------------------------------------------
!       initial set for Gaussian broadening
!-----------------------------------------------------------------------
!--- use fermi_variables
implicit none
integer :: ntbfi
real*8,  dimension(0:ntbfi) :: tbfi ,tbfia
real*8 :: tbxe
real*8 :: tbmx
real*8 :: piru

!-----declare local variables
integer :: i
real*8  :: e, gfi, DERFNC
real*8  :: ererr = 1.0D-11
save ererr


piru  = 1.d0/dsqrt(acos(-1.d0))
tbmx  = 5.2d0
tbxe  = tbmx*2.d0/dble(ntbfi)
do i = 0, ntbfi
   e = -tbmx + tbxe*dble(i)
   if( abs(e).lt.ererr ) then
       gfi = 0.5d0
     else
       gfi = 0.5d0*(1.d0 - sign(1.d0,e)*(1.d0-DERFNC(abs(e),ererr)))
   end if
   tbfi(i) = gfi
end do
do i = 0, ntbfi - 2
   tbfia(i) = 2.D0*( tbfi(i) + tbfi(i+2) - 2.d0*tbfi(i+1) )
end do
tbfia(ntbfi-1) = 0.d0
tbfia(ntbfi)   = 0.d0


return
end
