



#ifndef GNU
#ifdef LIBFFTW3
module for_fftw3
!-----------------------------------------------------------------------
! type declaration of variables for electronic bands and plane waves
!-----------------------------------------------------------------------
use, intrinsic :: iso_c_binding
implicit none
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: cin, cout
real(C_DOUBLE), allocatable, dimension(:) :: rwork

integer :: lplnfftkd = 0
integer :: lplkfft0  = 0

save
end module




subroutine fftw3_alloc( nfile, myid, nodes,  &
& alloc_mem, nfftkd, pwscale, lgamma, lnoncollinear, kfft0 )
!-----------------------------------------------------------------------
!     allocate memory for the plane wave method
!-----------------------------------------------------------------------
use for_fftw3
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: nfftkd  ! array size for FFT variables
real*8  :: pwscale
logical :: lgamma, lnoncollinear
integer :: kfft0

!-----declare local variables
!integer :: lplnfftkd = 0
!integer :: lplkfft0  = 0
real*8  :: the_mem
integer :: status
!save lplnfftkd, lplkfft0


if( .not.lgamma .or. lnoncollinear ) then

    if(  &
&       nfftkd > lplnfftkd  &
&  .or. kfft0  > lplkfft0   ) then
        !-----if already allocated, deallocate arrays
        if( allocated(cin) ) then
            the_mem = 16.d0 * ( size(cin) + size(cout) ) + 8.d0 * size(rwork)
            !-----deallocate arrays
            deallocate( cin, cout, rwork, stat=status )
            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'fftw3_alloc', .true. )
        end if

        !-----allocate arrays
        lplnfftkd = nfftkd * pwscale
        lplkfft0  = kfft0  * pwscale
        allocate(  &
& cin(lplkfft0), cout(lplkfft0), rwork(lplnfftkd),  &
& stat=status )

        the_mem = 16.d0 * ( size(cin) + size(cout) ) + 8.d0 * size(rwork)

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'fftw3_alloc', .true. )
    end if

else

    if(  &
&       nfftkd  > lplnfftkd  ) then
        !-----if already allocated, deallocate arrays
        if( allocated(rwork) ) then
            the_mem = 8.d0 * size(rwork)
            !-----deallocate arrays
            deallocate( rwork, stat=status )
            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'fftw3_alloc', .true. )
        end if

        !-----allocate arrays
        lplnfftkd  = nfftkd * pwscale
        allocate(  &
& rwork(lplnfftkd),  &
& stat=status )

        the_mem = 8.d0 * size(rwork)

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'fftw3_alloc', .true. )
    end if

end if


return
end subroutine




subroutine rfft3( key, fft3x, fft3y, work,                         &
&                      kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    FFTW3
!
!         real-to-complex transformation
!
!    compile :
!         f77 (options) fft3r.F90 -L$HOME/FFTW3/lib -lfftw3 -lm
!
!
!       key = 1    : real    to complex transformation
!       key = else : complex to real    transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!
! See http://www.fftw.org/fftw3_doc/Introduction.html#Introduction
! ... the standard FFTW distribution works most efficiently for arrays 
! whose size can be factored into small primes (2, 3, 5, and 7), and 
! otherwise it uses a slower general-purpose routine. If you need efficient 
! transforms of other sizes, you can use FFTW's code generator, ...
!
!-----------------------------------------------------------------------
use, intrinsic :: iso_c_binding
use for_fftw3
implicit none
include 'fftw3.f03'
integer :: key
real*8  :: fft3x(*), fft3y(*)
complex*16 :: work(*)
integer :: kfft1, kfft2, kfft3, kfft0, ierror
!--- for FFTW3 ----------------------------------------------------------
type(C_PTR) :: plan
!complex(C_DOUBLE_COMPLEX) :: work(nfftk)
!real(C_DOUBLE)  :: 
!-----------------------------------------------------------------------

!-----declare local variables
integer :: kfft1h, ndt, nd1, nd2, nd3, ijk, ijk2, ijk3
integer :: k, md1, md2, md3, mijk, mijk2, mijk3
real*8  :: fintot


ierror = 0

kfft1h = kfft1/2+1
if(key.eq.1) then
   plan = fftw_plan_dft_r2c_3d(kfft3,kfft2,kfft1, rwork, work, FFTW_ESTIMATE)
else !if(key.eq.2) then
   plan = fftw_plan_dft_c2r_3d(kfft3,kfft2,kfft1, work, rwork, FFTW_ESTIMATE)
end if

if(key.eq.1) then
   rwork(1:kfft0) = fft3x(1:kfft0)
   call fftw_execute_dft_r2c(plan, rwork, work)

   fintot=1.d0/dble(kfft0)
   do k=1,kfft1h*kfft2*kfft3
      work(k)=work(k)*fintot
   end do

   ndt = 0
   do nd3 = 0, kfft3 - 1
      ijk3 = nd3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      ijk2 = nd2*kfft1
   do nd1 = 0, kfft1h - 1
      ndt = ndt + 1
      ijk = 1+nd1+ijk2+ijk3
      fft3x(ijk) =  dble( work(ndt) )
      fft3y(ijk) = dimag( work(ndt) )
   enddo
   enddo
   enddo
   do nd3 = 0, kfft3 - 1
      md3 = kfft3 - nd3
      if( md3.eq.kfft3 ) md3 = 0
      ijk3  = nd3*kfft2*kfft1
      mijk3 = md3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      md2 = kfft2 - nd2
      if( md2.eq.kfft2 ) md2 = 0
      ijk2  = nd2*kfft1
      mijk2 = md2*kfft1
   do nd1 = kfft1h, kfft1 - 1
      md1 = kfft1 - nd1
      ijk  = 1+nd1+ ijk2+ ijk3
      mijk = 1+md1+mijk2+mijk3
      fft3x(ijk) =  fft3x(mijk)
      fft3y(ijk) = -fft3y(mijk)
   enddo
   enddo
   enddo
else !if(key.eq.2) then
   ndt = 0
   do nd3 = 0, kfft3 - 1
      ijk3 = nd3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      ijk2 = nd2*kfft1
   do nd1 = 0, kfft1h - 1
      ndt = ndt + 1
      ijk = 1+nd1+ijk2+ijk3
      work(ndt) = dcmplx(fft3x(ijk),fft3y(ijk))
   enddo
   enddo
   enddo

   call fftw_execute_dft_c2r(plan, work, rwork)
   fft3x(1:kfft0) = rwork(1:kfft0)
   fft3y(1:kfft0) = 0.d0
end if


!-----deallocate memory
call fftw_destroy_plan(plan)


return
end




subroutine fft3( key, fft3x, fft3y, work,                          &
&                     kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    FFTW
!
!         complex-to-complex transformation
!
!    compile :
!         f77 (options) fft3.F90 -L$HOME/FFTW3/lib -lfftw3 -lm
!
!
!       key = 1    : forward  complex to complex transformation
!       key = else : backward complex to complex transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!
!
! See http://www.fftw.org/fftw3_doc/Introduction.html#Introduction
! ... the standard FFTW distribution works most efficiently for arrays 
! whose size can be factored into small primes (2, 3, 5, and 7), and 
! otherwise it uses a slower general-purpose routine. If you need efficient 
! transforms of other sizes, you can use FFTW's code generator, ...
!
!-----------------------------------------------------------------------
use, intrinsic :: iso_c_binding
use for_fftw3
implicit none
include 'fftw3.f03'
integer :: key
real*8  :: fft3x(*), fft3y(*)
complex*16 :: work(*)
integer :: kfft1, kfft2, kfft3, kfft0, ierror
!--- for FFTW3 ----------------------------------------------------------
type(C_PTR) :: plan
!complex(C_DOUBLE_COMPLEX) :: work(nfftk)
!real(C_DOUBLE)  :: 
!-----------------------------------------------------------------------

!-----declare local variables
real*8  :: fintot
real*8  :: the_mem
integer :: status


if( kfft0  > lplkfft0   ) then
    !-----if already allocated, deallocate arrays
    if( allocated(cin) ) then
        the_mem = 16.d0 * ( size(cin) + size(cout) )
        !-----deallocate arrays
        deallocate( cin, cout, stat=status )
        !------error trap
!        call check_dealloc_accum( nfile, myid, nodes,  &
!& status, the_mem, 'fft3', .true. )
    end if

    !-----allocate arrays
    lplkfft0  = kfft0
    allocate( cin(lplkfft0), cout(lplkfft0), stat=status )

    the_mem = 16.d0 * ( size(cin) + size(cout) )

    !------error trap
!    call check_alloc_accum( nfile, myid, nodes,  &
!& status, the_mem, 'fft3', .true. )
end if


ierror = 0

if(key.eq.1) then
   plan = fftw_plan_dft_3d(kfft3,kfft2,kfft1, cin, cout, FFTW_FORWARD, FFTW_ESTIMATE)
else !if(key.eq.2) then
   plan = fftw_plan_dft_3d(kfft3,kfft2,kfft1, cin, cout, FFTW_BACKWARD, FFTW_ESTIMATE)
end if


cin(1:kfft0) = dcmplx(fft3x(1:kfft0),fft3y(1:kfft0))

call fftw_execute_dft(plan, cin, cout)


!-----deallocate memory
call fftw_destroy_plan(plan)


if(key.eq.1) then
   fintot=1.d0/dble(kfft0)
   cout(1:kfft0)=cout(1:kfft0)*fintot
end if

fft3x(1:kfft0) =  dble( cout(1:kfft0) )
fft3y(1:kfft0) = dimag( cout(1:kfft0) )


return
end
#else




#ifdef LIBFFTW
subroutine rfft3( key, fft3x, fft3y, work,                         &
&                      kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    FFTW
!
!         real-to-complex transformation
!
!    compile :
!         f77 (options)   fft3r.F -L$HOME/FFTW/lib -lrfftw -lfftw -lm
!             -DPOINTER64   <-  if necessary.
!
!
!       key = 1    : real    to complex transformation
!       key = else : complex to real    transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!
!
!  FFTW is best at handling sizes of the form 
!  2^a 3^b 5^c 7^d 11^e 13^f, where e+f is either 0 or 1, 
!  and the other exponents are arbitrary. Other sizes are computed by
!  means of a slow, general-purpose routine (which nevertheless retains 
!  O(n lg n) performance, even for prime sizes). (It is possible to 
!  customize FFTW for different array sizes. See Section Installation 
!  and Customization, for more information.)
!  Transforms whose sizes are powers of 2 are especially fast. 
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)
!--- for FFTW ----------------------------------------------------------
include 'fftw_f77.i'
#ifdef POINTER64
integer*8    iplannml,iplaninv
integer*8    iwisnml,iwisinv
integer*8    iplan,iwisdom
#else
integer      iplannml,iplaninv
integer      iwisnml,iwisinv
integer      iplan,iwisdom
#endif
save         iplannml, iplaninv
save         iwisnml, iwisinv
!-----------------------------------------------------------------------


ierror = 0

kfft1h = kfft1/2+1
if(key.eq.1) then
   call rfftw3d_f77_create_plan(iplannml,kfft1,kfft2,kfft3,        &
&                      FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)
   if (iplannml .eq. 0) then
       ierror = 1
       return
!             write(6,*) 'FFTW: plan allocate error(nml)'
!             stop
   end if
   iplan   = iplannml
!         iwisdom = iwisnml
else if(key.eq.2) then
   call rfftw3d_f77_create_plan(iplaninv,kfft1,kfft2,kfft3,        &
&                     FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)
   if (iplaninv .eq. 0) then
       ierror = 1
       return
!             write(6,*) 'FFTW: plan allocate error(inv)'
!             stop
   end if
   iplan   = iplaninv
!         iwisdom = iwisinv
end if

if(key.eq.1) then
   call rfftwnd_f77_one_real_to_complex(iplan,fft3x,work)
else if(key.eq.2) then
   ndt = 0
   do nd3 = 0, kfft3 - 1
      ijk3 = nd3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      ijk2 = nd2*kfft1
   do nd1 = 0, kfft1h - 1
      ndt = ndt + 1
      ijk = 1+nd1+ijk2+ijk3
      work(ndt) = dcmplx(fft3x(ijk),fft3y(ijk))
   enddo
   enddo
   enddo

   call rfftwnd_f77_one_complex_to_real(iplan,work,fft3x)
end if


!-----deallocate memory
call rfftwnd_f77_destroy_plan(iplan)


if(key.eq.1) then
   fintot=1.d0/dble(kfft0)
   do k=1,kfft1h*kfft2*kfft3
      work(k)=work(k)*fintot
   end do

   ndt = 0
   do nd3 = 0, kfft3 - 1
      ijk3 = nd3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      ijk2 = nd2*kfft1
   do nd1 = 0, kfft1h - 1
      ndt = ndt + 1
      ijk = 1+nd1+ijk2+ijk3
      fft3x(ijk) =  dble( work(ndt) )
      fft3y(ijk) = dimag( work(ndt) )
   enddo
   enddo
   enddo
   do nd3 = 0, kfft3 - 1
      md3 = kfft3 - nd3
      if( md3.eq.kfft3 ) md3 = 0
      ijk3  = nd3*kfft2*kfft1
      mijk3 = md3*kfft2*kfft1
   do nd2 = 0, kfft2 - 1
      md2 = kfft2 - nd2
      if( md2.eq.kfft2 ) md2 = 0
      ijk2  = nd2*kfft1
      mijk2 = md2*kfft1
   do nd1 = kfft1h, kfft1 - 1
      md1 = kfft1 - nd1
      ijk  = 1+nd1+ ijk2+ ijk3
      mijk = 1+md1+mijk2+mijk3
      fft3x(ijk) =  fft3x(mijk)
      fft3y(ijk) = -fft3y(mijk)
   enddo
   enddo
   enddo
else
   do k=1,kfft0
!cc            fft3x(k) = fft3x(k)*dble(kfft0)
      fft3y(k) = 0.d0
   end do
end if


return
end




subroutine fft3( key, fft3x, fft3y, work,                          &
&                     kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    FFTW
!
!         complex-to-complex transformation
!
!    compile :
!         f77 (options)   fft3.F -L$HOME/FFTW/lib -lrfftw -lfftw -lm
!             -DPOINTER64   <-  if necessary.
!
!
!       key = 1    : forward  complex to complex transformation
!       key = else : backward complex to complex transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!
!
!  FFTW is best at handling sizes of the form 
!  2^a 3^b 5^c 7^d 11^e 13^f, where e+f is either 0 or 1, 
!  and the other exponents are arbitrary. Other sizes are computed by
!  means of a slow, general-purpose routine (which nevertheless retains 
!  O(n lg n) performance, even for prime sizes). (It is possible to 
!  customize FFTW for different array sizes. See Section Installation 
!  and Customization, for more information.)
!  Transforms whose sizes are powers of 2 are especially fast. 
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)
!--- for FFTW ----------------------------------------------------------
include 'fftw_f77.i'
#ifdef POINTER64
integer*8    iplannml,iplaninv
integer*8    iwisnml,iwisinv
integer*8    iplan,iwisdom
#else
integer      iplannml,iplaninv
integer      iwisnml,iwisinv
integer      iplan,iwisdom
#endif
save         iplannml, iplaninv
save         iwisnml, iwisinv
!-----------------------------------------------------------------------


ierror = 0

kft0tm = kfft0


if(key.eq.1) then
   call fftw3d_f77_create_plan(iplannml,kfft1,kfft2,kfft3,  &
&                      FFTW_FORWARD,FFTW_ESTIMATE + FFTW_IN_PLACE)
   if (iplannml .eq. 0) then
       ierror = 1
       return
!            write(6,*) 'FFTW: plan allocate error(nml)'
!            stop
   end if
   iplan   = iplannml
!         iwisdom = iwisnml
else if(key.eq.2) then
   call fftw3d_f77_create_plan(iplaninv,kfft1,kfft2,kfft3,  &
&                     FFTW_BACKWARD,FFTW_ESTIMATE + FFTW_IN_PLACE)
   if (iplaninv .eq. 0) then
       ierror = 1
       return
!            write(6,*) 'FFTW: plan allocate error(inv)'
!            stop
   end if
   iplan   = iplaninv
!         iwisdom = iwisinv
end if



do k = 1, kft0tm
   work(k) = dcmplx(fft3x(k),fft3y(k))
end do


call fftwnd_f77_one(iplan,work,0)


!-----deallocate memory
call fftwnd_f77_destroy_plan(iplan)


if(key.eq.1) then
   fintot=1.d0/dble(kft0tm)
   do k=1,kft0tm
      work(k)=work(k)*fintot
   end do
end if

do k=1,kft0tm
    fft3x(k) =  dble( work(k) )
    fft3y(k) = dimag( work(k) )
end do


!c--- check
!c         DO nd1 = 0, kfft1 - 1
!            nd1 = 0
!            K1 = 1 + kfft2*kfft3*nd1
!c         DO nd2 = 0, kfft2 - 1
!            nd2 = 0
!            K2 = K1 + kfft3*nd2
!         DO nd3 = 0, kfft3 - 1
!            M = K2 + nd3
!            ijk = 1+nd1+(nd2+nd3*kfft2)*kfft1
!            WRITE(*,1000) fft3x(ijk), fft3y(ijk)
! 1000    FORMAT(3(2X,2F16.10))
!      enddo
!c      enddo
!c      enddo



return
end
#else




#ifdef SSL2VP
subroutine rfft3( key, fft3x, fft3y, work,                         &
&                      kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    Real-to-Complex 3D Fast Fourier Transformation
!
!    SSL2 library for VPP
!
!       key = 1    : real    to complex transformation
!       key = else : complex to real    transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)


ierror = 0

call rfft3a( key, fft3x, fft3y, work(1), work(1),                  &
&            work((kfft0+2*kfft2*kfft3)/2+1),                      &
&                 kfft1, kfft2, kfft3, kfft0, ierror )

return
end




subroutine rfft3a( key, fft3x, fft3y, datar, datai, work1,         &
&                      kft1tm, kft2tm, kft3tm, kft0tm, ierror )
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
dimension datar(kft1tm+2,kft2tm,kft3tm)
dimension datai((kft1tm+2)/2,kft2tm,kft3tm,2)
dimension work1(kft1tm+2,kft2tm,kft3tm)
!cc      equivalence (datar,datai)


if( key.eq.1 ) then
    isign =  1
  else
    isign = -1
endif

nl = kft1tm
nm = kft2tm
no = kft3tm
kfft1 = kft1tm
kfft2 = kft2tm
kfft3 = kft3tm
kfft0 = kft0tm

if( key.eq.1 ) then
    do k = 1,no
    do j = 1,nm
    do i = 1,nl
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       datar(i,j,k) = fft3x(ijk)
    enddo
    enddo
    enddo
  else
    do k=1,no
    do j=1,nm
    do i=1,nl/2 + 1
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       datai(i,j,k,1) = fft3x(ijk)
       datai(i,j,k,2) = fft3y(ijk)
    enddo
    enddo
    enddo
endif

call dvrpf3(datar, nl+2,nm,no,isign,work1, icon)

if( key.eq.1 ) then
    do k=1,no
    do j=1,nm
    do i=1,nl/2+1
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       fft3x(ijk) = datai(i,j,k,1)
       fft3y(ijk) = datai(i,j,k,2)
    enddo
    enddo
    enddo
    do k=1,no
    do j=1,nm
       nd2 = kfft2 - (j-1)
       if( nd2.eq.kfft2 ) nd2 = 0
       nd3 = kfft3 - (k-1)
       if( nd3.eq.kfft3 ) nd3 = 0
    do i=2,nl/2
       nd1 = kfft1 - (i-1)
       ijk = 1+nd1+(nd2+nd3*kfft2)*kfft1
       fft3x(ijk) =  datai(i,j,k,1)
       fft3y(ijk) = -datai(i,j,k,2)
    enddo
    enddo
    enddo
    fintot=1.d0/dble(nl*nm*no)
    do ijk = 1, kfft0
       fft3x(ijk) = fft3x(ijk) * fintot
       fft3y(ijk) = fft3y(ijk) * fintot
    enddo
  else
    do k=1,no
    do j=1,nm
    do i=1,nl
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       fft3x(ijk) = datar(i,j,k)
       fft3y(ijk) = 0.d0
    enddo
    enddo
    enddo
endif

return
end




subroutine fft3( key, fft3x, fft3y, work,                         &
&                     kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transformation
!
!    SSL2 library for VPP
!
!       key = 1    : forward  complex to complex transformation
!       key = else : backward complex to complex transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)


ierror = 0

if( mod(kfft0,2) == 0 ) then
    kfft0h = kfft0/2
  else
    kfft0h = ( kfft0 + 1 )/2
end if
call fft3a( key, fft3x, fft3y, work(1), work(kfft0h+1),  &
&                kfft1, kfft2, kfft3, kfft0, ierror )

return
end




subroutine fft3a( key, fft3x, fft3y, work1, work2,  &
&                      kft1tm, kft2tm, kft3tm, kft0tm, ierror )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transform for wave functions
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(*), fft3y(*)
dimension work1(kft1tm,kft2tm,kft3tm)
dimension work2(kft1tm,kft2tm,kft3tm)

if( key.eq.1 ) then
    isign =  1
  else
    isign = -1
endif

nl = kft1tm
nm = kft2tm
no = kft3tm
kfft1 = kft1tm
kfft2 = kft2tm
kfft3 = kft3tm
kfft0 = kft0tm

!          do k = 1,no
!          do j = 1,nm
!          do i = 1,nl
!             nd1 = i-1
!             nd2 = j-1
!             nd3 = k-1
!             ijk = 1+nd1+(nd2+nd3*kfft2)*kfft1
!             datar(i,j,k) = fft3x(ijk)
!             datai(i,j,k) = fft3y(ijk)
!          enddo
!          enddo
!          enddo
!
!      call dvcpf3(datar, datai, kfft1,kfft2,kfft3,isign,
!     &            work1, work2,icon)
!
!          do k=1,no
!          do j=1,nm
!          do i=1,nl
!             nd1 = i-1
!             nd2 = j-1
!             nd3 = k-1
!             ijk = 1+nd1+(nd2+nd3*kfft2)*kfft1
!             fft3x(ijk) = datar(i,j,k)
!             fft3y(ijk) = datai(i,j,k)
!          enddo
!          enddo
!          enddo

call dvcpf3(fft3x, fft3y, kfft1,kfft2,kfft3,isign,  &
&           work1, work2,icon)

if( key.eq.1 ) then
    fintot=1.d0/dble(nl*nm*no)
    do ijk = 1, kfft0
       fft3x(ijk) = fft3x(ijk) * fintot
       fft3y(ijk) = fft3y(ijk) * fintot
    enddo
endif

return
end
#else




#ifdef ASLSX
subroutine rfft3( key, fft3x, fft3y, work,                         &
&                      kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    Real-to-Complex 3D Fast Fourier Transformation
!
!    ASL library for SX
!
!       key = 1    : real    to complex transformation
!       key = else : complex to real    transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)


ierror = 0

!if( mod(kfft1,2) /= 0 ) then
!    lx = kfft1 + 1
!    ly = kfft2
!    lz = kfft3
!  else
!    lx = kfft1 + 2
!    ly = kfft2
!    lz = kfft3
!end if

!--      kfftlx is even    --
!        & kfftlx/2, kfftly, kfftlz = odd --
if( mod(kfft1,4) == 0 ) then
    lx = kfft1 + 2
  else if( mod(kfft1,4) == 1 ) then
    lx = kfft1 + 1
  else if( mod(kfft1,4) == 2 ) then
    lx = kfft1 + 4
  else
    lx = kfft1 + 3
endif
if( mod(kfft2,2) == 0 ) then
    ly = kfft2 + 1
  else
    ly = kfft2
endif
if( mod(kfft3,2) == 0 ) then
    lz = kfft3 + 1
  else
    lz = kfft3
endif

call rfft3a( key, fft3x, fft3y, work(1),  &
& work(lx*ly*lz/2+1), work(lx*ly*lz+1),  &
& kfft1, kfft2, kfft3, kfft0, lx, ly, lz, ierror )

return
end




subroutine rfft3a( key, fft3x, fft3y, datar, wk, trigs,  &
& kft1tm, kft2tm, kft3tm, kft0tm, lx, ly, lz, ierror )
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(*), fft3y(*)
dimension datar(lx,ly,lz)
dimension wk(lx*ly*lz)
dimension trigs(kft1tm+2*(kft2tm+kft3tm))
dimension ifax(60)


!-- initialization of datar --
datar(1:lx,1:ly,1:lz) = 0.d0


if( key.eq.1 ) then
    isign =  1
  else
    isign = -1
endif

nl = kft1tm
nm = kft2tm
no = kft3tm
kfft1 = kft1tm
kfft2 = kft2tm
kfft3 = kft3tm
kfft0 = kft0tm

if( key.eq.1 ) then
    do k = 1,no
    do j = 1,nm
    do i = 1,nl
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       datar(i,j,k) = fft3x(ijk)
    enddo
    enddo
    enddo
  else
    do k=1,no
    do j=1,nm
    do i=1,nl/2 + 1
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       datar(2*i-1,j,k) = fft3x(ijk)
       datar(2*i,  j,k) = fft3y(ijk)
    enddo
    enddo
    enddo
endif

call dfr3fb(nl,nm,no,datar,lx,ly,lz,isign,ifax,trigs,wk,ierr)

if( key.eq.1 ) then
    do k=1,no
    do j=1,nm
    do i=1,nl/2+1
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       fft3x(ijk) = datar(2*i-1,j,k)
       fft3y(ijk) = datar(2*i,  j,k)
    enddo
    enddo
    enddo
    do k=1,no
    do j=1,nm
       nd2 = kfft2 - (j-1)
       if( nd2.eq.kfft2 ) nd2 = 0
       nd3 = kfft3 - (k-1)
       if( nd3.eq.kfft3 ) nd3 = 0
    do i=2,nl/2
       nd1 = kfft1 - (i-1)
       ijk = 1+nd1+(nd2+nd3*kfft2)*kfft1
       fft3x(ijk) =  datar(2*i-1,j,k)
       fft3y(ijk) = -datar(2*i,  j,k)
    enddo
    enddo
    enddo
    fintot=1.d0/dble(nl*nm*no)
    do ijk = 1, kfft0
       fft3x(ijk) = fft3x(ijk) * fintot
       fft3y(ijk) = fft3y(ijk) * fintot
    enddo
  else
    do k=1,no
    do j=1,nm
    do i=1,nl
       ijk = i+((j-1)+(k-1)*kfft2)*kfft1
       fft3x(ijk) = datar(i,j,k)
       fft3y(ijk) = 0.d0
    enddo
    enddo
    enddo
endif

return
end




subroutine fft3( key, fft3x, fft3y, work,                         &
&                     kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transformation
!
!    ASL library for SX
!
!       key = 1    : forward  complex to complex transformation
!       key = else : backward complex to complex transformation
!
!       ierror = 0    : normal return
!       ierror = else : error
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)

ierror = 0

!if( mod(kfft0,2) == 0 ) then
!    kfft0h = kfft0/2
!  else
!    kfft0h = ( kfft0 + 1 )/2
!end if

if( mod(kfft1,2) == 1 .and. mod(kfft2,2) == 1  &
&                     .and. mod(kfft3,2) == 1) then
    call fft3a( key, fft3x, fft3y, work, work(kfft1*kfft2*kfft3+1),  &
&     kfft1, kfft2, kfft3, kfft0, kfft1, kfft2, kfft3, ierror )
else
    if( mod(kfft1,2) == 1 ) then
        lx = kfft1
    else
        lx = kfft1+1
    end if
    if( mod(kfft2,2) == 1 ) then
        ly = kfft2
    else
        ly = kfft2+1
    end if
    if( mod(kfft3,2) == 1 ) then
        lz = kfft3
    else
        lz = kfft3+1
    end if
    kfft0h = lx*ly*lz
    if( mod(kfft0h,2) == 0 ) then
        kfft0h = kfft0h/2
      else
        kfft0h = ( kfft0h + 1 )/2
    end if

    call fft3aw( fft3x, fft3y, work(1), work(kfft0h+1),  &
&          kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
    call fft3awr( fft3x, fft3y, work(1), work(kfft0h+1),  &
&          kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
    call fft3a( key, fft3x, fft3y, work, work(lx*ly*lz+1),  &
&          kfft1, kfft2, kfft3, kfft0, lx, ly, lz, ierror )
    call fft3awr( work(1), work(kfft0h+1), fft3x, fft3y,  &
&          kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
    call fft3aw2( fft3x, fft3y, work(1), work(kfft0h+1),  &
&          kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
end if

return
end




subroutine fft3aw( fft3x, fft3y, work1, work2,  &
&    kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transform for wave functions
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(kfft1,kfft2,kfft3), fft3y(kfft1,kfft2,kfft3)
dimension work1(lx,ly,lz)
dimension work2(lx,ly,lz)


    do k = 1,kfft3
    do j = 1,kfft2
    do i = 1,kfft1
       work1(i,j,k) = fft3x(i,j,k)
       work2(i,j,k) = fft3y(i,j,k)
    enddo
    enddo
    enddo

return
end




subroutine fft3aw2( fft3x, fft3y, work1, work2,  &
&    kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transform for wave functions
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(kfft1,kfft2,kfft3), fft3y(kfft1,kfft2,kfft3)
dimension work1(lx,ly,lz)
dimension work2(lx,ly,lz)


    do k = 1,kfft3
    do j = 1,kfft2
    do i = 1,kfft1
       fft3x(i,j,k) = work1(i,j,k)
       fft3y(i,j,k) = work2(i,j,k)
    enddo
    enddo
    enddo

return
end




subroutine fft3awr( fft3x, fft3y, work1, work2,  &
&    kfft1, kfft2, kfft3, kfft0, lx, ly, lz )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transform for wave functions
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(lx,ly,lz), fft3y(lx,ly,lz)
dimension work1(lx,ly,lz)
dimension work2(lx,ly,lz)


    do k = 1,kfft3
    do j = 1,kfft2
    do i = 1,kfft1
       fft3x(i,j,k) = work1(i,j,k)
       fft3y(i,j,k) = work2(i,j,k)
    enddo
    enddo
    enddo

return
end




subroutine fft3a( key, fft3x, fft3y, work, trigs, &
&    kfft1, kfft2, kfft3, kfft0, lx, ly, lz, ierror )
!-----------------------------------------------------------------------
!    Complex-to-Complex 3D Fast Fourier Transform for wave functions
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension fft3x(lx,ly,lz), fft3y(lx,ly,lz)
dimension work(2*lx*ly*lz)
dimension ifax(60)
dimension trigs(2*(kfft1+kfft2+kfft3))

if( key.eq.1 ) then
    isign =  1
  else
    isign = -1
endif

call dfc3fb(kfft1,kfft2,kfft3,fft3x,fft3y, &
&           lx,ly,lz,isign,ifax,trigs, work, ierror)

if( key.eq.1 ) then
    fintot=1.d0/dble(kfft1*kfft2*kfft3)
    do k = 1,kfft3
    do j = 1,kfft2
    do i = 1,kfft1
       fft3x(i,j,k) = fft3x(i,j,k) * fintot
       fft3y(i,j,k) = fft3y(i,j,k) * fintot
    enddo
    enddo
    enddo
endif

return
end
#else
#endif
#endif
#endif
#endif
#else
!---for GNU
subroutine rfft3( key, fft3x, fft3y, work,                         &
&                      kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)

return
end




subroutine fft3( key, fft3x, fft3y, work,                          &
&                     kfft1, kfft2, kfft3, kfft0, ierror )
!-----------------------------------------------------------------------
implicit real*8(a-h, o-z)
dimension  fft3x(*), fft3y(*)
complex*16 work(*)

return
end
#endif
