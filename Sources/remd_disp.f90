



module displacement_variables
!-----------------------------------------------------------------------
!  Variables for calculation with given atomic displacement
!-----------------------------------------------------------------------
implicit none

real(8) :: displacement = 0d0
integer :: nofconfig = 0
integer, allocatable :: nofatom(:), idofatom(:,:), disp(:,:,:)
real(8), allocatable :: dispofatom(:,:,:)

real(8), allocatable :: x0(:,:)

integer :: counter = 0

save

end module




subroutine set_nstop_from_remd_disp( nfile, myid, nodes, nstop )
!-----------------------------------------------------------------------
! calculation with given atomic displacement
!-----------------------------------------------------------------------
use displacement_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstop

nstop = nofconfig

return
end subroutine




subroutine out_dispdata( nfile )
!-----------------------------------------------------------------------
! out data
!-----------------------------------------------------------------------
use displacement_variables
implicit none
integer :: nfile

!-----declare local variables
integer :: i, n
character(30) :: fmt


write(nfile,'(1x,a,i6)')    '  # of atomic configurations :', nofconfig
write(nfile,'(1x,a,f11.6)') '  displacement in a.u.       :', displacement

write(nfile,'(8x,a,19x,a)') 'atom id', '   displacements in a.u.'
do i = 1, nofconfig
   write(nfile,'(i7,a)',advance='no') i, ':'
   do n = 1, nofatom(i)
      if( n == 1 ) then
          fmt = '(i7,a,3i5,a,3f11.6)'
      else
          fmt = '(8x,i7,a,3i5,a,3f11.6)'
      end if
      write(nfile,fmt) idofatom(n,i), ' (', disp(1:3,n,i), ' )', dispofatom(1:3,n,i)
   end do
end do


return
end subroutine




subroutine updtdisp( nfile, myid, nodes,  &
& x, is, n, h )
!-----------------------------------------------------------------------
!     Coordinate update
!-----------------------------------------------------------------------
use displacement_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: n
real*8,  dimension(3,n) :: x
integer, dimension(n) :: is
real*8,  dimension(3,3) :: h

!-----declare local variables
integer :: na, i, ix, ii
real(8) :: hi(3,3), vtemp, dsp(3)
real(8) :: the_mem
integer :: status


counter = counter + 1
if( counter == 1 ) then
    !---allocate memory
    allocate( x0(3,n), stat=status )
    the_mem = 8d0*( size(x0) )

    !------error trap
    call check_alloc_accum( nfile, myid, nodes, status, the_mem, 'updtdisp', .true. )

    !---set original configuration
    x0(1:3,1:n) = x(1:3,1:n)
end if

!-----transpose matrix of hi = inverse of h
CALL RCIPRL( h, hi, vtemp )

!---recover the original configuration
x(1:3,1:n) = x0(1:3,1:n)

do na = 1, nofatom(counter)
   i = idofatom(na,counter)
   dsp(1:3) = 0d0
   do ii = 1, 3
   do ix = 1, 3
      dsp(ii) = dsp(ii) + hi(ix,ii)*dispofatom(ix,na,counter)
   end do
   end do
!   x(1:3,i) = x0(1:3,i) + dsp(1:3)
   x(1:3,i) = x(1:3,i) + dsp(1:3)
end do


return
end




subroutine read_atomic_displacement( nfile, myid, nodes, iunit )
!-----------------------------------------------------------------------
! calculation with given atomic displacement
!-----------------------------------------------------------------------
use displacement_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit

!-----declare local variables
character(200) :: string
character(100) :: filename
integer :: i1, iend, i, n, nsum, nmax, nconfig
integer :: myiunit
real(8) :: the_mem
integer :: istat, status


myiunit = iunit

nconfig = 0
readdo: do

   read(myiunit,'(a)',iostat=istat) string
   !----- error trap
   if( istat /= 0 ) then
       if( myiunit == iunit ) then
           call fstop( nfile, myid, nodes,  &
&                               'error-0001 in read_atomic_displacement' )
       else
           exit readdo
       end if
   end if

   iend = len_trim(string)
   i1 = 1
   do
      if (string(i1:i1) /= ' ') exit
      i1 = i1 + 1
   end do
   if( string(i1:i1) == '#' ) cycle readdo
   if( string(i1:i1+2) == 'end' ) exit readdo

   !---set displacement
   if( string(i1:i1+2) == 'dis' ) then
       call read_d1( string, displacement )
       cycle readdo
   end if

   !---set displacement
   if( string(i1:i1+2) == 'fil' ) then
       call read_c1( string, filename )
       !-----open file
       call allocate_unit_number( myiunit )
       open(myiunit,iostat=istat,file=filename,status='old',form='formatted')
       if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                   'error-0002 in read_atomic_displacement' )
       if( myid == 0 ) then
          write(nfile(1),*) 'open file: ', trim(filename)
          write(nfile(2),*) 'open file: ', trim(filename)
       endif

       read(myiunit,'(a)',iostat=istat) string
       if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                   'error-0003 in read_atomic_displacement' )
!       cycle readdo
   end if

   !---set the number of displacements
   if( string(i1:i1+2) == 'Bas' ) then
       nofconfig = 0
       nsum = 0
       nmax = 0
       readdo2: do

          read(myiunit,'(a)',iostat=istat) string
          !----- error trap
          if( istat /= 0 ) then
              if( myiunit == iunit ) then
                  call fstop( nfile, myid, nodes,  &
&                             'error-0004 in read_atomic_displacement' )
              else
                  ! Let's suppose that there is no blank line in the last of the file
!                  backspace(myiunit)
                  exit readdo2
              end if
          end if

          iend = len_trim(string)
          i1 = 1
          do
             if (string(i1:i1) /= ' ') exit
             i1 = i1 + 1
          end do
          if( string(i1:i1) == '#' ) cycle readdo2
          if( string(i1:i1+2) == 'end' ) exit readdo2

          nofconfig = nofconfig + 1
          call read_i1( string, n )
          nsum = nsum + n
          nmax = max( nmax, n )
          do i = 1, n
             read(myiunit,'(a)',iostat=istat) string
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                         'error-0005 in read_atomic_displacement' )
          end do

       end do readdo2

       do i = 0, nofconfig + nsum
          backspace(myiunit)
       end do
       !---allocate memory
       allocate( nofatom(nofconfig), idofatom(nmax,nofconfig), disp(3,nmax,nofconfig),  &
& dispofatom(3,nmax,nofconfig), stat=status )
       the_mem = 4d0*( size(nofatom) + size(idofatom) + size(disp) )  &
&              + 8d0*( size(dispofatom) )

       !------error trap
       call check_alloc_accum( nfile, myid, nodes, status, the_mem, 'read_atomic_displacement', .true. )

       cycle readdo

   end if

   !---set displacement
   call read_i1( string, n )
   nconfig = nconfig + 1
   nofatom(nconfig) = n
   do i = 1, n
      read(myiunit,'(a)',iostat=istat) string
      !----- error trap
      if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                  'error-0004 in read_atomic_displacement' )
      read(string,*) idofatom(i,nconfig), disp(1:3,i,nconfig)
   end do

end do readdo
if( nconfig /= nofconfig ) call fstop( nfile, myid, nodes,  &
&                          'nconfig /= nofconfig in read_atomic_displacement' )

if( myiunit /= iunit ) then
    close(myiunit)
    call deallocate_unit_number( myiunit )
end if

!---set displacement in a.u.
do i = 1, nofconfig
   do n = 1, nofatom(i)
      dispofatom(1:3,n,i) = dble(disp(1:3,n,i)) * displacement
   end do
end do


return
end subroutine




subroutine read_d1( string, d1 )
!-----------------------------------------------------------------------
! read a real(8) variable from string
!-----------------------------------------------------------------------
implicit none
character(*) :: string
real(8) :: d1

!-----declare local variables
integer :: i1, i2, iend


iend = len_trim(string)

i1 = 1
do
   if (string(i1:i1) == '=' .or. string(i1:i1) == ':') exit
   i1 = i1 + 1
end do
i1 = i1 + 1

do
   if (string(i1:i1) /= ' ') exit
   i1 = i1 + 1
end do

i2 = i1 + 1
do
   if (string(i2:i2) == ' ') exit
   if( i2 >= iend ) exit
   i2 = i2 + 1
end do

read(string(i1:i2), *) d1


return
end subroutine




subroutine read_i1( string, d1 )
!-----------------------------------------------------------------------
! read an integer variable from string
!-----------------------------------------------------------------------
implicit none
character(*) :: string
integer :: d1

!-----declare local variables
integer :: i1, i2, iend


iend = len_trim(string)

i1 = 1
do
   if (string(i1:i1) == '=' .or. string(i1:i1) == ':') exit
   i1 = i1 + 1
end do
i1 = i1 + 1

do
   if (string(i1:i1) /= ' ') exit
   i1 = i1 + 1
end do

i2 = i1 + 1
do
   if (string(i2:i2) == ' ') exit
   if( i2 >= iend ) exit
   i2 = i2 + 1
end do

read(string(i1:i2), *) d1


return
end subroutine




subroutine read_c1( string, d1 )
!-----------------------------------------------------------------------
! read characters from string
!-----------------------------------------------------------------------
implicit none
character(*) :: string
character(*) :: d1

!-----declare local variables
integer :: i1, i2, iend


iend = len_trim(string)

i1 = 1
do
   if (string(i1:i1) == '=' .or. string(i1:i1) == ':') exit
   i1 = i1 + 1
end do
i1 = i1 + 1

do
   if (string(i1:i1) == '''' .or. string(i1:i1) == '"') exit
   i1 = i1 + 1
end do
i1 = i1 + 1

i2 = i1 + 1
do
   if (string(i2:i2) == '''' .or. string(i2:i2) == '"') exit
   if( i2 >= iend ) exit
   i2 = i2 + 1
end do
i2 = i2 - 1

d1 = ''
d1(1:i2-i1+1) = string(i1:i2)


return
end subroutine
