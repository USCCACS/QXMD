



module outfile
!-----------------------------------------------------------------------
! type declaration of variables to write data to files
!-----------------------------------------------------------------------
implicit none

logical :: loutfile(2) = .false.   ! for QM nodes
! loutfile(1) == .true. : write data to the standard output
! loutfile(2) == .true. : write data to the log file

save

end module




module allocated_memory
!-----------------------------------------------------------------------
! type declaration of variables for memory allocation
!-----------------------------------------------------------------------
implicit none

real*8  :: alloc_mem = 0.d0    ! counter for allocated memory

save

end module
