



subroutine get_nknod( nknod_ )
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: nknod_

nknod_  = 1 !nknod

return
end subroutine




subroutine get_nknodx( nknodx_ )
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: nknodx_

nknodx_  = 1 !nknodx

return
end subroutine




subroutine get_lgamma( lgamma_ )
!-----------------------------------------------------------------------
!use kp_variables
implicit none
logical :: lgamma_

lgamma_ = .true. !lgamma

return
end subroutine




subroutine get_pkscale( pkscale_ )
!-----------------------------------------------------------------------
!use kp_variables
implicit none
real*8  :: pkscale_

!if( lgamma ) then
    pkscale_ = 1.d0
!  else
!!    pkscale_ = 2.d0 * pkscale
!    if( .not.lnoncollinear ) then
!        pkscale_ = 2.d0 * pkscale
!    else
!        pkscale_ = 1.d0 * pkscale
!    end if
!end if

return
end subroutine




subroutine get_kpsi( lgamma_, nkpnt_, pkscale_ )
!-----------------------------------------------------------------------
!use kp_variables
implicit none
logical :: lgamma_
integer :: nkpnt_
real*8  :: pkscale_


lgamma_ = .true. !lgamma
nkpnt_  = 1 !nkpnt

!if( lgamma ) then
    pkscale_ = 1.d0
!  else
!!    pkscale_ = 2.d0 * pkscale
!    if( .not.lnoncollinear ) then
!        pkscale_ = 2.d0 * pkscale
!    else
!        pkscale_ = 1.d0 * pkscale
!    end if
!end if

return
end subroutine




subroutine get_lgammak( lgammak_ )
!-----------------------------------------------------------------------
!use kp_variables
!use kvector
implicit none
logical :: lgammak_

lgammak_ = .true.

return
end subroutine




subroutine unify_sumn( a, n, buff )
!-----------------------------------------------------------------------
!     unify by taking summation
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: n
real*8,  dimension(n) :: a, buff

return
end subroutine




subroutine unify_sum1( a )
!-----------------------------------------------------------------------
!     unify by taking summation
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: n
real*8  :: a

return
end subroutine




subroutine unify_max1( a )
!-----------------------------------------------------------------------
!     unify by taking max
!-----------------------------------------------------------------------
!use kp_variables
implicit none
real*8  :: a

return
end subroutine




subroutine unify_imax1( a )
!-----------------------------------------------------------------------
!     unify by taking max
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer  :: a

return
end subroutine




subroutine unify_by_dbcast( a, n )
!-----------------------------------------------------------------------
!     unify by broadcasting
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: n
real*8  :: a(n)

return
end subroutine




subroutine unify_by_lbcast( a, n )
!-----------------------------------------------------------------------
!     unify by broadcasting
!-----------------------------------------------------------------------
!use kp_variables
implicit none
integer :: n
logical :: a(n)

return
end subroutine




subroutine unify_by_dbcast1( a )
!-----------------------------------------------------------------------
!     unify by broadcasting
!-----------------------------------------------------------------------
!use kp_variables
implicit none
real*8  :: a

return
end subroutine




subroutine get_kvec( kvec_ )
!-----------------------------------------------------------------------
!     set k vector
!-----------------------------------------------------------------------
!use kvector
implicit none
integer :: kvec_

kvec_ = 1

return
end
