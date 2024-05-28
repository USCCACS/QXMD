!==================================================================================!
!start module
!==================================================================================!

module params 

!mxatom is total no of atom allowed for the program 
!ini_band - starting no of band 
!end_band - end no of band 
!suppose if you want to plot all the cube file for 80 to 89 put
! ini_band=80 and end_band=89
!nskip_step= how frequently have you printed wave function
! if you want to plot last point; put nskip to total no of run
! fixbox switch enable have volume change
!natomic(6) provide which atom. This values is arbitraty and you can choose any no of atom 
!root_dir= give the relative path of data file
!

integer, parameter:: mxatom = 500
!=====set up initial band to plot cube file 
integer, parameter:: ini_band= 4
integer, parameter:: end_band= 5
!=====set up final band to plot cube file 
integer, parameter:: nskip_step=10
logical, parameter:: fixbox=.true.
integer, parameter:: ncmax = 20000
integer, parameter:: nstop =300
integer, parameter:: mxtype=10
real(8), parameter:: btoA=0.529177249
real(8), parameter:: pi= acos(-1.0)
integer, parameter:: natomic(6)= (/8,1, 6,0, 0, 0/)

character(100):: root_dir,out_dir
integer:: natom(mxtype)
real(8):: egvl(2,3000), occ(2,3000)
real(8):: r(3,mxatom)
integer:: is(mxatom)
real(8) :: HH(3,3)
integer:: status, ntype, status1, status2, status3

integer:: nx, ny, nz, icount
real(8)::ival
real(8), allocatable:: val(:,:,:)
real(8):: mesh_factor
end module params

program create_cube_file
use params 
implicit none

integer:: nstep, ntot
integer:: i, j, kk, it, ia, ib, ii, idummy, ik
integer:: natm1, natm2, nband
real:: factor
real:: v1, v2, v3
integer:: ix, iy, iz, itot
character(100):: trj, file_eig
character(100):: band_num, step_num
character(200):: filename1, filename2, filename3, filename4
integer:: ntemp1, ntemp2

!====set here to path of data file!
root_dir= '../../data'
out_dir='./'
!=====end data file setup 

trj='/qm_ion.d'
filename1= trim(root_dir)//trim(trj)
open(1, file=filename1, status='old', iostat=status)
    if (status /=0) print *, "could not open file:" , filename1
    print *, "Open:" , filename1 


    ! read dummy lien 
    read(1,*)

    do kk=1, nstop
        read(1,*, iostat=status1) nstep, ntype, (natom(it), it=1, ntype)
        if(status1 /=0) stop
        !print *, nstep, ntype, (natom(it), it=1, ntype)
        read(1,*) factor
        !reset ntot
        ntot=0

        do it=1, ntype
            ntot= ntot+natom(it)
            !print *, ntot, natom(it)
        enddo 

        read(1,'(9f8.5)') ( ( r(i,ia), i = 1, 3 ), ia = 1, ntot )
        
        ii=0
        do it=1, ntype
            do i= 1, natom(it)
                ii=ii+1
                is(ii)= it   
                r(1,ii)= r(1, ii)* factor
                r(2,ii)= r(2, ii)* factor
                r(3,ii)= r(3, ii)* factor
            enddo 
        enddo
    !get box size 
        if (fixbox) then  
            if (kk==1) call celledg(nstep)
        else 
            call celledg(nstep)
        endif
        !print *, HH(1:3, 1:3)
    ! compute real Coordinate from scaled coordinate 
        do ia = 1, ntot
            v1 = r(1,ia)
            v2 = r(2,ia)
            v3 = r(3,ia)
            r(1,ia) = HH(1,1)*v1 + HH(1,2)*v2 + HH(1,3)*v3
            r(2,ia) = HH(2,1)*v1 + HH(2,2)*v2 + HH(2,3)*v3
            r(3,ia) = HH(3,1)*v1 + HH(3,2)*v2 + HH(3,3)*v3
        enddo
        if (mod(nstep, nskip_step) /=0) cycle  
            do ii= ini_band, end_band
               !if (mod(ii, nskip_band) /=0) cycle
                write(band_num,'(I2.2)') ii
                write(step_num,'(I6.6)') nstep
                filename2= trim(root_dir)//"/qm_eigv.d."//trim(band_num)//"."//trim(step_num)
                filename3= trim(out_dir)//"/state."//trim(band_num)//"."//trim(step_num)//".cube"
                open(10, file=filename2, status='old', action='read')
                open(11, file=filename3, action='write')
                print *, "Open: ", filename3
                read(10, *) !nx, ny, nz
                read(10, *) nx, ny, nz
                !print *, nx, ny, nz

                ! allocate data
                allocate(val(nx, ny, nz))

                read(10, *) mesh_factor
                read(10,*) icount
                read(10,*) ival
                i = 0
                do iz = 1, nz
                do iy = 1, ny
                do ix = 1, nx
                    i = i + 1
                    if( i > icount ) then
                    read(10,*) icount
                    read(10,*) ival
                    i = 1
                    end if
                    val(ix,iy,iz) = ival*mesh_factor
                !    print *, val(ix, iy, iz)
                end do
                end do
                end do
                close(10)

                write(11,'(a)') "CPMD CUBE FILE."
                write(11,'(a)') "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
                write(11,'(i4,3f12.6)') ntot, 0.d0, 0.d0, 0.d0
                write(11,'(i4,3f12.6)') nx,HH(1,1:3)/nx
                write(11,'(i4,3f12.6)') ny,HH(2,1:3)/ny
                write(11,'(i4,3f12.6)') nz,HH(3,1:3)/nz

                do ik=1,ntot
                write(11,'(i4,4f12.5)') natomic(is(ik)),0.d0,r(1:3,ik)
                enddo
                itot=0
                do ix=1,nx
                do iy=1,ny
                do iz=1,nz
                write(11,'(e14.5 $)') val(ix,iy,iz)*val(ix,iy,iz)
                itot=itot+1
                if(mod(itot,6)==0) write(11,*)
                enddo; enddo; enddo
                close(11)
                deallocate(val)

            enddo 
    enddo ! ending main loop

    ! end main loop
end program create_cube_file

!=================================================================================!
!find cell dimension from qm_box.d
!================================================================================!
subroutine celledg(nstep)
  use params
  implicit none
  integer:: nstep, i, ntemp
  integer:: MAXSTEP=100000
  character(100):: box
  character(200)::filename
  real(8):: la, lb, lc, lalpha, lbeta, lgamma
  real:: E1X, E1Y, E1Z, E2X, E2Y, E2Z, E3X, E3Y, E3Z

  !print* , nstep

  box='/qm_box.d'

  filename=trim(root_dir)//trim(box)
  !!print *, filename
  open(10, file=filename, status='old')
  read(10,*)
  read(10,*)
  do i= 0 , MAXSTEP
  read(10,*) ntemp, la, lb, lc, lalpha, lbeta, lgamma
     if (ntemp==nstep) exit
  enddo
  !print *, nstep, la, lb, lc, lalpha, lbeta, lgamma
  close(10)
  !! convert length to Angstrom and angle to radian
  la= la!*btoA
  lb= lb!*btoA
  lc= lc!*btoA
  lalpha= lalpha*pi/180
  lbeta = lbeta*pi/180
  lgamma = lgamma*pi/180
   !print *, la, lb, lc, lalpha, lbeta, lgamma

  !----unit vectors
  E1X = 1.0
  E1Y = 0.0
  E1Z = 0.0
  E2X = COS(lgamma)
  E2Y = SIN(lgamma)
  E2Z = 0.0
  E3X = COS(lbeta)
  E3Y = COS(lalpha) - E3X*E2X
  E3Y = E3Y/E2Y
  E3Z = 1.0 - E3X*E3X - E3Y*E3Y
  E3Z = SQRT(E3Z)
  !---compute matrix

HH(1,1) = la*E1X
HH(2,1) = la*E1Y
HH(3,1) = la*E1Z
HH(1,2) = lb*E2X
HH(2,2) = lb*E2Y
HH(3,2) = lb*E2Z
HH(1,3) = lc*E3X
HH(2,3) = lc*E3Y
HH(3,3) = lc*E3Z
!write(*,*) 'cell edges in [A]'
!write(*,'(3f10.5)') HH(1,1), HH(2,1), HH(3,1)
!write(*,'(3f10.5)') HH(1,2), HH(2,2), HH(3,2)
!write(*,'(3f10.5)') HH(1,3), HH(2,3), HH(3,3)
!print *, nstep, la, lb, lc, lalpha, lbeta, lgamma

end subroutine
!
