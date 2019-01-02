!-----------------------------------------------------------------------------!
!---  Pick up the configuration and velocity files                         ---!
!---                                                            2016/11/06 ---!
!---  Data are read from qm_ion.d, qm_box.d, qm_cel.d and md_vel.d.        ---!
!---  When you want to pick the configuration at last step, you set "npick"---!
!---  to "-1". And if lname is true, you can set the output file name to   ---!
!---  "basename".config and "basename".veloc.                              ---!
!---  This program is basen on "pick_config.f" and "pick_config_box.f".    ---!
!---                                                        By H. Kumazoe  ---!
!---  Update on 2016/12/04  by H. Kumazoe                                  ---!
!---    Impose PBC                                                         ---!
!-----------------------------------------------------------------------------!

program pick_config

   implicit none
   integer :: i, ia, j
   integer :: ios, istat, err
   integer :: iunit(1:4), ounit(0:3), ifile, nfile
   character(len=100) :: directory(5), dir, filebox
   integer :: npick, nstep(1:4), ntot, ifvel
   real(8) :: audang, dltca, dltcb, dltcc
   real(8) :: angalf, angbet, anggam, rba(1:3,1:3)
   integer, allocatable :: is(:)
   real(8), allocatable :: rscl(:,:), vscl(:,:), vel(:,:)
   real(8) :: rfactor, vfactor, vmax
   logical :: lname
   character(len=50) :: configname, velocname, cellname, basename
   character*200::argv
!-----------------------------------------------------------------------------!

   nfile = 1
!    directory(1) = ""
   directory(1) = "data"
   directory(2) = ""
   directory(3) = ""
   directory(4) = ""
   directory(5) = ""

   npick = -1
   npick = 200000
      do i=1, command_argument_count()
         call get_command_argument(i,argv)
         select case(adjustl(argv))
          case ("-h", "--help")
                print'(a)', "./PDBcell -d ${output_directory} -n $(TIMESTEP)"
                stop
          case ("-n", "-nstep")
           call get_command_argument(i+1,argv); read(argv,*)npick
          case ("-d", "-data")
             call get_command_argument(i+1,argv)
             directory(1)=adjustl(argv)
          case default
         end select
      enddo


   configname = "IN.CONFIG"
   velocname  = "IN.VELOC"
   cellname   = "Lattice_vector.dat"

   lname = .true.
   lname = .false.

   basename   = "test"
   if( lname ) then
       configname = trim( adjustl( basename))//".config"
       velocname  = trim( adjustl( basename))//".veloc"
   end if

!-----------------------------------------------------------------------------!


   audang = 0.529177249d0

   iunit(1) = 11 !  qm_ion.d
   iunit(2) = 12 !  qm_box.d
   iunit(3) = 13 !  qm_cel.d
   iunit(4) = 14 !  md_vel.d
   ounit(0) =  6 !  Standard output
   ounit(1) = 21 !  Velocity.dat
   ounit(2) = 22 !  Config.dat
   ounit(3) = 23 !  Lattice_vector.dat

   filebox = trim( adjustl( directory(1)))
   filebox = trim(filebox)//'/md_spc.d'
   open( unit=iunit(1), file=filebox, &
         iostat=ios, status="old", action="read" )
   if( ios/=0 ) stop "Error opening file md_spc.d"
   print*, "open  :", filebox

   read( unit=iunit(1), fmt="(a)", iostat=istat )  !  header
   read( unit=iunit(1), fmt="(i7)", iostat=istat ) !  ntype
   read( unit=iunit(1), fmt="(2i7)", iostat=istat ) nstep(1), ntot

   allocate( is(ntot), stat=err )
   if( err /= 0 ) print *, "is: Allocation request denied"
   allocate( rscl(1:3,1:ntot), stat=err )
   if( err /= 0 ) print *, "rscl: Allocation request denied"
   allocate( vscl(1:3,1:ntot), stat=err )
   if( err /= 0 ) print *, "vscl: Allocation request denied"
   allocate( vel(1:3,1:ntot), stat=err )
   if( err /= 0 ) print *, "vel: Allocation request denied"

   read( unit=iunit(1), fmt="(36i2)", iostat=istat ) &
         ( is(i), i = 1, ntot )

   close( unit=iunit(1), iostat=ios )

   filedo: do ifile = 1, nfile, 1
      dir = trim( adjustl( directory(ifile)))

      filebox = trim(dir)//'/qm_ion.d'
      open( unit=iunit(1), file=filebox, &
            iostat=ios, status="old", action="read" )
      if( ios/=0 ) stop "Error opening file qm_ion.d"
      print*, "open  :", filebox
      read( unit=iunit(1), fmt="(a)", iostat=istat ) !  header

      filebox = trim(dir)//'/qm_box.d'
      open( unit=iunit(2), file=filebox, &
            iostat=ios, status="old", action="read" )
      if( ios/=0 ) stop "Error opening file qm_box.d"
      print*, "open  :", filebox
      read( unit=iunit(2), fmt="(a)", iostat=istat ) !  header
      read( unit=iunit(2), fmt="(a)", iostat=istat ) !  header

      filebox = trim(dir)//'/qm_cel.d'
      open( unit=iunit(3), file=filebox, &
            iostat=ios, status="old", action="read" )
      if( ios/=0 ) stop "Error opening file qm_cel.d"
      print*, "open  :", filebox
      read( unit=iunit(3), fmt="(a)", iostat=istat ) !  header
      read( unit=iunit(3), fmt="(a)", iostat=istat ) !  header

      filebox = trim(dir)//'/md_vel.d'
      open( unit=iunit(4), file=filebox, &
            iostat=ios, status="old", action="read" )
      ifvel = ios
      if( ifvel==0 ) then
          print*, "open  :", filebox
          read( unit=iunit(4), fmt="(a)", iostat=istat ) !  header
      end if

      nstep(2) = -1
      readdo: do
         read( unit=iunit(1), fmt="(i7)", iostat=istat )  &
&              nstep(1)
         if( istat /= 0 ) npick = nstep(1)

         read( unit=iunit(1), fmt="(es14.7)", iostat=istat ) rfactor
         read( unit=iunit(1), fmt="(9f8.5)", iostat=istat )  &
&              (( rscl(i,ia), i = 1, 3), ia = 1, ntot )

         if( nstep(2)<nstep(1) ) then
             read( unit=iunit(2), fmt="(i7,3es14.7,3f11.6)", iostat=istat )  &
&                  nstep(2), dltca, dltcb, dltcc, angalf, angbet, anggam
             read( unit=iunit(3), fmt="(i7,9es15.7)", iostat=istat )  &
&                  nstep(3), rba(1:3,1:3)
             if( nstep(2)/=nstep(3) ) stop "nstep(2)/=nstep(3)"
         end if

         if( ifvel==0 ) then
             read( unit=iunit(4), fmt="(2i7)", iostat=istat ) nstep(4), ntot
             read( unit=iunit(4), fmt="(es14.7)", iostat=istat ) vfactor
             read( unit=iunit(4), fmt="(9f8.5)", iostat=istat )  &
&                  (( vscl(i,ia), i = 1, 3), ia = 1, ntot )
         end if

         if( nstep(1)==npick ) then
!              rscl(1:3,1:ntot) = rscl(1:3,1:ntot)*rfactor
             do ia = 1, ntot, 1
                do j = 1, 3, 1
                   rscl(j,ia) = rscl(j,ia)*rfactor
                   if( rscl(j,ia)>=1.d0 ) rscl(j,ia) = rscl(j,ia) - 1.d0
                   if( rscl(j,ia)< 0.d0 ) rscl(j,ia) = rscl(j,ia) + 1.d0
                end do
             end do
             rba(1:3,1:3) = rba(1:3,1:3)*audang
             dltca = dltca * audang
             dltcb = dltcb * audang
             dltcc = dltcc * audang

             if( ifvel==0 ) then
                 vscl(1:3,1:ntot) = vscl(1:3,1:ntot)*vfactor
                 vel(1:3,1:ntot) = matmul( rba(1:3,1:3), vscl(1:3,1:ntot) )
                 vel(1:3,1:ntot) = vel(1:3,1:ntot)/audang
                 vmax = maxval( abs( vel(1:3,1:ntot)))
                 vel(1:3,1:ntot) = vel(1:3,1:ntot)/vmax

                 write( unit=ounit(0), fmt="(a,i7,a)", iostat=ios )  &
&  "Pick up the configuration & velocity files at", nstep(1), " step"

                 open( unit=ounit(1), file=velocname, &
                       iostat=ios, status="replace", action="write" )

                 write( unit=ounit(1), fmt="(i7)", iostat=ios ) ntot
                 write( unit=ounit(1), fmt="(es15.7)", iostat=ios ) vmax
                 write( unit=ounit(1), fmt="(i2,3f10.6)", iostat=ios )  &
&                ( is(i), vel(1:3,i), i = 1, ntot )

                 close( unit=ounit(1), iostat=ios )
             else
                 write( unit=ounit(0), fmt="(a,i7,a)", iostat=ios )  &
&  "Pick up the configuration file at", nstep(1), " step"
             end if

             open( unit=ounit(2), file=configname, &
                   iostat=ios, status="replace", action="write" )
             write( unit=ounit(2), fmt="(i7)", iostat=ios ) ntot
             write( unit=ounit(2), fmt="(i2,3f9.6)", iostat=ios )  &
&            ( is(i), rscl(1:3,i), i = 1, ntot )
             close( unit=ounit(2), iostat=ios )

             open( unit=ounit(3), file=cellname, &
                   iostat=ios, status="replace", action="write" )

             write( unit=ounit(0), fmt="(a,i7,a)", iostat=ios )  &
&  "BOX at", nstep(2), " step [A]"
             write( unit=ounit(3), fmt="(a,i7,a)", iostat=ios )  &
&  "BOX at", nstep(2), " step [A]"
             write( unit=ounit(0), fmt="(f8.4,a,f8.4,a,f8.4,2x,a)", iostat=ios )  &
&  dltca,", ", dltcb,", ", dltcc, ":  angles between cell vec. in [deg.]"
             write( unit=ounit(3), fmt="(f8.4,a,f8.4,a,f8.4,2x,a)", iostat=ios )  &
&  dltca,", ", dltcb,", ", dltcc, ":  angles between cell vec. in [deg.]"
             write( unit=ounit(0), fmt="(f8.4,a,f8.4,a,f8.4,2x,a)", iostat=ios )  &
&  angalf,", ", angbet,", ", anggam, ":  lengths of cell vectors"
             write( unit=ounit(3), fmt="(f8.4,a,f8.4,a,f8.4,2x,a)", iostat=ios )  &
&  angalf,", ", angbet,", ", anggam, ":  lengths of cell vectors"
             write( unit=ounit(0), fmt="(a,i7,a)", iostat=ios )  &
&  "CELL at", nstep(2), " step [A]"
             write( unit=ounit(3), fmt="(a,i7,a)", iostat=ios )  &
&  "CELL at", nstep(2), " step [A]"
             write( unit=ounit(0), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,1), ": super cell vector L1"
             write( unit=ounit(3), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,1), ": super cell vector L1"
             write( unit=ounit(0), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,2), ": super cell vector L2"
             write( unit=ounit(3), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,2), ": super cell vector L2"
             write( unit=ounit(0), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,3), ": super cell vector L3"
             write( unit=ounit(3), fmt="(3f9.5,3x,a)", iostat=ios )  &
&  rba(1:3,3), ": super cell vector L3"

             close( unit=ounit(3), iostat=ios )
             exit readdo
         end if
      end do readdo
      close( unit=iunit(1), iostat=ios )
      close( unit=iunit(2), iostat=ios )
      close( unit=iunit(3), iostat=ios )
      close( unit=iunit(4), iostat=ios )
   end do filedo

   if( allocated(is) ) deallocate( is, stat=err )
   if( err /= 0 ) print *, "is: Deallocation request denied"
   if( allocated(rscl) ) deallocate( rscl, stat=err )
   if( err /= 0 ) print *, "rscl: Deallocation request denied"
   if( allocated(vscl) ) deallocate( vscl, stat=err )
   if( err /= 0 ) print *, "vscl: Deallocation request denied"
   if( allocated(vel) ) deallocate( vel, stat=err )
   if( err /= 0 ) print *, "vel: Deallocation request denied"

end program pick_config
