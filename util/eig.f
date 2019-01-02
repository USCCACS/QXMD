c---------------------------------------------------------------
c     plot eigenvalues
c---------------------------------------------------------------
      implicit real*8(a-h, o-z)
      dimension egvl(1000)
      dimension occ(1000)
      character*80 ofile, odirctry(20), native, argv
      character dummy
      character num(0:9)
      data num / '0','1','2','3','4','5','6','7','8','9' /


c------------------------------------------------
c--- No. of files to open
      nopen = 1
      odirctry(1) = '.'
c------------------------------------------------


      do i=1, command_argument_count()
         call get_command_argument(i,argv)
         select case(adjustl(argv))
          case ("-h", "--help")
                print'(a)', "./PDBcell -d ${output_directory}"
                stop
          case ("-n", "-nfile")
           call get_command_argument(i+1,argv); read(argv,*)noffil
            case ("-d", "-data")
             call get_command_argument(i+1,argv)
             odirctry(1)=adjustl(argv)
          case default
         end select
      enddo


      nini  = 0
      nskip = 1



      open( 10, file='EIG.dat' )



      ii1 = 0
      count = 0.0
      pcount = 0.0
      ncstp1 = 0
      ifop = 1
c------------------------------------------------
    2 continue
      nfile = ifop

      native = '/qm_eig.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 1, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      native = '/qm_fer.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 2, err=997, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      iffer = 1
      go to 887
  998 continue
      write(*,*) ofile,'does not exist.'
      stop
  997 continue
      iffer = 0
      write(*,*) ofile,'does not exist.'
  887 continue
      read(1,'(a1)') dummy
      if( iffer.eq.1 ) read(2,'(a1)') dummy


    1 continue
      read(1,*,end=999) ncstp, idummy, nband
      do ib = 1, nband
         read(1,*) ii, egvl(ib), occ(ib)
      enddo
      if( iffer.eq.1 ) read(2,*,end=999) ncstp, idummy, fermie 


      if( mod(ncstp,nskip).eq.0 .and. ncstp.ge.nini ) then

          write(*,*) ncstp

          if( iffer.ne.1 ) then
              do ib = 1, nband - 1
                 if( occ(ib).ge.1.0 .and. occ(ib+1).lt.1.0 ) then
                     fermie = ( egvl(ib) + egvl(ib+1) )*0.5
                 end if
              enddo
          end if

          count = count + 1.0
          do ib = 1, nband
             write(10,'(i6,e12.4)') ncstp, egvl(ib)-fermie
          enddo

      endif
      go to 1
  999 continue
      close(1)
      ifop = ifop + 1
      if( ifop.le.nopen ) then
          go to 2
      endif
c-------------------------------------------------------
      endfile(10)
      close(10)

      end




      subroutine getfname( odirctry, ofile, native )
      character*80 ofile, odirctry, native

      do i = 1, 80
         if( odirctry(i:i).ne.' ' ) then
             io1 = i
             go to 1
         endif
      enddo
    1 continue
      do i = 80, 1, -1
         if( odirctry(i:i).ne.' ' ) then
             io2 = i
             go to 2
         endif
      enddo
    2 continue
      do i = 1, 80
         if( native(i:i).ne.' ' ) then
             in1 = i
             go to 3
         endif
      enddo
    3 continue
      do i = 80, 1, -1
         if( native(i:i).ne.' ' ) then
             in2 = i
             go to 4
         endif
      enddo
    4 continue
      ofile = ( odirctry(io1:io2)//native(in1:in2) )
c      write(*,*) '1:', ofile
c      write(*,*) '2:', odirctry(io1:io2)
c      write(*,*) '3:', native(in1:in2)

      return
      end

