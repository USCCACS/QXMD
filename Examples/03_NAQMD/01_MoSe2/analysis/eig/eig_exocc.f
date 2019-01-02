c---------------------------------------------------------------
c     plot eigenvalues
c   extract the highest occupied energy
c---------------------------------------------------------------
      implicit real*8(a-h, o-z)
      dimension egvl(2,1000)
      dimension occ(2,1000)
      character*80 ofile, odirctry(20), native, argv
      character dummy
      character num(0:9)
      data num / '0','1','2','3','4','5','6','7','8','9' /


c------------------------------------------------
c--- No. of files to open
      nopen = 1
      odirctry(1) = '../../data'
      odirctry(2) = 'data2'

      nini  = 0
      nskip = 1

c------------------------------------------------
!      write(*,*) 'No. of extract occupied bands?'
!      read(*,*) nocc


      open( 10, file='EIG.dat' )
      open( 11, file='EIG_occ-two.dat' )
      open( 12, file='EIG_occ-one.dat' )

      do i=1, command_argument_count()
         call get_command_argument(i,argv)
         select case(adjustl(argv))
          case ("-h", "--help")
            print'(a)', "Calculate Eigen function and Occupancy"
                print'(a)', "./eig -d ${data_
     &directory}"
           print'(a)', "-d : Define path of data directory"
                stop
            case ("-d", "-data")
             call get_command_argument(i+1,argv)
             odirctry(1)=adjustl(argv)
          case default
         end select
      enddo



      ii1 = 0
      count = 0.0
      pcount = 0.0
      ncstp1 = 0
      ifop = 1
      evryd = 13.60568
c------------------------------------------------
    2 continue
      nfile = ifop

      native = '/qm_td_eig.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 1, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      native = '/qm_fer.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 2, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      go to 887
  998 continue
      write(*,*) ofile,'does not exist.'
      stop
  887 continue
      read(1,'(a1)') dummy
      read(2,'(a1)') dummy


    1 continue
      read(1,*,end=999) ncstp, idummy, nband
      do ib = 1, nband
         read(1,*) ii, egvl(1,ib), occ(1,ib) !, jj, egvl(2,ib), occ(2,ib)
         egvl(1,ib) = egvl(1,ib) * evryd
      enddo
      read(2,*,end=999) ncstp, idummy, fermie 
      fermie = fermie*evryd


      if( mod(ncstp,nskip).eq.0 .and. ncstp.ge.nini ) then

          write(*,*) ncstp

          count = count + 1.0
          do ib = 1, nband
             eneup = egvl(1,ib)
             write(10,'(i6,2e12.4)') ncstp, eneup-fermie
             if( 2.1 > occ(1,ib) .and. occ(1,ib) > 1.5 ) then
                 write(11,'(i6,2e12.4)') ncstp, eneup-fermie
             else if( 1.5 >= occ(1,ib) .and. occ(1,ib) > 0.5 ) then
                 write(12,'(i6,2e12.4)') ncstp, eneup-fermie
             end if
          end do

      endif
      go to 1
  999 continue
      close(1)
      ifop = ifop + 1
      if( ifop.le.nopen ) then
          go to 2
      endif
c-------------------------------------------------------
      close(10)
      close(11)
      close(12)

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

