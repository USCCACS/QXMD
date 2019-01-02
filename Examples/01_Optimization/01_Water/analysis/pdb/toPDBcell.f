      program main
! --- Make PDB file. This program can be used for NPT-MD.
! --- Information of lattice constants will be included in PDB file.
! --- You can set "element" and "name" differently. (selection in VMD) 
! --- This program was made based on "toreal.f" by Masaaki MISAWA.
      parameter( namax=1000 )
      parameter( mxtype=6 )
      implicit real*8 (a-h,o-z)
      real*8, dimension(3,namax) :: x, r
      integer, dimension(namax) :: is, nhk
      character*80 ofile,odirctry(50), native
      character*1 dummy
      dimension number(mxtype)
      CHARACTER*2   ANAME
      COMMON/NATOM/ ANAME(103)
      real:: xxx,yyy,zzz,aaa,bbb,ccc
      character*4:: codnam(10)
      character*200::argv
      logical:: usecodnam
c----------------------------------------------------------------------
      noffil = 1
      odirctry(1) = '../../data'

c --- atomic number ---
      !number(1) = 8
      !number(2) = 1
      !number(3) = 1
      !number(4) = 1

c --- atomic code-name (up to 4 characters) ---
      usecodnam = .false. !--- false: same as element symbol.
      codnam(1) = 'Cu'
      codnam(2) = 'As'
      codnam(3) = 'O '
      codnam(4) = 'H '  !--- code-name is selectable as "name" in VMD


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


      nini = 0
      nskip = 1
      nstop = 1440000

      open(10,file='config.pdb',status='unknown')

c----------------------------------------------------------------------
      ncstp = -1
      nfile = 0
  100 continue
      nfile = nfile + 1
      
      native = ' /qm_ion.d'
      call getfname( odirctry(nfile), ofile, native )
      write(*,*) 'open : ', ofile
      open( 1, err=998, file=ofile, status='old' )

      native = '/qm_box.d'
      call getfname( odirctry(nfile), ofile, native)
      write(*,*) 'open : ', ofile
      open ( 2, err=998,file=ofile, status='old' )

      native = '/md_spc.d'
      call getfname(odirctry(nfile), ofile, native)
      write(*,*) 'open : ', ofile
      open ( 3, err=998,file=ofile, status='old' )
 
 2001 read(3,'(a1)') dummy
      if( dummy=='#' ) go to 2001
      backspace 3
      read(3,*,end=9999)
     &               ntype, (number(it), it = 1, ntype)

 1001 read(1,'(a1)') dummy
      if( dummy=='#' ) go to 1001
      backspace 1
    1 read(1,'(100i7)',end=9999) nstep,
     &               ntype, (nhk(it), it = 1, ntype)
      ntot = 0
      do it = 1, ntype
         ntot = ntot + nhk(it)
      end do
      read(1,*) atconfx
      read(1,'(9f8.5)')
     &         ( ( x(i,ia), i = 1, 3 ), ia = 1, ntot )

      ii = 0
      do it = 1, ntype
        do i = 1, nhk(it)
           ii = ii + 1
           is(ii) = it
           do ix = 1, 3
              x(ix,ii) = x(ix,ii)*atconfx
           enddo
        enddo
      enddo

c---shift      
!      do i = 1, ntot
!        x(2,i)=x(2,i)+0.3       ! shift in y vector of 0.3
!        if(x(2,i)>=1.d0)then
!            x(2,i)=x(2,i)-1.d0
!        end if
!      end do
c---end shift

 1000 read(2,'(a1)',end=996) dummy
      if( dummy=='#' ) go to 1000
    3  backspace 2
c-----------------------------------------------------------------------
    2 if(nstep > ncstp)then
        read(2,*,end=997) ncstp,
     &           dltca, dltcb, dltcc, angalf, angbet, anggam
      end if
      if(nstep == ncstp)then
        bohr = 5.29177d-1
        pi = acos(-1.d0)

        xxx=dltca*bohr
        yyy=dltcb*bohr
        zzz=dltcc*bohr
        aaa=angalf
        bbb=angbet
        ccc=anggam

        angalf = angalf * pi / 180.d0 
        angbet = angbet * pi / 180.d0 
        anggam = anggam * pi / 180.d0
        dltca2 = dltca 
        dltcb2 = dltcb
        dltcc2 = dltcc
        angalf2 = angalf
        angbet2 = angbet
        anggam2 = anggam

!      a = dltca + dltcb*cos(anggam) + bltcc * cos(angbet) 
!      b =         dltcb*sin(anggam) + bltcc * sin(angbet) * sin(angalf)
!      c =                             bltcc * sin(angbet) * cos(angalf)
        r3x = dltcc * cos(angbet)
        r3y = (dltcc * cos(angalf) - dltcc * cos(angbet)
     &                  * cos(anggam)) /sin(anggam)
!        r3y = dltcc * cos(angalf) * sin(anggam)
        if(dltcc**2-r3x**2-r3y**2 <0) write(*,*)r3x,r3y,dltcc
        r3z = sqrt(dltcc**2 - r3x**2 - r3y**2)
      end if

      ii = 0
      do it = 1, ntype
        do i = 1, nhk(it)
           ii = ii + 1
           is(ii) = it
           r(1,ii) = (dltca2*x(1,ii) + dltcb2*cos(anggam2)*x(2,ii) 
     &                              + r3x*x(3,ii)) * bohr
           r(2,ii) = (                dltcb2*sin(anggam2)*x(2,ii)
     &                              + r3y*x(3,ii)) * bohr
           r(3,ii) = (                r3z*x(3,ii)) * bohr
        enddo
        enddo
      if( mod(nstep,nskip).eq.0 )then
      if( nstep >= nini)then
      write(*,*) nstep
      write(10,'(a17,i5)') 'COMPND      STEP:',nstep
!      write(10,'(a37)') 'AUTHOR    GENERATED BY FUYUKI SHIMOJO'
      write(10,'(a6,3f9.3,3f7.2,12x,i4)')
     & 'CRYST1',xxx,yyy,zzz,aaa,bbb,ccc,ntot
      do ii = 1, ntot
c--- coordinates
         sx = r(1,ii)
         sy = r(2,ii)
         sz = r(3,ii)
       if(usecodnam .eqv. .false.)then
         write(10,'(a6,i5,1x,a2,3x,a3,i6,4x,3f8.3,2f6.2,10x,a2,2x)')
     & 'HETATM',ii,aname(number(is(ii))),'UNK',1, sx, sy, sz, 1.0, 0.0,
     & aname(number(is(ii)))
       end if
       if(usecodnam .eqv. .true.)then
         write(10,'(a6,i5,1x,a4,x,a3,i6,4x,3f8.3,2f6.2,10x,a2,2x)')
     & 'HETATM',ii,codnam(is(ii)),'UNK',1, sx, sy, sz, 1.0, 0.0,
     & aname(number(is(ii)))
       end if

      end do
!      write(10,'(a6,4x,12i5)') 'MASTER',0,0,0,0,0,0,0,0,ntot,0,ntot,0
      write(10,'(a3)') 'END'
      end if
      end if
 
      if( nstep.ge.nstop) stop
      go to 1

 9999 continue
      close(1)
      close(2)
      if( nfile < noffil ) go to 100
      stop
  996 continue
      backspace 2
      go to 3
  997 continue
      write(*,*) 'do not much : '
      stop
  998 continue
      write(*,*) 'cannot find : ', ofile
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


      BLOCK DATA SETATM

      CHARACTER*2   ANAME
      COMMON/NATOM/  ANAME(103)


      DATA  ANAME/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     & 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',
     & 'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     & 'Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     & 'Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce',
     & 'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     & 'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu',
     & 'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/

      END

