      parameter(  mxatom=700,    mxtype=5,   nkmax=1   )
      parameter(  ncmax = mxtype*mxtype   )
      parameter(  ndatgr = ncmax*11 )
      implicit real*8(a-h, o-z)
      character*80 ofile, odirctry(60), native
      dimension natom(mxtype)
      dimension r(3,mxatom),v(3,mxatom),r0(3,mxatom)
      dimension rba(3,3)!,rba0(3,3)
      !dimension pairt(0:10,mxtype)
      !dimension pair(0:10,ncmax)
      !dimension irecngb(mxatom)
      !dimension distc2(ncmax)
      !dimension neighbtmp1(20),neighbtmp2(mxatom,ncmax)
      !dimension distce(mxtype,mxtype)
      !dimension dis(mxatom,mxatom), pre_dis(mxatom,mxatom)
      character num(0:9), dummy
      !real SDR1,SDR2,SDR3
      data num / '0','1','2','3','4','5','6','7','8','9' /

      !real*8    :: ma0(mxtype)
      !real*8    :: m0
      !real*8    :: MV1,MV2,MV3
      !integer   :: ob1,ob2
      !character*3 :: cob1,cob2
      !character*20 :: filename 

      integer   :: nini
      integer   :: nfin

      integer   :: is(mxatom)

      integer :: atn(mxtype)!,it(mxatom)
      real*8  :: tmp(mxtype+1)!,it(mxatom)
      real*8  :: amass(mxtype),tmpmas
      real*8  :: VV, VX1, VY1, VZ1
      real*8  :: SDV1,SDV2,SDV3
      real*8  :: avogadro
      real*8  :: amu,m0,audang,e0,boltz,t0

c------------------------------------------------
      noffil = 1
      odirctry(1) = '../data01'
      odirctry(2) = '../data02'

      nini = 0
      nfin = 1000000

      avogadro = 6.022136736d23
      m0 = 9.10393d-28 !g
      audang =  0.529177249d0
      e0 = 4.35975d-18 !(J)
      boltz=1.38065d-23 !(J/K)
      t0 = 2.41888d-17 !(s)
      amu = 1.d0/avogadro/m0
      !write(*,*)"amu=",amu
      !stop

      open(20,file="tempsp.dat")

c------------------------------------------------
c    get supercell edges
      call celledg( noffil, ofile, odirctry, rba )
c------------------------------------------------
      count = 0.d0

      !iflag = 0
      nfile = 1
      !ifib = 0
   10 continue
      native = '/qm_ion.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 1, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile
      read(1,'(a1)') dummy

!      native = '/qm_box.d'
      native = '/qm_cel.d'
      call getfname( odirctry(nfile), ofile, native)
      write(*,*) 'open : ', ofile
      open ( 2, err=998,file=ofile, status='old' )
      read(2,'(a1)') dummy
      read(2,'(a1)') dummy

      native = '/md_vel.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 3, err=998, file=ofile, status='old', action='read' )
      write(*,*) 'open : ', ofile
      read(3,'(a1)') dummy

      native = '/md_spc.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 4, err=998, file=ofile, action ='read', status='old' )
      read(4,'(a1)') dummy
      read(4,*)ntype,(atn(k),k=1,ntype)
      close(4)

c----------------------------------------------------------
    1 read(1,'(10i7)',end=999) ncstp, ntype, (natom(it), it = 1, ntype)
      ntot = 0
      do it = 1, ntype
         ntot = ntot + natom(it)
      enddo
      read(1,*) fact
      read(1,'(9f8.5)') ( ( r(i,ia), i = 1, 3 ), ia = 1, ntot )

      do i = 1, ntot
      do j = 1, 3
         r(j,i) = r(j,i) * fact
      end do
      end do

      !definition of index of atom
      ntc = 0
      do it=1,ntype
         do j=1,natom(it)
            ntc = ntc + 1
            is(ntc) = it
         end do
      end do

      do it = 1, ntype
         call getatmas( atn(it), tmpmas )
         amass(it) = tmpmas/avogadro/m0
         !amass(it) = tmpmas
         !write(*,*)amass(it)
      end do
      !amass = amass/avogadro/m0
      !write(*,*)amass
      !stop

      read(3,'(11i7)',end=999) ncstp3, ntot
      read(3,*) atconfx
      read(3,'(9f8.5)')
     &         ( ( v(i,ia), i = 1, 3 ), ia = 1, ntot )
      do i = 1, ntot
      do j = 1, 3
         v(j,i) = v(j,i) * atconfx
      end do
      end do


!      if(ifib == 1)then
!         go to 7777
!      end if
    
!    2 read(2,*,end=777) ncstp2,
!!     &           dltca, dltcb, dltcc, angalf, angbet, anggam
!     &           rba(1:3,1:3)
!      rba = rba*0.529177210d0
!!      call cellget( dltca, dltcb, dltcc, angalf, angbet,
!!     &              anggam, rba)
!  777 continue
      
! 7777 if( ncstp  <  ncstp2)then
!          ifib = 1
!          go to 1
!      else if(ncstp > ncstp2)then
!          go to 2 
!      else
!          ifib = 0
!      end if

      if( ncstp.lt.nini ) go to 1
c      if( ncstp.ge.nend ) go to 1
      if( mod(ncstp,10).eq.0 ) write(*,*) ncstp
c-------------------------------------------------

!      if( iflag.eq.0 ) then
!          iflag = 1
!          ncm = 0
!c         the same species
!          do it = 1, ntype
!             ncm = ncm + 1
!c             write(*,*) 'distance in A for pair of ', it, it
!c             read(*,*) distce
!             distc2(ncm) = distce(it,it)*distce(it,it)
!          enddo
!c         different species
!          do it1 = 1, ntype
!          do it2 = 1, ntype
!          if( it1.ne.it2 ) then
!             ncm = ncm + 1
!c             write(*,*) 'distance in A for pair of ', it1, it2
!c             read(*,*) distce
!             distc2(ncm) = distce(it1,it2)*distce(it1,it2)
!          end if
!          end do
!          end do
!      end if

c   calculation of velocity distribution
      if(nini <= ncstp .and. ncstp <= nfin)then
      !dmax = 5 
      !ndis = 100
      !ddis = dmax/dble(ndis)
      !vdis = 0.0d0
      !velm = 0.d0
      tmp = 0.d0
      do i=1,ntot
         SDV1 = v(1,i)
         SDV2 = v(2,i)
         SDV3 = v(3,i)
         VX1  = rba(1,1)*SDV1+rba(1,2)*SDV2+rba(1,3)*SDV3
         VY1  = rba(2,1)*SDV1+rba(2,2)*SDV2+rba(2,3)*SDV3
         VZ1  = rba(3,1)*SDV1+rba(3,2)*SDV2+rba(3,3)*SDV3
         !VX1  = VX1*1d-10/t0 !(m/s)!/audang
         !VY1  = VY1*1d-10/t0 !(m/s)!/audang
         !VZ1  = VZ1*1d-10/t0 !(m/s)!/audang
         !VV   = dsqrt(VX1**2+VY1**2+VZ1**2)*1d-3 !(km/s)
         VV   = VX1**2+VY1**2+VZ1**2

         tmp(is(i)) = tmp(is(i)) + amass(is(i))*VV

!               write(1818,'(I6,I4,4(1x,F8.3))')
!     &         ncstp,i,VV,VX1*1d-3,VY1*1d-3,VZ1*1d-3

!         do j=1,ndis
!            if(ddis*dble(j-1) <= VV .and.
!     &         VV < ddis*dble(j))then
!               vdis(is(i),0) = vdis(is(i),0) + 1.d0
!               vdis(is(i),j) = vdis(is(i),j) + 1.d0

!               write(1818,'(I6,I4,4(1x,F8.3))')
!     &         ncstp,i,VV,VX1*1d-3,VY1*1d-3,VZ1*1d-3

!            if(ddis*dble(0) <= VV .and.
!     &         VV < ddis*dble(1))then    
!               write(1818,'(I6,I4,4(1x,F8.3))')
!     &         ncstp,i,VV,VX1*1d-3,VY1*1d-3,VZ1*1d-3
!            end if

!               exit
!            end if
!         end do
         !velm = max(VV,velm)
      end do 

      !tmp(ntype+1) = sum(tmp(1:ntype))
      do it=1,ntype
         tmp(ntype+1) = tmp(ntype+1) + tmp(it)
         tmp(it) = tmp(it)/dble(natom(it))*e0/boltz/3.d0
         !write(*,*)ncstp,natom(it)
      end do 
      tmp(ntype+1) = tmp(ntype+1)/dble(ntot)*e0/boltz/3.d0
      write(20,'(I7,10(1x,F8.1))') ncstp,(tmp(it),it=1,ntype+1)

      !write(*,*)ncstp,vdis(1:ntype,0),ntot
      !do it=1,ntype
      !   sum = 0.d0
      !   do j=1,ndis
      !      sum = sum + ddis*vdis(it,j)
      !   end do
      !   do j=1,ndis
      !      vdis(it,j) = vdis(it,j)/sum
      !   end do
      !end do
      !!write(9796,*)ncstp,velm
      !do it=1,ntype
      !   do j=1,ndis
      !      write(9000+it,*)ddis*dble(j-1),vdis(it,j)
      !   end do
      !end do
      end if

      go to 1
 999  continue
      close(1)
      nfile = nfile + 1
      if( nfile.le.noffil ) go to 10

! -- close files--
        close(20)
        !close(50)
       
 998  continue
      end




      subroutine celledg( noffil, ofile, odirctry, rba )
c------------------------------------------------
c    get supercell edges
c------------------------------------------------
      implicit real*8(a-h, o-z)
      parameter( audang =  0.529177249d0 )
      character*80 ofile, odirctry(20), native
      character dummy
      dimension rba(3,3)

      count = 0.d0
!      ava = 0.d0
!      avb = 0.d0
!      avc = 0.d0
!      avbc = 0.d0
!      avca = 0.d0
!      avab = 0.d0
      rba = 0.d0
      nfile = 1
!   10 continue
!      native = '/qm_box.d'
      native = '/qm_cel.d'
      call getfname( odirctry(nfile), ofile, native )
      write(*,*) 'open : ', ofile
      open( 1, err=998, file=ofile, status='old' )
      read(1,'(a1)') dummy
      read(1,'(a1)') dummy
c----------------------------------------------------------
!    1 read(1,*,end=999) ncstp,
    1 read(1,*) ncstp,
!     &           dltca, dltcb, dltcc, angalf, angbet, anggam
     &           rba(1:3,1:3)
      !rba = rba*audang
      count = count + 1.d0
!      ava = ava + dltca
!      avb = avb + dltcb
!      avc = avc + dltcc
!      avbc = avbc + angalf
!      avca = avca + angbet
!      avab = avab + anggam

! 999  continue
      close(1)
!      nfile = nfile + 1
!      if( nfile.le.noffil ) go to 10
c----------------------------------------------------------
!      dltca = ava/count * audang
!      dltcb = avb/count * audang
!      dltcc = avc/count * audang
!      avbc = avbc/count
!      avca = avca/count
!      avab = avab/count
!      angcon = acos(-1.d0)/180.d0
!      angalf = avbc * angcon
!      angbet = avca * angcon
!      anggam = avab * angcon
!C                           --- unit vectors parallel to cell vectors ---
!                                   E1X = 1.0
!                                   E1Y = 0.0
!                                   E1Z = 0.0
!                                   E2X = DCOS(anggam)
!                                   E2Y = DSIN(anggam)
!                                   E2Z = 0.0
!                                   E3X = DCOS(angbet)
!                                   E3Y = DCOS(angalf) - E3X*E2X
!                                   E3Y = E3Y/E2Y
!                                   E3Z = 1.0 - E3X*E3X - E3Y*E3Y
!                                   E3Z = DSQRT(E3Z)

!                                   rba(1,1) = dltca*E1X
!                                   rba(2,1) = dltca*E1Y
!                                   rba(3,1) = dltca*E1Z
!                                   rba(1,2) = dltcb*E2X
!                                   rba(2,2) = dltcb*E2Y
!                                   rba(3,2) = dltcb*E2Z
!                                   rba(1,3) = dltcc*E3X
!                                   rba(2,3) = dltcc*E3Y
!                                   rba(3,3) = dltcc*E3Z

      write(*,*) 'cell edges in [bohr]'
      write(*,'(3f10.5)') rba(1,1), rba(2,1), rba(3,1)
      write(*,'(3f10.5)') rba(1,2), rba(2,2), rba(3,2)
      write(*,'(3f10.5)') rba(1,3), rba(2,3), rba(3,3)

      return
 998  continue
      write(*,*) 'missing fort.47'
      stop
      end




      subroutine cellget( dltca, dltcb, dltcc, angalf,
     &                    angbet, anggam, rba )
c------------------------------------------------
c    get supercell edges
c------------------------------------------------
      implicit real*8(a-h, o-z)
      parameter( audang =  0.529177210d0 )
      character*80 ofile, odirctry(20), native
      dimension rba(3,3)
      character dummy

c      count = 0.d0
      count = 1.d0
      ava = 0.d0
      avb = 0.d0
      avc = 0.d0
      avbc = 0.d0
      avca = 0.d0
      avab = 0.d0
c      nfile = 1
c   10 continue
c      native = '/qm_box.d'
c      call getfname( odirctry(nfile), ofile, native )
c      write(*,*) 'open : ', ofile
c      open( 1, err=998, file=ofile, status='old' )
c      read(1,'(a1)') dummy
c      read(1,'(a1)') dummy
c----------------------------------------------------------
c    1 read(1,*,end=999) ncstp,
c     &           dltca, dltcb, dltcc, angalf, angbet, anggam
c      count = count + 1.d0
      ava = ava + dltca
      avb = avb + dltcb
      avc = avc + dltcc
      avbc = avbc + angalf
      avca = avca + angbet
      avab = avab + anggam

c 999  continue
c      close(1)
c      nfile = nfile + 1
c      if( nfile.le.noffil ) go to 10
c----------------------------------------------------------
      dltca = ava/count * audang
      dltcb = avb/count * audang
      dltcc = avc/count * audang
      avbc = avbc/count
      avca = avca/count
      avab = avab/count
      angcon = acos(-1.d0)/180.d0
      angalf = avbc * angcon
      angbet = avca * angcon
      anggam = avab * angcon
C                           --- unit vectors parallel to cell vectors
C                           ---
                                   E1X = 1.0
                                   E1Y = 0.0
                                   E1Z = 0.0
                                   E2X = DCOS(anggam)
                                   E2Y = DSIN(anggam)
                                   E2Z = 0.0
                                   E3X = DCOS(angbet)
                                   E3Y = DCOS(angalf) - E3X*E2X
                                   E3Y = E3Y/E2Y
                                   E3Z = 1.0 - E3X*E3X - E3Y*E3Y
                                   E3Z = DSQRT(E3Z)

                                   rba(1,1) = dltca*E1X
                                   rba(2,1) = dltca*E1Y
                                   rba(3,1) = dltca*E1Z
                                   rba(1,2) = dltcb*E2X
                                   rba(2,2) = dltcb*E2Y
                                   rba(3,2) = dltcb*E2Z
                                   rba(1,3) = dltcc*E3X
                                   rba(2,3) = dltcc*E3Y
                                   rba(3,3) = dltcc*E3Z

c      write(*,*) 'cell edges in [A]'
c      write(*,'(3f10.5)') rba(1,1), rba(2,1), rba(3,1)
c      write(*,'(3f10.5)') rba(1,2), rba(2,2), rba(3,2)
c      write(*,'(3f10.5)') rba(1,3), rba(2,3), rba(3,3)

      return
c 998  continue
c      write(*,*) 'missing qm_box.d'
c      stop
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

      return
      end



      subroutine getatmas(atn, tmpmas)
      integer :: atn
      real*8 :: tmpmas
      real*8, dimension(103) :: dmassn =                             
     & (/ 1.0079400d0,   4.0026020d0,   6.9410000d0,   9.0121820d0,  
     &  10.8110000d0,  12.0100000d0,  14.0067400d0,  15.9994000d0,      
     &  18.9984032d0,  20.1797000d0,  22.9897680d0,  24.3050000d0,      
     &  26.9815390d0,  28.0855000d0,  30.9737620d0,  32.0660000d0,      
     &  35.4527000d0,  39.9480000d0,  39.0983000d0,  40.0780000d0,      
     &  44.9559100d0,  47.8800000d0,  50.9415000d0,  51.9961000d0,      
     &  54.9380500d0,  55.8470000d0,  58.9332000d0,  58.6900000d0,      
     &  63.5460000d0,  65.3900000d0,  69.7230000d0,  72.6100000d0,      
     &  74.9215900d0,  78.9600000d0,  79.9040000d0,  83.8000000d0,      
     &  85.4678000d0,  87.6200000d0,  88.9058500d0,  91.2240000d0,      
     &  92.9063800d0,  95.9400000d0,  98.9100000d0, 101.0700000d0,      
     & 102.9055000d0, 106.4200000d0, 107.8682000d0, 112.4110000d0,      
     & 114.8200000d0, 118.7100000d0, 121.7500000d0, 127.6000000d0,      
     & 126.9044700d0, 131.2900000d0, 132.9054300d0, 137.3270000d0,      
     & 138.9055000d0, 140.1150000d0, 140.9076500d0, 144.2400000d0,      
     & 145.0000000d0, 150.0600000d0, 151.9650000d0, 157.2500000d0,      
     & 158.9253400d0, 162.5000000d0, 164.9303200d0, 167.2600000d0,      
     & 168.9342100d0, 173.0400000d0, 174.9670000d0, 178.4900000d0,      
     & 180.9479000d0, 183.8500000d0, 186.2070000d0, 190.2000000d0,      
     & 192.2200000d0, 195.0800000d0, 196.9665400d0, 200.5900000d0,      
     & 204.3833000d0, 207.2000000d0, 208.9803700d0, 210.0000000d0,      
     & 210.0000000d0, 222.0000000d0, 223.0000000d0, 226.0000000d0,      
     & 227.0000000d0, 232.0381000d0, 231.0358800d0, 238.0289000d0,      
     & 237.0500000d0, 244.0000000d0, 243.0000000d0, 247.0000000d0,      
     & 247.0000000d0, 251.0000000d0, 254.0000000d0, 257.0000000d0,      
     & 256.0000000d0, 254.0000000d0, 257.0000000d0 /)
     
      tmpmas = dmassn(atn)
     
      return
      end

