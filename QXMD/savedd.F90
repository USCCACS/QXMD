



module savedata
!-----------------------------------------------------------------------
! type declaration of variables for data save
!-----------------------------------------------------------------------
implicit none

logical :: lsreal8
logical :: lsifreal
integer*2, allocatable, dimension(:)   :: saveigv
integer*2, allocatable, dimension(:,:) :: sav2slm

save

end module




subroutine set_lsreal8( nfile, myid, nodes, lsreal8_ )
!-----------------------------------------------------------------------
! set lsreal8 in savedata
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
logical :: lsreal8_

lsreal8 = lsreal8_

return
end




subroutine rdacon( nfile, iogpsz, ct0,  &
& fname1, lstart, prevr, nstepCG, nstepMD,  &
& ratm, frc, ntype, nhk1, nhk2, nhk1r, nhk2r, lclust, ifmd,  &
& natom, hcell )
!-----------------------------------------------------------------------
!    restore atomic configuration
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension   nfile(*)
integer :: iogpsz
character(*) :: fname1
dimension   prevr(3,natom,3)
dimension   ratm(3,natom)
dimension   frc(3,natom)
dimension   nhk1(ntype), nhk2(ntype)
dimension   nhk1r(ntype), nhk2r(ntype)
real*8, dimension(3,3) :: hcell
real*8, dimension(3,3) :: hcellr = 0.d0
logical     lstart, lclust, lclusr

!-----declare local variables
character(80) :: fname
integer :: myid, nodes, myid_qm_un, nodes_qm_un, digit
integer :: iunit
integer :: nfiles, ii, isnd, nsd, ibuf(2)


!-----set communicator
call get_worldqm( myid, nodes )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


ifmdr   = 0
do j = 1, 3
do i = 1, natom
do ix = 1, 3
   prevr(ix,i,j) = 0.d0
end do
end do
end do
if( lstart ) then
    ierror = 0
    if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

        fname = trim(fname1)

        if( myid == 0 ) then
            if(loutfile(1)) write(nfile(1),*) 'open file(in rdacon): ', trim(fname)
            if(loutfile(2)) write(nfile(2),*) 'open file(in rdacon): ', trim(fname)
        end if

        call allocate_unit_number( iunit )

        open(iunit,iostat=istat,file=fname,status='old', form='unformatted')

        if( istat == 0 ) read(iunit,iostat=istat) nversion
        !-----keep compatibility with old version
        if( istat == 0 ) then
            if( nversion > 0 ) then
                ntyper = nversion
              else
                read(iunit,iostat=istat) ntyper
            end if
        end if
        if( ntyper.ne.ntype ) then
            write(nfile(1),*) 'ntyper.ne.ntype', ntyper, ntype
            write(nfile(2),*) 'ntyper.ne.ntype', ntyper, ntype
            istat = 1
        end if
        if( istat == 0 ) read(iunit,iostat=istat)  &
&                       ( nhk1r(it), nhk2r(it), it = 1, ntype )
        do it = 1, ntype
           if( nhk1r(it).ne.nhk1(it) .or. nhk2r(it).ne.nhk2(it) ) then
               write(nfile(1),*) 'error in nhk1, nhk2'
               write(nfile(2),*) 'error in nhk1, nhk2'
               istat = 1
               exit
           end if
        end do
        if( istat == 0 ) read(iunit,iostat=istat) nstepCG, nstepMD
        if( istat == 0 ) read(iunit,iostat=istat) lclusr
        if(      lclusr .and. .not.lclust .or.  &
&           .not.lclusr .and. lclust            ) then
            write(nfile(1),*) 'error in logical lclusr'
            write(nfile(2),*) 'error in logical lclusr'
            istat = 1
        end if

        if( istat == 0 ) read(iunit,iostat=istat) &
&                      ((ratm(ix,i), ix = 1, 3), i = 1, nhk2(ntype))

        if( istat == 0 ) read(iunit,iostat=istat) ifmdr
        if( ifmdr.ne.ifmd ) then
            write(nfile(1),*) 'warning: ifmdr.ne.ifmd', ifmdr,ifmd
            write(nfile(2),*) 'warning: ifmdr.ne.ifmd', ifmdr,ifmd
        end if

        if( istat == 0 ) read(iunit,iostat=istat) &
&                      ((frc(ix,i), ix = 1, 3), i = 1, nhk2(ntype))
        if( istat == 0 ) read(iunit,iostat=istat) &
&                     (((prevr(ix,i,j), ix = 1, 3), i = 1, nhk2(ntype)), j = 1, 3)

        if( nversion <= 0 ) then
            if( istat == 0 ) read(iunit,iostat=istat) &
&                        ( ( hcellr(ix,j), ix = 1, 3 ), j = 1, 3 )
        end if
        if( istat == 0 ) then
            !-----keep compatibility with nversion == -1
            if( nversion == -1 ) then
                if( .not.lclust ) then
                    do j = 1, 3
                    do i = 1, nhk2(ntype)
                    do ix = 1, 3
                       prevr(ix,i,j) = prevr(ix,i,j)/hcellr(ix,ix)
                    end do
                    end do
                    end do
                end if
            end if
        else
            ierror = 1
        end if
        close(iunit)

        call deallocate_unit_number( iunit )

        if( ifmd.ge.1 .and. ifmd.ne.ifmdr ) then
            if( ifmd.eq.1 ) then
                nstepMD = 0
              else
                nstepCG = 0
            end if
        end if

        !-----for divide-and-conquer MD
        do ii=1, nfiles-1
           ibuf(1) = nstepMD
           ibuf(2) = nstepCG
           call cisend(myid_qm_un+ii,ibuf,2,myid_qm_un+ii,0)
           call cdsend(myid_qm_un+ii,hcellr,9,myid_qm_un+ii,0)
        end do

  else ioif

        !-----for divide-and-conquer MD
        isnd=(myid_qm_un/iogpsz)*iogpsz
        call cirecvs(myid_qm_un,ibuf,2,isnd,0)
        call cdrecvs(myid_qm_un,hcellr,9,isnd,0)
        nstepMD = ibuf(1)
        nstepCG = ibuf(2)

  end if ioif


    end if

    !-----set communicator
    call get_worldqm( myid, nodes )

    call gimax(ierror)
    if( ierror.gt.0 ) then
        call fstop( nfile, myid, nodes, 'error : in file ion000' )
    end if

    ntot  = nhk2(ntype)
    ntot3 = 3*nhk2(ntype)
    call dbcast(ratm,ntot3,0)

    call gimax(nstepCG)
    call gimax(nstepMD)
    if( myid.eq.0 ) then
        write(nfile(1),*) ' nstepCG, nstepMD =', nstepCG, nstepMD
        write(nfile(2),*) ' nstepCG, nstepMD =', nstepCG, nstepMD
    end if

    if( ifmd.ge.1 .and. max(nstepCG,nstepMD).gt.0 ) then
        call dbcast(frc,ntot3,0)
        call dbcast(prevr(1,1,1),ntot3,0)
        call dbcast(prevr(1,1,2),ntot3,0)
        call dbcast(prevr(1,1,3),ntot3,0)
    end if

    !-----check supercell size in the previous calculation
    call dbcast(hcellr,9,0)
    epsilon = 1.d-12
    if( abs(hcell(1,1)-hcellr(1,1)) > epsilon .or.  &
&       abs(hcell(2,2)-hcellr(2,2)) > epsilon .or.  &
&       abs(hcell(3,3)-hcellr(3,3)) > epsilon .or.  &
&       abs(hcell(1,2)-hcellr(1,2)) > epsilon .or.  &
&       abs(hcell(2,1)-hcellr(2,1)) > epsilon .or.  &
&       abs(hcell(2,3)-hcellr(2,3)) > epsilon .or.  &
&       abs(hcell(3,2)-hcellr(3,2)) > epsilon .or.  &
&       abs(hcell(3,1)-hcellr(3,1)) > epsilon .or.  &
&       abs(hcell(1,3)-hcellr(1,3)) > epsilon      ) then
        lstart = .false.
        if( myid.eq.0 ) then
         write(nfile(1),*) ' *** The supercell size is different',  &
&                          ' from the pervious one.'
         write(nfile(2),*) ' *** The supercell size is different',  &
&                          ' from the pervious one.'
        end if
    end if

    ct = timecnt()
    if( myid.eq.0 ) then
        if(loutfile(1)) write(nfile(1),*) ' set atomic config. : cpu-time :', ct - ct0
        if(loutfile(2)) write(nfile(2),*) ' set atomic config. : cpu-time :', ct - ct0
    end if
    ct0 = ct
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




subroutine rdwfns( nfile, iogpsz, ct0,  &
& fname1, lstart, cgjr, eig, egvlud, occ, wegud, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, lspin, gdcrsv, nspnod, lcgjsv, lwfrand, ncscale, wvratio )
!-----------------------------------------------------------------------
!     set initial wave functions decomposed by bands : cgjr
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
integer :: iogpsz
character(*) :: fname1
dimension cgjr(nplwex*ncscale,*)
dimension gdcr(nplwex*ncscale,*)
dimension nbndsp(*), nbncnt(*)
dimension eig(*), occ(*)
logical   lstart, lsifreal
dimension egvlud(*), wegud(*)
dimension gdcrsv(*)
logical   lspin, lcgjsv
logical :: lwfrand
real*8 :: seed = 25.d0
integer :: myid, nodes, nkd
logical :: lreadocc
integer :: ncscale
real*8  :: wvratio(*)

integer :: nplwoex2
real*8  :: dummy

!-----declare local variables
character(80) :: fname
integer :: iunit, digit
logical :: lsetnspin12
integer :: nfiles, ii, isnd, nsd, ibuf(7)
real(8), allocatable :: abuf(:,:)


!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then
nkdif: if( nkd == 0 ) then

if( myid.eq.0 ) then
    if(loutfile(1)) write(nfile(1),*) ' initial wave functions'
    if(loutfile(2)) write(nfile(2),*) ' initial wave functions'
end if

lsetnspin12 = .false.
startif: if( lstart ) then
!=======================================================================

ierror = 0
if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    fname = trim(fname1)

    if( myid == 0 ) then
        if(loutfile(1)) write(nfile(1),*) 'open file(in rdwfns): ', trim(fname)
        if(loutfile(2)) write(nfile(2),*) 'open file(in rdwfns): ', trim(fname)
    end if

    call allocate_unit_number( iunit )

    open(iunit,iostat=istat,file=fname,status='old',form='unformatted')
    ncscaleo = 1
    if( istat == 0 ) read(iunit,iostat=istat) nbando
    lreadocc = nbando < 0
    if( istat == 0 .and. lreadocc ) then
        ncscaleo = abs(nbando)
        if( istat == 0 ) read(iunit,iostat=istat) nbando
    end if
    if( istat == 0 ) read(iunit,iostat=istat) nplwo
    if( istat == 0 ) read(iunit,iostat=istat) lsifreal
    if( istat == 0 ) read(iunit,iostat=istat) nnspin

    if( istat == 0 ) then
        if( nbando.ne.nband .and. myid.eq.0 ) then
            write(nfile(1),*) 'warning - different nband in ', trim(fname)
            write(nfile(2),*) 'warning - different nband in ', trim(fname)
        end if
        if( ncscaleo.ne.ncscale .and. myid.eq.0 ) then
            write(nfile(1),*) 'warning - different ncscale in ', trim(fname)
            write(nfile(2),*) 'warning - different ncscale in ', trim(fname)
            if( ncscaleo > ncscale ) then
                ierror = 1
            end if
        else
            if( nplwo.ne.nplw .and. myid.eq.0 ) then
                ierror = 1
                write(nfile(1),*) 'error - different nplw in ', trim(fname)
                write(nfile(2),*) 'error - different nplw in ', trim(fname)
            end if
        end if
        nbandrec = nbando
        nbando = min( nbando, nband )
        nplwo = min( nplwo, nplw )
        nplwoex = 2*( nplwo + 1 )
    else
        ierror = 1
    end if

    !-----for divide-and-conquer MD
    do ii=1, nfiles-1
       ibuf(1) = nbando
       ibuf(2) = nplwo
       if( lsifreal ) then
           ibuf(3) = 1
       else
           ibuf(3) = 0
       end if
       ibuf(4) = nnspin
       ibuf(5) = nbandrec
       ibuf(6) = ncscaleo
       ibuf(7) = ierror
       call cisend(myid_qm_un+ii,ibuf,7,myid_qm_un+ii,0)
    end do

  else ioif

    !-----for divide-and-conquer MD
    isnd=(myid_qm_un/iogpsz)*iogpsz
    call cirecvs(myid_qm_un,ibuf,7,isnd,0)
    nbando   = ibuf(1)
    nplwo    = ibuf(2)
    lsifreal = ibuf(3) == 1
    nnspin   = ibuf(4)
    nbandrec = ibuf(5)
    ncscaleo = ibuf(6)
    ierror   = ibuf(7)
    nplwoex = 2*( nplwo + 1 )

  end if ioif

end if

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!--- error trap
call gimax(ierror)
if( ierror.gt.0 ) go to 199

call ibcast(nnspin,1,0)
call ibcast(ncscaleo,1,0)

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )

!-----read eigenvalues
if( lspin ) then

    !--- spin-polarized case
    if( myid.eq.0 ) then

  ioif2: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
      nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

        spindo: do nspin = 1, nnspin
           read(iunit,iostat=istat) eig(1:min(nbando,nband))
           if( istat /= 0 ) exit spindo
           do i = nbando + 1, nband
              eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
           end do
           if( lreadocc ) then
               read(iunit,iostat=istat) occ(1:min(nbando,nband))
               if( istat /= 0 ) exit spindo
               do i = nbando + 1, nband
                  occ(i) = 0.d0
               end do
           end if
           call sveig( nspin, eig, egvlud, nband )
           call sveig( nspin, occ, wegud, nband )
        end do spindo
        if( istat == 0 ) then
            if( nnspin == 1 ) then
                nspin = 2
                call sveig( nspin, eig, egvlud, nband )
                call sveig( nspin, occ, wegud, nband )
            end if
        else
            ierror = 1
        end if

        !-----for divide-and-conquer MD
        if( nfiles > 1 ) then
            allocate(abuf(nbando*2,nnspin))
            abuf(:,:) = 0d0
            do ii=1, nfiles-1
               if( istat == 0 ) read(iunit,iostat=istat) abuf(1:nbando*2,1:nnspin)
               call cdsend(myid_qm_un+ii,abuf,nbando*2*nnspin,myid_qm_un+ii,0)
            end do
            deallocate(abuf)
        end if

        if( istat /= 0 ) ierror = 1

  else ioif2

        !-----for divide-and-conquer MD
        allocate(abuf(nbando*2,nnspin))
        isnd=(myid_qm_un/iogpsz)*iogpsz
        call cdrecvs(myid_qm_un,abuf,nbando*2*nnspin,isnd,0)

        do nspin = 1, nnspin
           eig(1:min(nbando,nband)) = abuf(1:min(nbando,nband),nspin)
           occ(1:min(nbando,nband)) = abuf(nbando+1:nbando+min(nbando,nband),nspin)
           do i = nbando + 1, nband
              eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
              occ(i) = 0.d0
           end do
           call sveig( nspin, eig, egvlud, nband )
           call sveig( nspin, occ, wegud, nband )
        end do
        if( nnspin == 1 ) then
            nspin = 2
            call sveig( nspin, eig, egvlud, nband )
            call sveig( nspin, occ, wegud, nband )
        end if

        deallocate(abuf)

  end if ioif2

    end if

    !-----set communicator
    call get_worldkd( myid, nodes, nkd )

    !--- error trap
    call gimax(ierror)
    if( ierror.gt.0 ) go to 199
    call dbcast(egvlud,nband*2,0)
    call dbcast(wegud,nband*2,0)

else

    !--- spin-unpolarized case
    if( myid.eq.0 ) then

  ioif22: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
      nfiles=min(iogpsz,nodes_qm_un-myid_qm_un)  ! always   1    for no divide-and-conquer MD

        if( ncscale == 1 ) then
            do nspin = 1, nnspin
               if( istat == 0 ) read(iunit,iostat=istat) eig(1:min(nbando,nband))
            end do
            if( istat == 0 ) then
                do i = nbando + 1, nband
                   eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
                end do
            end if
            if( lreadocc ) then
            do nspin = 1, nnspin
               if( istat == 0 ) read(iunit,iostat=istat) occ(1:min(nbando,nband))
            end do
            do i = nbando + 1, nband
               occ(i) = 0.d0
            end do
            end if
        else
            !---noncollinear magnetism
            if( nnspin == 2 ) then
                !---from LSDA calculation
                ! nspin = 1
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( eig(2*i-1), i = 1, min(nbando,nband/2) )
                ! nspin = 1
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( occ(2*i-1), i = 1, min(nbando,nband/2) )
                ! nspin = 2
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( eig(2*i  ), i = 1, min(nbando,nband/2) )
                ! nspin = 2
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( occ(2*i  ), i = 1, min(nbando,nband/2) )

                if( istat == 0 ) then
                    do i = nbando*2 + 1, nband
                       eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
                    end do
                end if
                do i = nbando*2 + 1, nband
                   occ(i) = 0.d0
                end do
            else if( ncscaleo == 1 ) then
                !---from LDA calculation
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( eig(2*i-1), i = 1, min(nbando,nband/2) )
                do i = 1, min(nbando,nband/2)
                   eig(2*i) = eig(2*i-1)
                end do
                if( istat == 0 ) then
                    do i = nbando*2 + 1, nband
                       eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
                    end do
                end if
                if( istat == 0 ) read(iunit,iostat=istat)  &
&                                ( occ(2*i-1), i = 1, min(nbando,nband/2) )
                do i = 1, min(nbando,nband/2)
                   occ(2*i) = occ(2*i-1)
                end do
                do i = nbando*2 + 1, nband
                   occ(i) = 0.d0
                end do
                do i = 1, nband
                   occ(i) = occ(i) * 0.5d0
                end do
            else
                !---from noncollinear magnetism calculation
                if( istat == 0 ) read(iunit,iostat=istat) eig(1:min(nbando,nband))
                if( istat == 0 ) then
                    do i = nbando + 1, nband
                       eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
                    end do
                end if
                if( istat == 0 ) read(iunit,iostat=istat) occ(1:min(nbando,nband))
                do i = nbando + 1, nband
                   occ(i) = 0.d0
                end do
            end if
        end if

        !-----for divide-and-conquer MD
        if( nfiles > 1 ) then
            allocate(abuf(nbando*2,nnspin))
            abuf(:,:) = 0d0
            do ii=1, nfiles-1
               if( istat == 0 ) read(iunit,iostat=istat) abuf(1:nbando*2,1:nnspin)
               call cdsend(myid_qm_un+ii,abuf,nbando*2*nnspin,myid_qm_un+ii,0)
            end do
            deallocate(abuf)
        end if

        if( istat /= 0 ) ierror = 1

  else ioif22

        !-----for divide-and-conquer MD
        allocate(abuf(nbando*2,nnspin))
        isnd=(myid_qm_un/iogpsz)*iogpsz
        call cdrecvs(myid_qm_un,abuf,nbando*2*nnspin,isnd,0)

        if( ncscale == 1 ) then
            eig(1:min(nbando,nband)) = abuf(1:min(nbando,nband),1)
            occ(1:min(nbando,nband)) = abuf(nbando+1:nbando+min(nbando,nband),1)
            do i = nbando + 1, nband
               eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
               occ(i) = 0.d0
            end do
        else
            !---noncollinear magnetism
            if( nnspin == 2 ) then
                !---from LSDA calculation
                eig(1:2*min(nbando,nband/2)-1:2) = abuf(1:min(nbando,nband/2),1)
                occ(1:2*min(nbando,nband/2)-1:2) = abuf(nbando+1:nbando+min(nbando,nband/2),1)
                eig(2:2*min(nbando,nband/2)  :2) = abuf(1:min(nbando,nband/2),2)
                occ(2:2*min(nbando,nband/2)  :2) = abuf(nbando+1:nbando+min(nbando,nband/2),2)
                do i = nbando*2 + 1, nband
                   eig(i) = eig(i-1) + ( eig(2*nbando) - eig(1) )/dble(2*nbando)
                   occ(i) = 0.d0
                end do
            else if( ncscaleo == 1 ) then
                !---from LDA calculation
                eig(1:2*min(nbando,nband/2)-1:2) = abuf(1:min(nbando,nband/2),1)
                occ(1:2*min(nbando,nband/2)-1:2) = abuf(nbando+1:nbando+min(nbando,nband/2),1)
                do i = 1, min(nbando,nband/2)
                   eig(2*i) = eig(2*i-1)
                   occ(2*i) = occ(2*i-1)
                end do
                do i = nbando*2 + 1, nband
                   eig(i) = eig(i-1) + ( eig(2*nbando) - eig(1) )/dble(2*nbando)
                   occ(i) = 0.d0
                end do
                do i = 1, nband
                   occ(i) = occ(i) * 0.5d0
                end do
            else
                !---from noncollinear magnetism calculation
                eig(1:min(nbando,nband)) = abuf(1:min(nbando,nband),1)
                occ(1:min(nbando,nband)) = abuf(nbando+1:nbando+min(nbando,nband),1)
                do i = nbando + 1, nband
                   eig(i) = eig(i-1) + ( eig(nbando) - eig(1) )/dble(nbando)
                   occ(i) = 0.d0
                end do
            end if
        end if

        deallocate(abuf)

  end if ioif22

    end if

    !-----set communicator
    call get_worldkd( myid, nodes, nkd )

    !--- error trap
    call gimax(ierror)
    if( ierror.gt.0 ) go to 199
    call dbcast(eig,nband,0)
    call dbcast(occ,nband,0)

end if

if( ncscaleo == 2 ) then

    !-----set communicator
    call get_worldqmun( myid_qm_un, nodes_qm_un )

    if( ncscale == 1 ) then
        if( myid.eq.0 ) then

  ioif3: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
      nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

            read(iunit,iostat=istat) dummy

            !-----for divide-and-conquer MD
            if( nfiles > 1 ) then
                do ii=1, nfiles-1
                   if( istat == 0 ) read(iunit,iostat=istat) dummy
                end do
            end if

            if( istat /= 0 ) ierror = 1

  end if ioif3

        end if

        !-----set communicator
        call get_worldkd( myid, nodes, nkd )

        !--- error trap
        call gimax(ierror)
        if( ierror.gt.0 ) go to 199

    else

        !---noncollinear magnetism
        if( myid.eq.0 ) then

  ioif32: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
      nfiles=min(iogpsz,nodes_qm_un-myid_qm_un)  ! always   1    for no divide-and-conquer MD

            read(iunit,iostat=istat) wvratio(1:min(nbando,nband)*2)

            !-----for divide-and-conquer MD
            if( nfiles > 1 ) then
                allocate(abuf(nbando,2))
                abuf(:,:) = 0d0
                do ii=1, nfiles-1
                   if( istat == 0 ) read(iunit,iostat=istat) abuf(1:nbando,2)
                   call cdsend(myid_qm_un+ii,abuf,nbando*2,myid_qm_un+ii,0)
                end do
                deallocate(abuf)
            end if

            if( istat /= 0 ) ierror = 1

  else ioif32

            !-----for divide-and-conquer MD
            allocate(abuf(nbando,2))
            isnd=(myid_qm_un/iogpsz)*iogpsz
            call cdrecvs(myid_qm_un,abuf,nbando*2,isnd,0)

            do i = 1, 2
               wvratio(1:min(nbando,nband)*i) = abuf(1:min(nbando,nband),i)
            end do

            deallocate(abuf)

  end if ioif32

        end if

        !-----set communicator
        call get_worldkd( myid, nodes, nkd )

        !--- error trap
        call gimax(ierror)
        if( ierror.gt.0 ) go to 199
        call dbcast(wvratio,nband*2,0)
    end if
end if

csh0 = timecnt()
if( ncscale == 1 ) then
    if( lspin ) then
        nnspin = 2
    else
        nnspin = 1
    end if
    do nspin = 1, nnspin
       call rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& gdcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex, nbando, nbandrec, .true., seed, ierror )

       if( ierror.gt.0 ) then
           if( nspin.eq.1 ) then
               go to 199
             else
               !--- initial w.f. for down-spin
               ierror = 0
               if( .not.lwfrand ) then
                   !--- same as up-spin
                   !--- copy gdcrsv to cgjr ---
                   call cprhcg( nfile, myid, nodes,  &
&                           cgjr, gdcrsv, nplwex, nbnod )
                 else
                   !--- by random numbers
                   lsetnspin12 = .true.
               end if
           end if
       end if
       if( lspin .and. .not.lsetnspin12 ) then
           !--- save wavefunctions decomposed by bands
           call stpsud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
       end if
    end do
else
    !---noncollinear magnetism
    if( nnspin == 2 ) then
        !---from LSDA calculation
        call rdwfn2_from_LSDA( nfile, myid, nodes, iunit, &
& gdcr, nplwex*ncscale, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex, nbando, nbandrec, .true., seed, ierror )
    else if( ncscaleo == 1 ) then
        !---from LDA calculation
        call rdwfn2_from_LDA( nfile, myid, nodes, iunit, &
& gdcr, nplwex*ncscale, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex, nbando, nbandrec, .true., seed, ierror )
    else
        !---from noncollinear magnetism calculation
        if( myid == 0 ) nplwoex2 = nplwoex*ncscale
        call rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& gdcr, nplwex*ncscale, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex2, nbando, nbandrec, .false., seed, ierror )
    end if
end if

199 if( myid.eq.0 ) then
  ioif4: if( mod(myid_qm_un,iogpsz) == 0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
  end if ioif4
end if
if( ierror == 0 ) then
    csh = timecnt()
    if( myid.eq.0 ) then
        if(loutfile(1)) write(nfile(1),*) '                read & distrib. w.f.',  &
    &                         ' : cpu-time :', csh-csh0
        if(loutfile(2)) write(nfile(2),*) '                read & distrib. w.f.',  &
    &                         ' : cpu-time :', csh-csh0
    end if
    csh0 = csh
else
    lstart = .false.
    if( myid.eq.0 ) then
        if(loutfile(1)) write(nfile(1),*) 'error : in files eig000 ...'
        if(loutfile(2)) write(nfile(2),*) 'error : in files eig000 ...'
    end if
end if

!=======================================================================
end if startif


!-----set communicator
call get_worldkd( myid, nodes, nkd )


!---  initial vectors by random numbers --------------------------------
if( .not.lstart .or. lsetnspin12 ) then
if( lspin ) then
    if( .not.lsetnspin12 ) then
        if( .not.lwfrand ) then
            !--- up-spin w.f. are set by random numbers
            nspin1 = 1
            nspin2 = 1
          else
            !--- both w.f. are set by random numbers
            nspin1 = 1
            nspin2 = 2
        end if
      else
        nspin1 = 2
        nspin2 = 2
    end if
  else
    nspin1 = 1
    nspin2 = 1
end if
do nspin = nspin1, nspin2
if( lspin ) then
    if( myid.eq.0 ) then
        if( nspin == 1 ) then
 if(loutfile(1)) write(nfile(1),*) ' initial w.f. by random numbers for up-spin'
 if(loutfile(2)) write(nfile(2),*) ' initial w.f. by random numbers for up-spin'
          else
 if(loutfile(1)) write(nfile(1),*) ' initial w.f. by random numbers for down-spin'
 if(loutfile(2)) write(nfile(2),*) ' initial w.f. by random numbers for down-spin'
        end if
    end if
else
    if( myid.eq.0 ) then
    if(loutfile(1)) write(nfile(1),*) ' initial wave functions by random numbers'
    if(loutfile(2)) write(nfile(2),*) ' initial wave functions by random numbers'
    end if
end if

do idst = 0, nodes - 1

   if( myid.eq.0 ) then
!--- set No. of sending data (bands)
       is1 = nbndsp(idst+1) + 1
       is2 = nbndsp(idst+1) + nbncnt(idst+1)
       is3  = is2 - is1 + 1
       do n = 1, is3
          do i = 1, nplwex*ncscale
             CALL rnd00( rr1, seed )
             gdcr(i,n) = rr1 - 0.5d0
          end do
          gdcr(2,n) = 0.d0
       end do
       if( idst.eq.0 ) then
           do n = 1, is3
              do i = 1, nplwex*ncscale
                 cgjr(i,n) = gdcr(i,n)
              end do
           end do
         else
           nsd = is3*nplwex*ncscale
           call cdsend(100,gdcr,nsd,idst,0)
       end if
   else if( myid.eq.idst ) then
       nrc = nbnod*nplwex*ncscale
       call cdrecv(100,cgjr,nrc,0)
   end if
   call gsync

end do

   if( lspin ) then
!--- save wavefunctions decomposed by bands
       call stpsud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
   end if
end do
if( lspin .and. .not.lwfrand ) then
    !--- same as up-spin
    !--- copy gdcrsv to cgjr ---
    call cprhcg( nfile, myid, nodes,  &
&                       cgjr, gdcrsv, nplwex, nbnod )
    nspin = 2
    call stpsud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
end if
!-----------------------------------------------------------------------
end if


ct = timecnt()
if( myid.eq.0 ) then
    if(loutfile(1)) write(nfile(1),*) '                   set wavefunctions',  &
&                         ' : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '                   set wavefunctions',  &
&                         ' : cpu-time :', ct-ct0
end if
ct0 = ct

end if nkdif


!-----set communicator
call get_worldun( myid, nodes )


call dbcast(cgjr,nbnod*nplwex*ncscale,0)
call dbcast(eig,nband,0)
if( lspin ) then
    call dbcast(egvlud,nband*2,0)
    call dbcast(gdcrsv,nspnod,0)
    call lbcast(lcgjsv,1,0)
end if
if( ncscale == 2 ) then
    !---noncollinear magnetism
    call dbcast(wvratio,nband*2,0)
end if

end if npwif


if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync

    call dbcast(cgjr,nbnod*nplwex*ncscale,0)
    call dbcast(eig,nband,0)
    if( lspin ) then
        call dbcast(egvlud,nband*2,0)
        call dbcast(gdcrsv,nspnod,0)
        call lbcast(lcgjsv,1,0)
    end if
    if( ncscale == 2 ) then
        !---noncollinear magnetism
        call dbcast(wvratio,nband*2,0)
    end if
end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync

    call dbcast(cgjr,nbnod*nplwex*ncscale,0)
    call dbcast(eig,nband,0)
    if( lspin ) then
        call dbcast(egvlud,nband*2,0)
        call dbcast(gdcrsv,nspnod,0)
        call lbcast(lcgjsv,1,0)
    end if
    if( ncscale == 2 ) then
        !---noncollinear magnetism
        call dbcast(wvratio,nband*2,0)
    end if
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




subroutine rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, lsifreal_, nplwoex, nbando, nbandrec, lgammak, seed, ierror )
!-----------------------------------------------------------------------
!     read wave functions : rhcr
!
!    read data on node 0, and distribute to other nodes
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iogpsz
integer :: iunit
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: cgjr(nplwex,nbnod)
real*8  :: rhcr(nplwex,nbnod)
logical :: lsifreal_, lgammak
integer :: nplwoex, nbando, nbandrec
real*8  :: seed
integer :: ierror

!-----declare local variables
integer :: idst, is1, is2, is3, iss3, n, i, ig, nrc, nsd
real*8  :: efactr, rr1
integer :: istat
integer :: nkd, myid_qm_un, nodes_qm_un
integer :: nfiles, ii, isnd, iogp


!-----set communicator
!-----to get myid_qm_un, nodes_qm_un
call get_worldqmun( myid_qm_un, nodes_qm_un )

!-----set communicator
call get_worldkd( myid, nodes, nkd )


if( myid == 0 ) lsifreal = lsifreal_
call lbcast(lsifreal,1,0) 
if( .not.lsifreal ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < nplwoex ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(nplwoex), stat=ierror )
        end if
    end if
    call ibcast(ierror,1,0) 
    if( ierror /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in rdwfn2', ierror
            write(nfile(2),*) '*** memory allocation error in rdwfn2', ierror
        end if
        return
    end if
end if


ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD

    do idst = 0, nodes - 1
       if( myid == 0 ) then
           !--- set No. of sending data (bands)
           is1 = nbndsp(idst+1) + 1
           is2 = nbndsp(idst+1) + nbncnt(idst+1)
           is3  = is2 - is1 + 1
           iss3 = min( is2, nbando ) - is1 + 1
           cgjr(1:nplwex,1:is3) = 0.d0

           if( lsifreal ) then
               !--- read data in real*8
               do n = 1, iss3
                  read(iunit,iostat=istat) cgjr(1:nplwoex,n)
                  if( istat /= 0 ) exit
               end do
           else
               !--- read data in integer*2
               do n = 1, iss3
                  read(iunit,iostat=istat) efactr
                  if( istat /= 0 ) exit
                  read(iunit,iostat=istat) saveigv(1:nplwoex)
                  if( istat /= 0 ) exit
                  cgjr(1:nplwoex,n) = saveigv(1:nplwoex) * efactr
               end do
           end if
           if( istat == 0 ) then
               if( nbando.lt.is2 ) then
                   do n = max( is1, nbando+1 ) - is1 + 1, is3
                      do ig = 1, nplwex
                         CALL rnd00( rr1, seed )
                         cgjr(ig,n) = rr1 - 0.5d0
                      end do
                      if( lgammak ) cgjr(2,n) = 0.d0
                   end do
               end if
           else
               ierror = 1
           end if

           if( idst == 0 ) then
               rhcr(1:nplwex,1:is3) = cgjr(1:nplwex,1:is3)
           else
               nsd = is3*nplwex
               if( nsd > 0 ) call cdsend(100,cgjr,nsd,idst,0)
           end if

       else if( myid == idst ) then
           nrc = nbnod*nplwex
           if( nrc > 0 ) call cdrecv(100,rhcr,nrc,0)
       end if

!   call gsync
       call gimax(ierror)
!   if( ierror /=0 ) return
       if( ierror /=0 ) exit
    end do

    if( myid == 0 ) then
        if( lsifreal ) then
            !--- read data in real*8
            do n = nbando + 1, nbandrec
               read(iunit,iostat=istat) efactr
               if( istat /= 0 ) exit
            end do
        else
            !--- read data in integer*2
            do n = nbando + 1, nbandrec
               read(iunit,iostat=istat) efactr
               if( istat /= 0 ) exit
               read(iunit,iostat=istat) saveigv(1:nplwoex)
               if( istat /= 0 ) exit
            end do
        end if
    end if

end if ioif


return
end




subroutine rdwfn2( nfile, myid, nodes, iunit, &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, lsifreal_, nplwoex, nbando, nbandrec, lgammak, seed, ierror )
!-----------------------------------------------------------------------
!     read wave functions : rhcr
!
!    read data on node 0, and distribute to other nodes
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: cgjr(nplwex,nbnod)
real*8  :: rhcr(nplwex,nbnod)
logical :: lsifreal_, lgammak
integer :: nplwoex, nbando, nbandrec
real*8  :: seed
integer :: ierror

!-----declare local variables
integer :: idst, is1, is2, is3, iss3, n, i, ig, nrc, nsd
real*8  :: efactr, rr1
integer :: istat


if( myid == 0 ) lsifreal = lsifreal_
call lbcast(lsifreal,1,0) 
if( .not.lsifreal ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < nplwoex ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(nplwoex), stat=ierror )
        end if
    end if
    call ibcast(ierror,1,0) 
    if( ierror /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in rdwfn2', ierror
            write(nfile(2),*) '*** memory allocation error in rdwfn2', ierror
        end if
        return
    end if
end if

do idst = 0, nodes - 1
   if( myid == 0 ) then
!--- set No. of sending data (bands)
       is1 = nbndsp(idst+1) + 1
       is2 = nbndsp(idst+1) + nbncnt(idst+1)
       is3  = is2 - is1 + 1
       iss3 = min( is2, nbando ) - is1 + 1
       cgjr(1:nplwex,1:is3) = 0.d0

       if( lsifreal ) then
!--- read data in real*8
           do n = 1, iss3
              read(iunit,iostat=istat) cgjr(1:nplwoex,n)
              if( istat /= 0 ) exit
           end do
         else
!--- read data in integer*2
           do n = 1, iss3
              read(iunit,iostat=istat) efactr
              if( istat /= 0 ) exit
              read(iunit,iostat=istat) saveigv(1:nplwoex)
              if( istat /= 0 ) exit
              cgjr(1:nplwoex,n) = saveigv(1:nplwoex) * efactr
           end do
       end if
       if( istat == 0 ) then
           if( nbando.lt.is2 ) then
               do n = max( is1, nbando+1 ) - is1 + 1, is3
                  do ig = 1, nplwex
                     CALL rnd00( rr1, seed )
                     cgjr(ig,n) = rr1 - 0.5d0
                  end do
                  if( lgammak ) cgjr(2,n) = 0.d0
               end do
           end if
       else
           ierror = 1
       end if

       if( idst == 0 ) then
           rhcr(1:nplwex,1:is3) = cgjr(1:nplwex,1:is3)
         else
           nsd = is3*nplwex
           if( nsd > 0 ) call cdsend(100,cgjr,nsd,idst,0)
       end if

   else if( myid == idst ) then
       nrc = nbnod*nplwex
       if( nrc > 0 ) call cdrecv(100,rhcr,nrc,0)
   end if

!   call gsync
   call gimax(ierror)
   if( ierror /=0 ) return
end do

if( myid == 0 ) then
    if( lsifreal ) then
!--- read data in real*8
        do n = nbando + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
        end do
    else
!--- read data in integer*2
        do n = nbando + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
           read(iunit,iostat=istat) saveigv(1:nplwoex)
           if( istat /= 0 ) exit
        end do
    end if
end if


return
end




subroutine rdwfn2_from_LDA( nfile, myid, nodes, iunit, &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, lsifreal_, nplwoex, nbando, nbandrec, lgammak, seed, ierror )
!-----------------------------------------------------------------------
!     read wave functions : rhcr from LDA
!
!    read data on node 0, and distribute to other nodes
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: cgjr(nplwex,nbnod)
real*8  :: rhcr(nplwex,nbnod)
logical :: lsifreal_, lgammak
integer :: nplwoex, nbando, nbandrec
real*8  :: seed
integer :: ierror

!-----declare local variables
integer :: idst, is1, is2, is3, iss3, n, i, ig, nrc, nsd, nplwo
integer :: igo, igor, igoi, igr, igi, nread
real*8  :: efactr, rr1
integer :: istat


if( myid == 0 ) lsifreal = lsifreal_
call lbcast(lsifreal,1,0) 
if( .not.lsifreal ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < nplwoex ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(nplwoex), stat=ierror )
        end if
    end if
    call ibcast(ierror,1,0) 
    if( ierror /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in rdwfn2', ierror
            write(nfile(2),*) '*** memory allocation error in rdwfn2', ierror
        end if
        return
    end if
end if

nread = 0
do idst = 0, nodes - 1
   if( myid == 0 ) then
!--- set No. of sending data (bands)
       is1 = nbndsp(idst+1) + 1
       if( mod(is1,2) == 0 ) then
           cgjr(nplwex/2+1:nplwex,1) = cgjr(1:nplwex/2,is3)
           cgjr(1:nplwex/2,1) = 0.d0
       else
           cgjr(1:nplwex,1) = 0.d0
       end if
       is2 = nbndsp(idst+1) + nbncnt(idst+1)
       is3  = is2 - is1 + 1
!       iss3 = min( is2, nbando ) - is1 + 1
       cgjr(1:nplwex,2:is3) = 0.d0

       !---read odd band index
!       do n = 1, iss3
       do n = is1+1-mod(is1,2), is1 + is3 - 1, 2
          nread = nread + 1
          if( nread <= nbandrec ) then
              if( lsifreal ) then
                  !--- read data in real*8
                  read(iunit,iostat=istat) cgjr(1:nplwoex,n-is1+1)
              else
                  !--- read data in integer*2
                  read(iunit,iostat=istat) efactr
                  if( istat /= 0 ) exit
                  read(iunit,iostat=istat) saveigv(1:nplwoex)
                  if( istat /= 0 ) exit
                  cgjr(1:nplwoex,n-is1+1) = saveigv(1:nplwoex) * efactr
              end if
          else
              do ig = 1, nplwoex
                 CALL rnd00( rr1, seed )
                 cgjr(ig,n-is1+1) = rr1 - 0.5d0
              end do
!              if( lgammak ) cgjr(2,n) = 0.d0
          end if
          if( istat /= 0 ) exit
          if( lgammak ) then
              nplwo = nplwoex/2 - 1
              do igo = nplwo, 2, -1
                 igor = 2*igo + 1
                 igoi = 2*igo + 2
                 ig  = 2*igo - 1
                 igr = 2*ig + 1
                 igi = 2*ig + 2
                 cgjr(igr:igi,n-is1+1) = cgjr(igor:igoi,n-is1+1)
              end do
              do igo = 1, nplw, 2
                 igor = 2*igo + 1
                 igoi = 2*igo + 2
                 ig  = igo + 1
                 igr = 2*ig + 1
                 igi = 2*ig + 2
                 cgjr(igr,n-is1+1) =  cgjr(igor,n-is1+1)
                 cgjr(igi,n-is1+1) = -cgjr(igoi,n-is1+1)
              end do
          end if
          if( n + 1 <= is2 ) then
              cgjr(1:nplwex,n-is1+2) = 0.d0
              cgjr(nplwex/2+1:nplwex,n-is1+2) = cgjr(1:nplwex/2,n-is1+1)
          end if
       end do
       if( istat /= 0 ) then
           ierror = 1
       end if

       if( idst == 0 ) then
           rhcr(1:nplwex,1:is3) = cgjr(1:nplwex,1:is3)
         else
           nsd = is3*nplwex
           if( nsd > 0 ) call cdsend(100,cgjr,nsd,idst,0)
       end if

   else if( myid == idst ) then
       nrc = nbnod*nplwex
       if( nrc > 0 ) call cdrecv(100,rhcr,nrc,0)
   end if

!   call gsync
   call gimax(ierror)
   if( ierror /=0 ) return
end do

if( myid == 0 ) then
    if( lsifreal ) then
!--- read data in real*8
        do n = nread + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
        end do
    else
!--- read data in integer*2
        do n = nread + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
           read(iunit,iostat=istat) saveigv(1:nplwoex)
           if( istat /= 0 ) exit
        end do
    end if
end if


return
end




subroutine rdwfn2_from_LSDA( nfile, myid, nodes, iunit, &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, lsifreal_, nplwoex, nbando, nbandrec, lgammak, seed, ierror )
!-----------------------------------------------------------------------
!     read wave functions : rhcr from LDA
!
!    read data on node 0, and distribute to other nodes
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: cgjr(nplwex,nbnod)
real*8  :: rhcr(nplwex,nbnod)
logical :: lsifreal_, lgammak
integer :: nplwoex, nbando, nbandrec
real*8  :: seed
integer :: ierror

!-----declare local variables
integer :: idst, is1, is2, is3, iss3, n, i, ig, nrc, nsd, nplwo
integer :: igo, igor, igoi, igr, igi, nnsd, noddeven, n0, nread
real*8  :: efactr, rr1
integer :: istat


if( myid == 0 ) lsifreal = lsifreal_
call lbcast(lsifreal,1,0) 
if( .not.lsifreal ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < nplwoex ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(nplwoex), stat=ierror )
        end if
    end if
    call ibcast(ierror,1,0) 
    if( ierror /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in rdwfn2', ierror
            write(nfile(2),*) '*** memory allocation error in rdwfn2', ierror
        end if
        return
    end if
end if

oddevendo: do noddeven = 1, 2

nread = 0
do idst = 0, nodes - 1
   !--- set No. of sending data (bands)
   is1 = nbndsp(idst+1) + 1
   is2 = nbndsp(idst+1) + nbncnt(idst+1)
   is3  = is2 - is1 + 1
!   iss3 = min( is2, nbando ) - is1 + 1
   if( noddeven == 1 ) then
       n0 = is1+1-mod(is1,2)
   else
       n0 = is1+mod(is1,2)
   end if
   if( myid == 0 ) then
       cgjr(1:nplwex,1:is3) = 0.d0

       !---read odd band index
!       do n = 1, iss3
       nnsd = 0
       do n = n0, is1 + is3 - 1, 2
          nnsd = nnsd + 1
          nread = nread + 1
          if( nread <= nbandrec ) then
              if( lsifreal ) then
                  !--- read data in real*8
                  read(iunit,iostat=istat) cgjr(1:nplwoex,nnsd)
              else
                  !--- read data in integer*2
                  read(iunit,iostat=istat) efactr
                  if( istat /= 0 ) exit
                  read(iunit,iostat=istat) saveigv(1:nplwoex)
                  if( istat /= 0 ) exit
                  cgjr(1:nplwoex,nnsd) = saveigv(1:nplwoex) * efactr
              end if
          else
              do ig = 1, nplwoex
                 CALL rnd00( rr1, seed )
                 cgjr(ig,nnsd) = rr1 - 0.5d0
              end do
!              if( lgammak ) cgjr(2,nnsd) = 0.d0
          end if
          if( istat /= 0 ) exit
          if( lgammak ) then
              nplwo = nplwoex/2 - 1
              do igo = nplwo, 2, -1
                 igor = 2*igo + 1
                 igoi = 2*igo + 2
                 ig  = 2*igo - 1
                 igr = 2*ig + 1
                 igi = 2*ig + 2
                 cgjr(igr:igi,nnsd) = cgjr(igor:igoi,nnsd)
              end do
              do igo = 1, nplw, 2
                 igor = 2*igo + 1
                 igoi = 2*igo + 2
                 ig  = igo + 1
                 igr = 2*ig + 1
                 igi = 2*ig + 2
                 cgjr(igr,nnsd) =  cgjr(igor,nnsd)
                 cgjr(igi,nnsd) = -cgjr(igoi,nnsd)
              end do
          end if
       end do
       if( istat /= 0 ) then
           ierror = 1
       end if

       if( idst == 0 ) then
           nnsd = 0
!           rhcr(:,:) = 0.d0
           do n = n0, is2, 2
              nnsd = nnsd + 1
              if( noddeven == 1 ) then
                  rhcr(         1:nplwex/2,n-is1+1) = cgjr(1:nplwex/2,nnsd)
                  rhcr(nplwex/2+1:nplwex,  n-is1+1) = 0.d0
              else
                  rhcr(         1:nplwex/2,n-is1+1) = 0.d0
                  rhcr(nplwex/2+1:nplwex,  n-is1+1) = cgjr(1:nplwex/2,nnsd)
              end if
           end do
         else
           call cisend(10,nnsd,1,idst,0)
           nsd = nnsd*nplwex
           if( nsd > 0 ) call cdsend(100,cgjr,nsd,idst,0)
       end if

   else if( myid == idst ) then
       call cirecv(10,nnsd,1,0)
       nrc = nnsd*nplwex
       if( nrc > 0 ) call cdrecv(100,cgjr,nrc,0)
       nnsd = 0
!       rhcr(:,:) = 0.d0
       do n = n0, is2, 2
          nnsd = nnsd + 1
          if( noddeven == 1 ) then
              rhcr(         1:nplwex/2,n-is1+1) = cgjr(1:nplwex/2,nnsd)
              rhcr(nplwex/2+1:nplwex,  n-is1+1) = 0.d0
          else
              rhcr(         1:nplwex/2,n-is1+1) = 0.d0
              rhcr(nplwex/2+1:nplwex,  n-is1+1) = cgjr(1:nplwex/2,nnsd)
          end if
       end do
   end if

!   call gsync
   call gimax(ierror)
   if( ierror /=0 ) return
end do

if( myid == 0 ) then
    if( lsifreal ) then
!--- read data in real*8
        do n = nread + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
        end do
    else
!--- read data in integer*2
        do n = nread + 1, nbandrec
           read(iunit,iostat=istat) efactr
           if( istat /= 0 ) exit
           read(iunit,iostat=istat) saveigv(1:nplwoex)
           if( istat /= 0 ) exit
        end do
    end if
end if
end do oddevendo


return
end




subroutine rdcdty( nfile, iogpsz, ct0,  &
& fname1, lstart, hcell, h_MD, lorthrhmbc, rho, rdelv, rdel, nel, multg,  &
& lclust, nd1v, ntype, nhk1, nhk2, zv, lclno, ratm, natom, mx1,  &
& mshnod, mshnx, mshny, mshnz, mshx1, mshy1, mshz1,  &
& mshx, mshy, mshz, mulpms, noddatx,  &
& lspin, nspnmx, rhoud, diffud, lfixud,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& rhgr, rinr, nplw5ex, nplw5, fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& nga, ngb, ngc, lnoncollinear, keypcd, ncprvmx )
!-----------------------------------------------------------------------
!     initial electron density : rho
!-----------------------------------------------------------------------
use outfile
use ncmagne_variables
implicit none
integer :: nfile(*)
integer :: iogpsz
real*8  :: ct0
character(*) :: fname1
logical :: lstart
real*8  :: hcell(3,3), h_MD(3,3)
logical :: lorthrhmbc
real*8  :: rho(*), rdelv, rdel(3)
integer :: nel, multg
logical :: lclust
integer, dimension(3) :: nd1v
integer :: ntype, natom
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: lclno
real*8,  dimension(3,natom)   :: ratm
integer :: mx1
integer :: mshnod(*), mshnx(*), mshny(*), mshnz(*), mshx1(*), mshy1(*), mshz1(*),  &
& mshx(*), mshy(*), mshz(*), mulpms(*), noddatx
logical :: lspin
integer :: nspnmx
real*8  :: rhoud(*), diffud
logical :: lfixud
real*8  :: glocal(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex,*), rinr(nplw5ex,*)
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
integer, dimension(0:nplw5) :: nga, ngb, ngc
logical :: lnoncollinear
integer :: keypcd, ncprvmx  ! for divide-and-conquer MD

!-----declare local variables
integer :: myid, nodes, myid_kd, nodes_kd, nkd,  &
& myid_lr, nodes_lr, myid_pw, nodes_pw
integer :: ig, ierror, ierro2, ierro3, nplw5o, nplw5oex, nspnmxo, m1
integer :: kfft0o, j
integer :: iunit, nversion, istat
logical :: lspinini
logical :: lnoncollinearo, lreadspd, lncmagini
real*8  :: telcor, dbuf1r
real*8  :: ct00, ct, timecnt


!-----set communicator
call get_worldqm( myid, nodes )


ct00 = ct0
if( myid == 0 ) then
    if(loutfile(1)) write(nfile(1),*) ' initial electron density'
    if(loutfile(2)) write(nfile(2),*) ' initial electron density'
end if

do ig = 1, nplw5ex
   rhgr(ig,1) = 0.d0
end do
nspnmxo = 0
lnoncollinearo = .false.
lreadspd = .false.
lstartif: if( lstart ) then


    ierror = 0
    ierro2 = 0
    ierro3 = 0
if( myid == 0 ) then

    call allocate_unit_number( iunit )
!=======================================================================
    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in rdcdty): ', trim(fname1)
        write(nfile(2),*) 'open file(in rdcdty): ', trim(fname1)
    end if

    open(iunit,iostat=istat,file=fname1,status='old',form='unformatted')
    if( istat == 0 ) read(iunit,iostat=istat) nplw5o
    nversion = 0
    if( istat == 0 .and. nplw5o < 0 ) then
        nversion = abs(nplw5o)
        if( istat == 0 ) read(iunit,iostat=istat) nplw5o
    end if
    if( nversion >= 3 ) then
        if( istat == 0 ) read(iunit,iostat=istat) nspnmxo
        if( istat == 0 ) read(iunit,iostat=istat) lnoncollinearo
    end if
    if( istat == 0 ) then
        if( nplw5o.ne.nplw5 .and. myid.eq.0 ) then
            ierror = 1
            write(nfile(1),*) 'error - different nplw5 in ', trim(fname1)
            write(nfile(2),*) 'error - different nplw5 in ', trim(fname1)
        end if

        nplw5o = min( nplw5o, nplw5 )
        nplw5oex = 2*( nplw5o + 1 )
        read(iunit,iostat=istat) ( rhgr(m1,1), m1 = 1, nplw5oex )
    end if

    if( istat == 0 ) then
        if( lspin .or. lnoncollinear .and. nspnmxo == 0  &
&                 .or. lnoncollinear .and. nspnmxo == 2 ) then
!            ierro2 = 0
            read(iunit,iostat=istat) ( rhgr(m1,2), m1 = 1, nplw5oex )
            lreadspd = istat == 0
            if( .not.lreadspd ) then
                ierro2 = 1
                !--- from LDA ---
                do m1 = 1, nplw5oex
                   rhgr(m1,2) = rhgr(m1,1) * diffud/dble(nel)
                end do
            end if
        end if
    else
       ierror = 1
    end if
!=======================================================================
end if

!--- error trap
    call gimax(ierror)

    !---noncollinear magnetism
!    call ibcast(nspnmxo,1,0)
    call lbcast(lnoncollinearo,1,0)
    call lbcast(lreadspd,1,0)
!    if( ierror == 0 ) then
!    if( lnoncollinear .and. lnoncollinearo ) then
!
!        !-----set communicator
!        call get_worldlr( myid_lr, nodes_lr )
!
!        !-----set communicator
!        call get_worldpw( myid_pw, nodes_pw )
!
!        !-----set communicator
!        call get_worldkd( myid_kd, nodes_kd, nkd )
!
!        nlrif: if( myid_lr == 0 ) then
!        npwif: if( myid_pw == 0 ) then
!        nkdif: if( nkd == 0 ) then
!
!           if( myid == 0 ) then
!               read(iunit,iostat=istat) kfft0o
!               if( istat == 0 ) then
!                   if( kfft0o /= kfft0 .and. myid.eq.0 ) then
!                       ierro3 = 1
!                       write(nfile(1),*) 'error - different kfft0 in ', trim(fname1)
!                       write(nfile(2),*) 'error - different kfft0 in ', trim(fname1)
!                   end if
!               else
!                   ierro3 = 1
!               end if
!           end if
!           !--- error trap
!           call gimax(ierro3)
!           ier0: if( ierro3 == 0 ) then
!
!           !---rhomx
!           call rd_local_data( nfile, myid, nodes, iunit, ierro3,  &
!& rhomx, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0, fft3x )
!
!           ier1: if( ierro3 == 0 ) then
!
!           !---rhomy
!           call rd_local_data( nfile, myid, nodes, iunit, ierro3,  &
!& rhomy, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0, fft3x )
!
!           ier2: if( ierro3 == 0 ) then
!
!           !---rhomz
!           call rd_local_data( nfile, myid, nodes, iunit, ierro3,  &
!& rhomz, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0, fft3x )
!
!           ier3: if( ierro3 == 0 ) then
!
!               if( myid == 0 ) then
!                   write(nfile(1),*) 'read noncollinear magnetism : successful'
!                   write(nfile(2),*) 'read noncollinear magnetism : successful'
!               end if
!
!           end if ier3
!           end if ier2
!           end if ier1
!           end if ier0
!
!        end if nkdif
!
!        !-----set communicator
!        call get_worldun( myid, nodes )
!
!        !--- error trap
!        call gimax(ierro3)
!        if( ierro3 == 0 ) then
!            call dbcast(rhomx,mshnod(1),0)
!            call dbcast(rhomy,mshnod(1),0)
!            call dbcast(rhomz,mshnod(1),0)
!        end if
!
!        end if npwif
!
!
!        if( nodes_pw > 1 ) then
!            !-----set communicator
!            call get_worldpw( myid, nodes )
!
!            !--- error trap
!            call gimax(ierro3)
!            if( ierro3 == 0 ) then
!                call dbcast(rhomx,mshnod(1),0)
!                call dbcast(rhomy,mshnod(1),0)
!                call dbcast(rhomz,mshnod(1),0)
!            end if
!
!        end if
!
!        end if nlrif
!
!        if( nodes_lr > 1 ) then
!            !-----set communicator
!            call get_worldlr( myid, nodes )
!
!            !--- error trap
!            call gimax(ierro3)
!            if( ierro3 == 0 ) then
!                call dbcast(rhomx,mshnod(1),0)
!                call dbcast(rhomy,mshnod(1),0)
!                call dbcast(rhomz,mshnod(1),0)
!            end if
!
!        end if
!
!        !-----set communicator
!        call get_worldqm( myid, nodes )
!
!    end if
!    end if

if( myid == 0 ) then
!=======================================================================
    close(iunit)
!=======================================================================
    call deallocate_unit_number( iunit )
end if

    if( ierror == 0 ) then
        ct = timecnt()
        if( myid.eq.0 ) then
        write(nfile(1),*) '          read c.d. : cpu-time :', ct - ct0
        write(nfile(2),*) '          read c.d. : cpu-time :', ct - ct0
        end if
        ct0 = ct

        !--- broadcast charge density to all nodes
        call dbcast(rhgr,nplw5ex*nspnmx,0)
        ct = timecnt()
        if( myid.eq.0 ) then
        write(nfile(1),*) '     broadcast c.d. : cpu-time :', ct - ct0
        write(nfile(2),*) '     broadcast c.d. : cpu-time :', ct - ct0
        end if
        ct0 = ct

        call gimax(ierro2)

        !-----set communicator
        call get_worldkd( myid, nodes, nkd )

        if( lspin ) then
        if( ierro2 /= 0 ) then
            !--- initial spin density
            !---if inispin /= 1
            call spindenini( nfile, myid, nodes,  &
& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
& lorthrhmbc, rhoud, natom, mx1,  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1), lspinini )

            if( lspinini ) then
                call chknud( nfile, myid, nodes,  &
& rhoud, mshnod(1), rdelv, diffud, lfixud, noddatx, .true. )

               !--- unify spin charge density
               call unifylc( nfile, myid, nodes,  &
& rhoud, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
               ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' unify spin c.d.    : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) ' unify spin c.d.    : cpu-time :', ct - ct0
               ct0 = ct

               !--- charge density transformation from r- to g-spaces
               call chgr2g( nfile, myid, nodes,  &
& glocal, rhgr(1,2), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
               ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' spn-cd.tr.(r to g) : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) ' spn-cd.tr.(r to g) : cpu-time :', ct - ct0
               ct0 = ct

            end if

        end if
        end if

    else
        lstart = .false.
        if( myid == 0 ) then
            write(nfile(1),*) 'error : in files cds000 ...'
            write(nfile(2),*) 'error : in files cds000 ...'
        end if
    end if

end if lstartif


lstartif2: if( .not.lstart ) then
    !-----set communicator
    call get_worldkd( myid, nodes, nkd )
    !=======================================================================
    if(loutfile(1)) write(nfile(1),*)  &
&   ' Initial density is calculated from atomic charge density.'
    if(loutfile(2)) write(nfile(2),*)  &
&   ' Initial density is calculated from atomic charge density.'
    !--- initial electron density
    call denini( nfile, myid, nodes,                                   &
& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
& lorthrhmbc, rho, natom, mx1,  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )

    if( lspin ) then
!        do m1 = 1, mshnod(1)
!           rhoud(m1) = rho(m1) * diffud/dble(nel)
!        end do
        !--- initial spin density
        !---if inispin == 1
        call spindenini1( nfile, myid, nodes,  &
& rhoud, rho, mshnod(1), diffud, nel )
        !---if inispin /= 1
        call spindenini( nfile, myid, nodes,  &
& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
& lorthrhmbc, rhoud, natom, mx1,  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1), lspinini )
    end if

    !--- check No. of electrons
    call chknel( nfile, myid, nodes,  &
& rho, mshnod(1), rdelv, nel, noddatx, .true. )

    if( lspin ) then
        call chknud( nfile, myid, nodes,  &
& rhoud, mshnod(1), rdelv, diffud, lfixud, noddatx, .true. )
    end if

    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '                    : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '                    : cpu-time :', ct - ct0
    ct0 = ct

    !--- unify charge density
    call unifylc( nfile, myid, nodes,  &
& rho, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' unify chg. density : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) ' unify chg. density : cpu-time :', ct - ct0
    ct0 = ct
end if lstartif2


lstartif3: if( .not.lstart ) then
    !--- charge density transformation from r- to g-spaces
    call chgr2g( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '    c.d.tr.(r to g) : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '    c.d.tr.(r to g) : cpu-time :', ct - ct0
    ct0 = ct

    if( lspin ) then
    !--- unify spin charge density
        call unifylc( nfile, myid, nodes,  &
& rhoud, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) ' unify spin c.d.    : cpu-time :', ct - ct0
        if(loutfile(2)) write(nfile(2),*) ' unify spin c.d.    : cpu-time :', ct - ct0
        ct0 = ct

    !--- charge density transformation from r- to g-spaces
        call chgr2g( nfile, myid, nodes,  &
& glocal, rhgr(1,2), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) ' spn-cd.tr.(r to g) : cpu-time :', ct - ct0
        if(loutfile(2)) write(nfile(2),*) ' spn-cd.tr.(r to g) : cpu-time :', ct - ct0
        ct0 = ct
    end if
!=======================================================================
end if lstartif3

!-----Check symmetry of charge density
!call csymcd( nfile, myid, nodes, .true.,  &
!& rhgr, nd1vks, nplw5, nplw5ex, kfft0, nga, ngb, ngc, ijkgd, 1 )

!--- charge density transformation from g- to r-spaces
call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '    c.d.tr.(g to r) : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) '    c.d.tr.(g to r) : cpu-time :', ct - ct0
ct0 = ct

!--- check No. of electrons and charge density
    call chkchg2( nfile, myid, nodes,  &
& glocal, rdelv, nel, ntotfd, telcor )
    do j = 1, nplw5ex
       rhgr(j,1)  = rhgr(j,1) * telcor
    end do

!--- store local charge density
call distlc( nfile, myid, nodes,  &
& rho, mshnod(1), glocal, mftdsp )

if( lspin ) then
    !-----Check symmetry of charge density
!    call csymcd( nfile, myid, nodes, .true.,  &
!& rhgr(1,2), nd1vks, nplw5, nplw5ex, kfft0, nga, ngb, ngc, ijkgd, 2 )
!--- spin charge density transformation from g- to r-spaces
    call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr(1,2), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' spn-cd.tr.(g to r) : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) ' spn-cd.tr.(g to r) : cpu-time :', ct - ct0
    ct0 = ct

!--- check No. of electrons and charge density
!          call chkcud( nfile, myid, nodes,  &
!     & glocal, rdelv, diffud, lfixud, ntotfd )
   call chkcud2( nfile, myid, nodes,  &
& glocal, rdelv, diffud, lfixud, ntotfd, telcor )
    if( lfixud .and. abs(telcor).gt.1.d-10 ) then
        telcor = dble(diffud)/telcor
        do j = 1, nplw5ex
           rhgr(j,2)  = rhgr(j,2) * telcor
        end do
    end if

!--- store local charge density
    call distlc( nfile, myid, nodes,  &
& rhoud, mshnod(1), glocal, mftdsp )
end if

!--- store rinr as input charge density
rinr(1:nplw5ex,1:nspnmx) = rhgr(1:nplw5ex,1:nspnmx)

!--- check No. of electrons
call chknel( nfile, myid, nodes,  &
& rho, mshnod(1), rdelv, nel, noddatx, .true. )

if( lspin ) then
    call chknud( nfile, myid, nodes,  &
& rhoud, mshnod(1), rdelv, diffud, lfixud, noddatx, .true. )
end if


!---set noncollinear magnetic moment
!if( lnoncollinear ) then
!if( .not.lstart .or. .not.lnoncollinearo .or. ierro3 /= 0 ) then
!    lncmagini = .false.
!    if( lreadspd ) then
!        !-----set communicator
!        call get_worldqm( myid, nodes )
!        !--- broadcast charge density to all nodes
!        call dbcast(rhgr(1,2),nplw5ex,0)
!        !-----set communicator
!        call get_worldkd( myid, nodes, nkd )
!
!        !---if inispin == 0, use spindensity as rhomz
!        call ncmagneini0( nfile, myid, nodes,  &
!& rhomx, rhomy, rhomz, mshnod(1), glocal, rhgr(1,2), nplw5ex, nplw5,  &
!& mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, lncmagini )
!    end if
!
!    if( .not.lncmagini ) then
!        !---if inispin == 0, use spindensity as rhomz
!        call ncmagneini( nfile, myid, nodes,  &
!& lclust, nd1v, ntype, nhk1, nhk2, zv, ratm, hcell, rdel,  &
!& lorthrhmbc, rhomx, rhomy, rhomz, natom, mx1,  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1), lncmagini )
!
!    end if
!end if
!
!    call fitmagnepw( nfile, myid, nodes, ct0,  &
!& rhomx, rhomy, rhomz, mshnod(1), glocal,  &
!& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& rhgr, nplw5ex, nplw5, fft3x, fft3y, fftwork, ijkgd, ngenh,  &
!& nga, ngb, ngc )
!
!    !---check
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomx( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomx', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomx', telcor * rdelv
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomy( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomy', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomy', telcor * rdelv
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomz( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomz', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomz', telcor * rdelv
!
!    !--- store rinr as input magnetic moment
!    rinr(1:nplw5ex,2:4) = rhgr(1:nplw5ex,2:4)
!
!    !---save initial magnetic moment
!!    call set_magini( nfile, myid, nodes, mshnod(1) )
!end if


ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '     set c.d. total : cpu-time :', ct - ct00
if(loutfile(2)) write(nfile(2),*) '     set c.d. total : cpu-time :', ct - ct00
ct0 = ct


return
end




subroutine rd_local_data( nfile, myid, nodes, iunit, ierro3,  &
& rhomx, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!-----------------------------------------------------------------------
!    read rhomx at node 0 and distribute it to other nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit, ierro3
integer :: mshnod
real*8  :: rhomx(mshnod), glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, kfft0d

!-----declare local variables
integer :: ift, ifd
integer :: istat


if( myid == 0 ) then
    read(iunit,iostat=istat) fft3x(1:kfft0d)
    ierro3 = abs(istat)
end if
!--- error trap
call gimax(ierro3)
if( ierro3 /= 0 ) return

call dbcast(fft3x,kfft0d,0)

do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glocal(ifd) = fft3x(ift)
end do
!--- store local charge density
call distlc( nfile, myid, nodes, rhomx, mshnod, glocal, mftdsp )


return
end




subroutine rdvhar( nfile,  &
& fname1, lstart, vhar, mshnod, ldouble_grid,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!-----------------------------------------------------------------------
!     load hartree potential
!-----------------------------------------------------------------------
use outfile
implicit none
integer, dimension(*) :: nfile
character(*) :: fname1
logical :: lstart
integer :: mshnod
real*8, dimension(mshnod) :: vhar
logical :: ldouble_grid
real*8  :: glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, kfft0d

!-----declare local variables
integer :: myid, nodes, nkd, myid_lr, nodes_lr, myid_pw, nodes_pw
logical :: ldouble_grid_
integer :: kfft0d_
real*8  :: ct, ct0, timecnt
integer :: iunit, ierror, istat


if( .not.lstart ) return


!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )

ct0 = timecnt()

!=======================================================================
ierror = 0
nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then
nkdif: if( nkd == 0 ) then
if( myid.eq.0 ) then

    call allocate_unit_number( iunit )

    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in rdvhar): ', trim(fname1)
        write(nfile(2),*) 'open file(in rdvhar): ', trim(fname1)
    end if

    open( iunit, file=fname1, status='old', action='read',         &
&         form='unformatted', iostat=ierror )

    if( ierror == 0 ) then
        read(iunit,iostat=istat) kfft0d_
        ierror = max( ierror, abs(istat) )

        if( ierror == 0 ) then
            !-----error trap
            if( kfft0d /= kfft0d_ ) ierror = 9
        end if

!        read(iunit,iostat=istat) ( glocal(m1), m1 = 1, ntotfd )
!        ierror = max( ierror, abs(istat) )
    end if

end if

!--- global maximum
call gimax(ierror)

if( ierror == 0 ) then
    call rd_local_data( nfile, myid, nodes, iunit, ierror,  &
& vhar, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
else
    vhar(1:mshnod) = 0.d0
end if

if( myid.eq.0 ) then

    if( ierror == 0 ) then
        read(iunit,iostat=istat) ldouble_grid_
        ierror = max( ierror, abs(istat) )

        if( ierror == 0 ) then
            !-----error trap
            if( .not.ldouble_grid .and.      ldouble_grid_ .or.  &
&                    ldouble_grid .and. .not.ldouble_grid_ ) ierror = 9
        end if

    end if
    if( ierror /= 0 ) then
        write(nfile(1),*) 'error in file:',  trim(fname1), ', err=',ierror
        write(nfile(2),*) 'error in file:',  trim(fname1), ', err=',ierror
    end if

end if
!--- global maximum
call gimax(ierror)
!if( ierror == 0 ) call dbcast(glocal,ntotfd,0)
end if nkdif

!-----set communicator
call get_worldun( myid, nodes )

call ibcast(ierror,1,0)

!if( ierror == 0 ) then
!if( ldouble_grid ) then
!    call rdvhar_dg( nfile, iunit )
!end if
!end if

!-----set communicator
call get_worldkd( myid, nodes, nkd )

if( nkd == 0 ) then
if( myid == 0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
end if
end if
!=======================================================================

!if( nkd == 0 ) then
!if( ierror == 0 ) then
!
!    !-----store local hartree potential
!    call distlc( nfile, myid, nodes,  &
!& vhar, mshnod, glocal, mftdsp )
!
!  else
!
!    vhar(1:mshnod) = 0.d0
!
!end if
!end if

!-----set communicator
call get_worldun( myid, nodes )
!-----internode synchronization
call gsync

call dbcast(vhar,mshnod,0)
!if( ldouble_grid ) call dbcastvhar_dg

end if npwif

if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync

    call dbcast(vhar,mshnod,0)
!    if( ldouble_grid ) call dbcastvhar_dg
end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync

    call dbcast(vhar,mshnod,0)
!    if( ldouble_grid ) call dbcastvhar_dg
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )


ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '   set hartree pot. : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) '   set hartree pot. : cpu-time :', ct - ct0
ct0 = ct


return
end




subroutine svwfns( nfile, iogpsz,  &
& fname1, cgjr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, eig, egvlud, occ, wegud,  &
& rhcr, lsreal8, lspin, nspnmx, gdcrsv, nspnod, lcgjsv, ncscale, wvratio )
!-----------------------------------------------------------------------
!     save wave functions : rhcr
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
integer :: iogpsz
character(*) :: fname1
dimension cgjr(*)
dimension rhcr(*)
dimension nbncnt(*), nbndsp(*)
real*8  :: eig(nband), occ(nband)
logical :: lsreal8
real*8  :: egvlud(*), wegud(*)
real*8  :: gdcrsv(*)
logical :: lspin, lcgjsv
integer :: ncscale
real*8  :: wvratio(*)

!-----declare local variables
character(80) :: fname
integer :: iunit, digit
integer :: nfiles, ii, isnd, nsd
real(8), allocatable :: abuf(:)


ct0 = timecnt()

!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then
if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    fname = trim(fname1)

    if( myid == 0 ) then
        if(loutfile(1)) write(nfile(1),*) 'open file(in svwfns): ', trim(fname)
        if(loutfile(2)) write(nfile(2),*) 'open file(in svwfns): ', trim(fname)
    end if

    call allocate_unit_number( iunit )

    open(iunit, file=fname, status='unknown', form='unformatted')
!    write(iunit) -1
    write(iunit) -ncscale
    write(iunit) nband
    write(iunit) nplw
    write(iunit) lsreal8
    write(iunit) nspnmx
    do nspin = 1, nspnmx
       if( lspin ) then
           call ldeig( nspin, eig, egvlud, nband )
           call ldocc( nspin, occ, wegud, nband )
       end if
       write(iunit) eig(1:nband)
       write(iunit) occ(1:nband)
    end do

    !-----for divide-and-conquer MD
    if( nfiles > 1 ) then
        allocate(abuf(nband*2))
        do ii=1, nfiles-1
        do nspin = 1, nspnmx
           call cdrecvs(myid_qm_un+ii,abuf,nband*2,myid_qm_un+ii,0)

           write(iunit) abuf(1:nband*2)
        end do
        end do
    end if

    if( ncscale == 2 ) then
        !-----noncollinear magnetism
        write(iunit) wvratio(1:nband*2)

        !-----for divide-and-conquer MD
        if( nfiles > 1 ) then
            do ii=1, nfiles-1
               call cdrecvs(myid_qm_un+ii,abuf,nband*2,myid_qm_un+ii,0)

               write(iunit) abuf(1:nband*2)
            end do
        end if
    end if

  else ioif

    !-----for divide-and-conquer MD
    allocate(abuf(nband*2))
    isnd=(myid_qm_un/iogpsz)*iogpsz
    do nspin = 1, nspnmx
       if( lspin ) then
           call ldeig( nspin, eig, egvlud, nband )
           call ldocc( nspin, occ, wegud, nband )
       end if
       abuf(1:nband) = eig(1:nband)
       abuf(nband+1:nband*2) = occ(1:nband)

       call cdsend(myid_qm_un,abuf,nband*2,isnd,0)
    end do

    if( ncscale == 2 ) then
        !-----noncollinear magnetism
        abuf(1:nband*2) =  wvratio(1:nband*2)
        call cdsend(myid_qm_un,abuf,nband*2,isnd,0)
    end if

  end if ioif

end if

do nspin = 1, nspnmx
   if( lspin ) then
       call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
   end if

   call svwfn22( nfile, myid, nodes, iogpsz, iunit,  &
& rhcr, nplwex*ncscale, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, cgjr )

end do
if( myid.eq.0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
end if
end if npwif

if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync
end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )
!---------Internode synchronization
!      call gsync

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '              svwfns :          :', ct- ct0
if(loutfile(2)) write(nfile(2),*) '              svwfns :          :', ct- ct0
ct0 = ct


return
end




subroutine svwfn22( nfile, myid, nodes, iogpsz, iunit,  &
& datasaved, ndat, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, work )
!-----------------------------------------------------------------------
!     save wave functions : datasaved
!
!    collect data to node 0, and save to file
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iogpsz
integer :: iunit
integer :: ndat, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: datasaved(ndat,nbnod), work(ndat,nbnod)

!-----declare local variables
integer :: idst, is1, is2, is3, n, i, nrc, nsd
real*8  :: eigmax, efact, efactr
integer :: status
integer :: nkd, myid_qm_un, nodes_qm_un
integer :: nfiles, ii, isnd, iogp


!-----set communicator
!-----to get myid_qm_un, nodes_qm_un
call get_worldqmun( myid_qm_un, nodes_qm_un )

!-----set communicator
call get_worldkd( myid, nodes, nkd )


if( .not.lsreal8 ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < ndat ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(ndat), stat=status )
        end if
    end if
    call ibcast(status,1,0) 
    if( status /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in svwfn2', status
            write(nfile(2),*) '*** memory allocation error in svwfn2', status
        end if
        return
    end if
end if


ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD

    do idst = 0, nodes - 1
       if( myid == 0 ) then
           !--- set No. of recieving data (bands)
           is1 = nbndsp(idst+1) + 1
           is2 = nbndsp(idst+1) + nbncnt(idst+1)
           is3  = is2 - is1 + 1
           if( idst == 0 ) then
               work(1:ndat,1:is3) = datasaved(1:ndat,1:is3)
             else
               nrc = is3*ndat
               if( nrc > 0 ) call cdrecv(100,work,nrc,0)
           end if

           if( lsreal8 ) then
               !--- save data in real*8
               do n = 1, is3
                  write(iunit) work(1:ndat,n)
               end do
           else
               !--- save data in integer*2
               do n = 1, is3
                  eigmax = 0.d0
                  do i = 1, ndat
                     eigmax = max( eigmax, abs(work(i,n)) )
                  end do
                  efact  = 32767.d0/eigmax
                  efactr = eigmax/32767.d0
                  saveigv(1:ndat) = efact*work(1:ndat,n)
                  write(iunit) efactr
                  write(iunit) saveigv(1:ndat)
               end do
           end if

       else if( myid == idst ) then
           nsd = nbnod*ndat
           if( nsd > 0 ) call cdsend(100,datasaved,nsd,0,0)
       end if

       call gsync
    end do

end if ioif


return
end




subroutine svwfn2( nfile, myid, nodes, iunit,  &
& datasaved, ndat, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, work )
!-----------------------------------------------------------------------
!     save wave functions : datasaved
!
!    collect data to node 0, and save to file
!-----------------------------------------------------------------------
use savedata
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ndat, nband, nbnod1, nbnod2, nbnod
integer :: nbncnt(nodes), nbndsp(nodes)
real*8  :: datasaved(ndat,nbnod), work(ndat,nbnod)

!-----declare local variables
integer :: idst, is1, is2, is3, n, i, nrc, nsd
real*8  :: eigmax, efact, efactr
integer :: status


if( .not.lsreal8 ) then
    if( myid == 0 ) then
        if( allocated(saveigv) ) then
            if( size(saveigv) < ndat ) deallocate( saveigv )
        end if
        if( .not.allocated(saveigv) ) then
            allocate( saveigv(ndat), stat=status )
        end if
    end if
    call ibcast(status,1,0) 
    if( status /= 0 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** memory allocation error in svwfn2', status
            write(nfile(2),*) '*** memory allocation error in svwfn2', status
        end if
        return
    end if
end if

do idst = 0, nodes - 1
   if( myid == 0 ) then
!--- set No. of recieving data (bands)
       is1 = nbndsp(idst+1) + 1
       is2 = nbndsp(idst+1) + nbncnt(idst+1)
       is3  = is2 - is1 + 1
       if( idst == 0 ) then
           work(1:ndat,1:is3) = datasaved(1:ndat,1:is3)
         else
           nrc = is3*ndat
           if( nrc > 0 ) call cdrecv(100,work,nrc,0)
       end if

       if( lsreal8 ) then
!--- save data in real*8
           do n = 1, is3
              write(iunit) work(1:ndat,n)
           end do
         else
!--- save data in integer*2
           do n = 1, is3
              eigmax = 0.d0
              do i = 1, ndat
                 eigmax = max( eigmax, abs(work(i,n)) )
              end do
              efact  = 32767.d0/eigmax
              efactr = eigmax/32767.d0
              saveigv(1:ndat) = efact*work(1:ndat,n)
              write(iunit) efactr
              write(iunit) saveigv(1:ndat)
           end do
       end if

   else if( myid == idst ) then
       nsd = nbnod*ndat
       if( nsd > 0 ) call cdsend(100,datasaved,nsd,0,0)
   end if

   call gsync
end do


return
end




subroutine svcdty( nfile,  &
& fname1, rhgr, nplw5ex, nplw5, nspnmx, lnoncollinear, mshnod,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!-----------------------------------------------------------------------
!    store electron density
!-----------------------------------------------------------------------
use ncmagne_variables
implicit none
integer :: nfile(*)
character(*) :: fname1
integer :: nplw5ex, nplw5, nspnmx
real*8  :: rhgr(nplw5ex,*)
logical :: lnoncollinear
integer :: mshnod
real*8  :: glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0d

!-----declare local variables
integer :: myid, nodes, myid_kd, nodes_kd, nkd
integer :: iunit, nversion
integer :: m1, nspin
real*8  :: ct0, ct, timecnt

!-----set communicator
call get_worldqm( myid, nodes )


ct0 = timecnt()

if( myid == 0 ) then
    call allocate_unit_number( iunit )
!=======================================================================
    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in svcdty): ', trim(fname1)
        write(nfile(2),*) 'open file(in svcdty): ', trim(fname1)
    end if

    open(iunit, file=fname1, status='unknown', form='unformatted')
    nversion = -3
    write(iunit) nversion
    write(iunit) nplw5
    write(iunit) nspnmx          ! <- nversion = -3
    write(iunit) lnoncollinear   ! <- nversion = -3
    do nspin = 1, nspnmx
       write(iunit) ( rhgr(m1,nspin), m1 = 1, nplw5ex )
    end do
!=======================================================================
end if


!---non-collinear magnetism
!if( lnoncollinear ) then
!    !-----set communicator
!    call get_worldkd( myid_kd, nodes_kd, nkd )
!
!    if( myid == 0 ) write(iunit) kfft0d
!
!    !---rhomx
!    call sv_local_data( nfile, myid, iunit,  &
!& rhomx, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!
!    !---rhomy
!    call sv_local_data( nfile, myid, iunit,  &
!& rhomy, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!
!    !---rhomz
!    call sv_local_data( nfile, myid, iunit,  &
!& rhomz, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!
!    !-----set communicator
!    call get_worldqm( myid, nodes )
!end if


if( myid == 0 ) then
!=======================================================================
    close(iunit)
!=======================================================================
    call deallocate_unit_number( iunit )
end if

!---------Internode synchronization
!      call gsync

ct = timecnt()
if( myid.eq.0 ) then
    write(nfile(1),*) '              svcdty :          :', ct- ct0
    write(nfile(2),*) '              svcdty :          :', ct- ct0
end if
ct0 = ct


!-----set communicator
call get_worldkd( myid_kd, nodes_kd, nkd )


return
end




subroutine sv_local_data( nfile, myid, iunit,  &
& rhomx, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!-----------------------------------------------------------------------
!    store rhomx by gathering it to node 0
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid
integer :: iunit
integer :: mshnod
real*8  :: rhomx(mshnod), glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, kfft0d

!-----declare local variables
integer :: ift, ifd


call dgatherv( rhomx, mshnod, glocal, mftnod, mftdsp, 0 )
if( myid == 0 ) then
    fft3x(1:kfft0d) = 0.d0
    do ifd = 1, ntotfd
       ift = mfd2ft(ifd)
       fft3x(ift) = glocal(ifd)
    end do
    write(iunit) fft3x(1:kfft0d)
end if


return
end




subroutine svvhar( nfile,  &
& fname1, vhar, mshnod, ldouble_grid,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!-----------------------------------------------------------------------
!     save hartree potential
!-----------------------------------------------------------------------
implicit none
integer, dimension(*) :: nfile
character(*) :: fname1
integer :: mshnod
real*8, dimension(mshnod) :: vhar
logical :: ldouble_grid
real*8  :: glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, kfft0d

!-----declare local variables
integer :: myid, nodes, nkd, myid_lr, nodes_lr, myid_pw, nodes_pw
real*8  :: ct, ct0, timecnt
integer :: iunit


!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )


nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then
nkdif: if( nkd == 0 ) then

ct0 = timecnt()

!--- unify hartree potential
!call unifylc( nfile, myid, nodes,  &
!& vhar, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )


!=======================================================================
if( myid == 0 ) then
    call allocate_unit_number( iunit )

    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in svvhar): ', trim(fname1)
        write(nfile(2),*) 'open file(in svvhar): ', trim(fname1)
    end if

    open(iunit, file=fname1, status='unknown', form='unformatted')
!    write(iunit) ntotfd
!    write(iunit) ( glocal(m1), m1 = 1, ntotfd )
    write(iunit) kfft0d
end if

call sv_local_data( nfile, myid, iunit,  &
& vhar, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )

if( myid == 0 ) then
    write(iunit) ldouble_grid
end if

!if( ldouble_grid ) then
!    call svvhar_dg( nfile, myid, nodes, iunit )
!end if

if( myid == 0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
end if
!=======================================================================


ct = timecnt()
if( myid == 0 ) then
    write(nfile(1),*) '              svvhar :          :', ct- ct0
    write(nfile(2),*) '              svvhar :          :', ct- ct0
end if
ct0 = ct

end if nkdif
end if npwif

if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync
end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




subroutine svacon( nfile, iogpsz,  &
& fname1, prevr, ifmd, nstepCG, nstepMD,  &
& ratm, frc, ntype, nhk1, nhk2, lclust, natom, hcell )
!-----------------------------------------------------------------------
!    store atomic configuration
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension   nfile(*)
character(*) :: fname1
dimension   prevr(3,natom,3)
dimension   ratm(3,natom)
dimension   frc(3,natom)
dimension   nhk1(ntype), nhk2(ntype)
real*8, dimension(3,3) :: hcell
logical     lclust

!-----declare local variables
character(80) :: fname
integer :: myid, nodes, myid_qm_un, nodes_qm_un, digit
integer :: iunit
integer :: nfiles, ii, isnd, nsd


!-----set communicator
call get_worldqm( myid, nodes )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


!-----set version, must be negative
nversion = -2

ct0 = timecnt()

if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    fname = trim(fname1)

    call allocate_unit_number( iunit )
!=======================================================================
    if( myid == 0 ) then
        if(loutfile(1)) write(nfile(1),*) 'open file(in svacon): ', trim(fname)
        if(loutfile(2)) write(nfile(2),*) 'open file(in svacon): ', trim(fname)
    end if

    open(iunit, file=fname, status='unknown', form='unformatted')
    write(iunit) nversion
    write(iunit) ntype
    write(iunit) ( nhk1(it), nhk2(it), it = 1, ntype )
    write(iunit) nstepCG, nstepMD
    write(iunit) lclust
    write(iunit) ((ratm(ix,i), ix = 1, 3), i = 1, nhk2(ntype))
    write(iunit) ifmd
    write(iunit) ((frc(ix,i), ix = 1, 3), i = 1, nhk2(ntype))
    write(iunit) (((prevr(ix,i,j), ix = 1, 3), i = 1, nhk2(ntype)), j = 1, 3)
    write(iunit) ( ( hcell(ix,j), ix = 1, 3), j = 1, 3 )
    close(iunit)
!=======================================================================
    call deallocate_unit_number( iunit )

  end if ioif
end if

ct = timecnt()
if( myid.eq.0 ) then
    if(loutfile(1)) write(nfile(1),*) '              svacon :          :', ct- ct0
    if(loutfile(2)) write(nfile(2),*) '              svacon :          :', ct- ct0
end if
ct0 = ct


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




subroutine rdhcell( nfile, hcell, lvshape )
!-----------------------------------------------------------------------
!     read hcell
!-----------------------------------------------------------------------
implicit none

integer, dimension(*) :: nfile
real*8,  dimension(3,3) :: hcell
logical :: lvshape

!-----declare local variables
integer :: iunit
integer :: i, ix, ierror, istat
integer :: myid, nodes, nkd


!-----set communicator
call get_worldqm( myid, nodes )


if( myid == 0 ) then

    call allocate_unit_number( iunit )

    call open_cellfiles( nfile, myid, nodes, iunit, ierror )

    istat = 0
    if( lvshape ) then
        !-----if (NPT) dynamics, read hcell
        read(iunit,'(3e23.15)',iostat=istat)  &
&                  ( ( hcell(ix,i), ix = 1, 3 ), i = 1, 3 )
    end if

    close(iunit)

    call deallocate_unit_number( iunit )

    if( istat /= 0 .or. .not.lvshape ) then
        do i = 1, 3
        do ix = 1, 3
           hcell(ix,i) = 0.d0
        end do
        end do
    end if

end if

call dbcast(hcell,9,0)


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




subroutine svhcell( nfile, hcell, lvshape )
!-----------------------------------------------------------------------
!     save hcell
!-----------------------------------------------------------------------
implicit none

integer, dimension(*) :: nfile
real*8,  dimension(3,3) :: hcell
logical :: lvshape

!-----declare local variables
integer :: iunit
integer :: i, ix, ierror
integer :: myid, nodes, nkd


!-----set communicator
call get_worldqm( myid, nodes )


if( myid == 0 ) then

    call allocate_unit_number( iunit )

    call open_cellfiles( nfile, myid, nodes, iunit, ierror )

    if( lvshape ) then
        !-----if (NPT) dynamics, write hcell
        write(iunit,'(3e23.15)')  &
&                  ( ( hcell(ix,i), ix = 1, 3 ), i = 1, 3 )

      else

        !-----if .not.(NPT) dynamics, write zeros
        write(iunit,'(3e23.15)') ( ( 0.d0, ix = 1, 3 ), i = 1, 3 )

    end if

    close(iunit)

    call deallocate_unit_number( iunit )

end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end
