



subroutine open_md_logfile( nfile, myid, nodes )
!-----------------------------------------------------------------------
!    Allocate I/0 files for log & scratch files
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: ierror = 0
character(50) :: dname = 'data/'


!loutfile(1) = .true.
loutfile(1) = myid == 0
loutfile(2) = myid == 0

!-----Allocate I/0 file
if( myid == 0 ) then

    open( nfile(2), file=dname(1:len_trim(dname)) // 'md_log',  &
&     status='unknown', action='write', iostat=ierror )

!  else
!
!    open( nfile(2), status='scratch', action='write',  &
!&     iostat=ierror )

end if

!----- error trap
ierror = abs(ierror)
call gimax(ierror)
if( ierror /= 0 ) call fstop( nfile, myid, nodes, 'cannot open file md_log' )


return
end subroutine




subroutine open_qm_logfile( nfile, myid, nodes,  &
& myid_qm, nodes_qm, myid_qm_un, nodes_qm_un, myid_kd, nodes_kd, nkd,  &
& myid_pw, nodes_pw, myid_lr, nodes_lr, ierror )
!-----------------------------------------------------------------------
!    Allocate I/0 files for log & scratch files
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: myid_qm, nodes_qm
integer :: myid_qm_un, nodes_qm_un
integer :: myid_kd, nodes_kd, nkd
integer :: myid_pw, nodes_pw
integer :: myid_lr, nodes_lr
integer :: ierror

!-----declare local variables
integer :: npk, digit
character(50) :: dname = 'data/'
character(1), dimension(0:9) :: num =  &
& (/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /)


npk = nodes_qm/nodes_kd/nodes_pw/nodes_lr

dname = trim(dname) // 'qm_log'
if( nodes_qm_un > 1 ) then
    digit = log10(dble(nodes_qm_un)+1.d-06) + 1.d0
    call get_dumpfname( dname, myid_qm_un, digit )
end if
if( npk > 1 ) then
    if( nodes_qm_un > 1 ) dname = trim(dname) // '_'
    digit = log10(dble(npk)+1.d-06) + 1.d0
    call get_dumpfname( dname, nkd, digit )
end if


loutfile(1) = nodes_qm_un == 1 .and. myid_qm == 0  &
&        .or. nodes_qm_un >  1 .and. myid    == 0
!if( myid_qm == 0 ) then

    nfile(1) = 06

!  else
!
!    !-----Allocate I/0 file
!    open( nfile(1), status='scratch', action='write',  &
!&     iostat=ierror )
!
!end if


loutfile(2) = myid_kd == 0 .and. myid_pw == 0 .and. myid_lr == 0 .and. myid_qm_un == 0
!-----Allocate I/0 file
if( loutfile(2) ) then

    open( nfile(2), file=trim(dname),  &
&         status='unknown', action='write', iostat=ierror )

!  else
!
!    open( nfile(2), status='scratch', action='write',  &
!&     iostat=ierror )

end if

!!----- error trap
!ierror = abs(ierror)
!call gimax(ierror)
!if( ierror /= 0 ) call fstop( nfile, myid, nodes,  &
!&                 'cannot open file qm_log' )


return
end subroutine




subroutine open_files( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!    Allocate I/0 files for QM nodes
!-----------------------------------------------------------------------
use outfile
use param, only:  &
& ifmd, lmulken, lspdatom, lpmchgat, lmdecmpdos,  &
& lsphexp, leda, lwannier, lacconduct, loutuni, lstress,  &
& lintchg, ltddft, ltddft_fssh, lfssh_gsscf, lrtddft, lnoncollinear, lefield
use pwlda_variables, only:  &
& fname_ion, fname_eig, fname_cds, fname_hrt, fname_pcds, fname_peig,  &
& fname_eigk, fname_peigk, fname_tddft,  &
& dname_wann, dname_eig, dname_cds, dname_pot, dname_potav
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
character(50) :: dname = 'data/'
character(50) :: fname
integer :: i


!-----set data-file names for the next run
fname_ion = dname(1:len_trim(dname)) // 'QM_ion'
fname_eig = dname(1:len_trim(dname)) // 'QM_eig'
fname_eigk= dname(1:len_trim(dname)) // 'QM_eigk'
fname_tddft = dname(1:len_trim(dname)) // 'QM_tddft'
fname_cds = dname(1:len_trim(dname)) // 'QM_cds'
fname_hrt = dname(1:len_trim(dname)) // 'QM_hrt'
fname_pcds= dname(1:len_trim(dname)) // 'QM_pcds'
fname_peig= dname(1:len_trim(dname)) // 'QM_peig'
fname_peigk= dname(1:len_trim(dname)) // 'QM_peigk'

dname_wann = dname(1:len_trim(dname)) // 'qm_wannier'
dname_eig  = dname(1:len_trim(dname)) // 'qm_eigv.d'
dname_cds  = dname(1:len_trim(dname)) // 'qm_cds.d'
dname_pot  = dname(1:len_trim(dname)) // 'qm_pot.d'
dname_potav  = dname(1:len_trim(dname)) // 'qm_potav.d'


ierror = 0
myidif: if( loutfile(1) ) then
!-----------------------------------------------------------------------

    fname = dname(1:len_trim(dname)) // 'qm_eig.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(3), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_ion.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(4), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_eng.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(5), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_zan.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(29), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_box.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(6), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_cel.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(27), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_frc.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(8), ierror )

    fname = dname(1:len_trim(dname)) // 'qm_fer.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(16), ierror )

    if( lintchg ) then
        fname = dname(1:len_trim(dname)) // 'qm_chg.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(9), ierror )
    end if

    if( lmulken ) then
       fname = dname(1:len_trim(dname)) // 'qm_mul.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(10), ierror )

       fname = dname(1:len_trim(dname)) // 'qm_pds.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(11), ierror )

       if( lspdatom ) then
           fname = dname(1:len_trim(dname)) // 'qm_pda.d'
           if( ierror == 0 )  &
&           call openf( nfile, myid, nodes, fname, nfile(18), ierror )
       end if

       if( lmdecmpdos ) then
           fname = dname(1:len_trim(dname)) // 'qm_pdm.d'
           if( ierror == 0 )  &
&           call openf( nfile, myid, nodes, fname, nfile(34), ierror )
       end if

       if( lpmchgat ) then
           fname = dname(1:len_trim(dname)) // 'qm_mum.d'
           if( ierror == 0 )  &
&           call openf( nfile, myid, nodes, fname, nfile(35), ierror )
       end if

       fname = dname(1:len_trim(dname)) // 'qm_ovp.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(12), ierror )

       if( lnoncollinear ) then
           fname = dname(1:len_trim(dname)) // 'qm_mag_mul.d'
           if( ierror == 0 )  &
&           call openf( nfile, myid, nodes, fname, nfile(32), ierror )
       end if
    end if

    if( lsphexp ) then
       fname = dname(1:len_trim(dname)) // 'qm_sph.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(7), ierror )
    end if

    if( leda ) then
       fname = dname(1:len_trim(dname)) // 'qm_eda.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(19), ierror )
    end if

    if( lwannier ) then
        fname = dname(1:len_trim(dname)) // 'qm_wan.d'
        if( ierror == 0 ) &
&       call openf( nfile, myid, nodes, fname, nfile(13), ierror )

        if( loutuni ) then
            fname = dname(1:len_trim(dname)) // 'qm_uni.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(15), ierror )
        end if

        fname = dname(1:len_trim(dname)) // 'qm_wsp.d'
        if( ierror == 0 ) &
&       call openf( nfile, myid, nodes, fname, nfile(17), ierror )
    end if

   if( lacconduct ) then
       fname = dname(1:len_trim(dname)) // 'qm_mmt.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(26), ierror )
    end if

    if( lstress ) then
       fname = dname(1:len_trim(dname)) // 'qm_str.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(14), ierror )
    end if

    if( ltddft .or. ltddft_fssh .and. lfssh_gsscf ) then
       fname = dname(1:len_trim(dname)) // 'qm_td_eig.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(20), ierror )
    end if

    if( lrtddft ) then
        fname = dname(1:len_trim(dname)) // 'qm_td_exe.d'
        if( ierror == 0 )  &
&           call openf( nfile, myid, nodes, fname, nfile(28), ierror )
    end if

!       if( ltddft_fssh .and. lfssh_gsscf ) then
!           !---exchange unit number
!           i         = nfile(3)
!           nfile(3)  = nfile(20)
!           nfile(20) = i
!       end if
!    end if

    if( ltddft ) then
       fname = dname(1:len_trim(dname)) // 'qm_td_off.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(21), ierror )

       fname = dname(1:len_trim(dname)) // 'qm_td_cef.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(22), ierror )

       fname = dname(1:len_trim(dname)) // 'qm_td_nrm.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(23), ierror )
    end if

    if( lrtddft ) then
       fname = dname(1:len_trim(dname)) // 'qm_lrtddft.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(24), ierror )

       fname = dname(1:len_trim(dname)) // 'qm_lrtddft_osc.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(25), ierror )
    end if

    if( lnoncollinear ) then
       fname = dname(1:len_trim(dname)) // 'qm_mag.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(30), ierror )

       fname = dname(1:len_trim(dname)) // 'qm_mag_atm.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(31), ierror )
    end if

    if( lefield ) then
       fname = dname(1:len_trim(dname)) // 'qm_pol.d'
       if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(33), ierror )
    end if
!-----------------------------------------------------------------------
!  else
!
!    do i = 3, 30
!       !-----Allocate I/0 file
!       open( nfile(i), status='scratch', action='write',  &
!&        iostat=ierror )
!    end do
!
!end if
end if myidif


!--- global maximum
call gimax(ierror)
!      !----- error trap
!      if( ierror /= 0 )  &
!     &    call fstop( nfile, myid, nodes, ' ' )


return
end subroutine




subroutine openf( nfile, myid, nodes, fname, iunit, ierror )
!-----------------------------------------------------------------------
!    Allocate an I/0 file
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
character(50) :: fname
integer :: iunit
integer :: ierror


if(loutfile(1)) write(nfile(1),*) 'open file: ', fname(1:len_trim(fname))
if(loutfile(2)) write(nfile(2),*) 'open file: ', fname(1:len_trim(fname))

open( iunit, file=fname, status='unknown',  &
&            action='readwrite', iostat=ierror )

!----- error trap
if( ierror /= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'myid=',myid,' cannot open file :'//fname
    if(loutfile(2)) write(nfile(2),*) 'myid=',myid,' cannot open file :'//fname
end if


return
end subroutine




subroutine open_a_file( nfile, myid, nodes, fname, iunit, ierror, lwrite )
!-----------------------------------------------------------------------
!     Allocate an I/0 file
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
character(50) :: fname
integer :: iunit, ierror
logical :: lwrite


if(loutfile(1)) write(nfile(1),*) 'open file: ', trim(fname)
if(loutfile(2)) write(nfile(2),*) 'open file: ', trim(fname)

if( lwrite ) then
    open( iunit, file=fname, status='unknown', form='unformatted',  &
&                action='write', iostat=ierror )
else
    open( iunit, file=fname, status='old', form='unformatted',  &
&                action='read', iostat=ierror )
end if

!----- error trap
if( ierror /= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'myid=',myid,' cannot open file :'//trim(fname)
    if(loutfile(2)) write(nfile(2),*) 'myid=',myid,' cannot open file :'//trim(fname)
end if


return
end subroutine




subroutine delete_a_file( nfile, myid, nodes, fname, ierror )
!-----------------------------------------------------------------------
!     Delete an unnecessary I/0 file
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
character(50) :: fname
integer :: ierror

!-----declare local variables
logical :: lexist
integer :: iunit


inquire( file=fname, exist=lexist )
if( lexist ) then
    call allocate_unit_number( iunit )
    call open_a_file( nfile, myid, nodes, fname, iunit, ierror, .true. )
    !------delete an unnecessary file
    close(iunit,status='delete')
    call deallocate_unit_number( iunit )
end if


return
end subroutine




subroutine open_cellfiles( nfile, myid, nodes, iunit, ierror )
!-----------------------------------------------------------------------
!    Allocate I/0 files for hcell in QM nodes
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit, ierror

!-----declare local variables
character(50) :: dname = 'data/'
character(50) :: fname


fname = dname(1:len_trim(dname)) // 'QM_cell'

if(loutfile(1)) write(nfile(1),*) 'open file: ', fname(1:len_trim(fname))
if(loutfile(2)) write(nfile(2),*) 'open file: ', fname(1:len_trim(fname))

open( iunit, file=fname, status='unknown',  &
&            iostat=ierror )

!----- error trap
if( ierror /= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'myid=',myid,' cannot open file :'//fname
    if(loutfile(2)) write(nfile(2),*) 'myid=',myid,' cannot open file :'//fname
end if


return
end subroutine




subroutine open_tddftfsshfiles( nfile, myid, nodes,  &
& myid_qm_un, nodes_qm_un, iunit, ierror )
!-----------------------------------------------------------------------
!    Allocate I/0 files for TDDFT in QM nodes
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: myid_qm_un, nodes_qm_un
integer :: iunit, ierror

!-----declare local variables
character(50) :: dname = 'data/'
character(50) :: fname
integer :: digit


fname = trim(dname) // 'QM_tddftfssh'

if(loutfile(1)) write(nfile(1),*) 'open file: ', trim(fname)
if(loutfile(2)) write(nfile(2),*) 'open file: ', trim(fname)

open( iunit, file=fname, status='unknown', form='unformatted', iostat=ierror )

!----- error trap
if( ierror /= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'myid=',myid,' cannot open file :'//trim(fname)
    if(loutfile(2)) write(nfile(2),*) 'myid=',myid,' cannot open file :'//trim(fname)
end if


return
end subroutine




subroutine title_in_files( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     title in files for QM nodes
!-----------------------------------------------------------------------
use outfile
use param, only:  &
& ifmd, lpaw, lmulken, lspdatom, lpmchgat, lmdecmpdos,  &
& leda, lwannier, lacconduct, loutuni, lstress, lvdw, &
& ltddft, ltddft_fssh, lfssh_gsscf, ldftd, lrtddft, lplusU, jhybrid, lnoncollinear, &
& lefield, lsawtooth, lconstraintD, lrela_full, lsoc_full
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!---Declare local variables
character(300) :: message
integer :: i


myidif: if( loutfile(1) ) then
!-----------------------------------------------------------------------
!--- eigenvalues ---
    if( .not.lnoncollinear ) then
        write(nfile(3),2100) 'Eigenvalues                   '
    else
        write(nfile(3),'(a)') '#       Eigenvalues  occ.   M_x    M_y    M_z'
    end if

!--- Fermi energy ---
    write(nfile(16),2100) 'Fermi energy                  '

!--- supercell (FFT cell) vectors ---
    write(nfile(6),'(a)') '#  supercell (FFT cell) vectors (lengths & angles)'
    write(nfile(6),2110) '        L_1         ',  &
&      '  L_2           L_3 ','         angle(2-3) ',  &
&      'angle(3-1) angle(1-2',')                   '

    write(nfile(27),'(a)') '#  supercell (FFT cell) vectors'
    write(nfile(27),'(a)') '#        L_(1:3,1:3)'

!--- atomic configuration ---
    write(nfile(4),2100) 'Atomic scaled coordinates     '
    write(nfile(8),2100) 'Atomic forces in [a.u.]       '

!--- total potential energy and energy parts ---
write(nfile(5),2100) 'Total potential energy and ene',  &
                   'rgy parts in [Ryd.] units     '
message ='#              Total(HF)         Total(KS)  &
&       Kinetic       External      Hartree       Exchange'
if( lvdw ) then
  message = message(1:len_trim(message)) &
            //'      Ec(LDA)       Ec(nl)        Entropy'
else
  message = message(1:len_trim(message)) &
            //'      Correlation   ------        Entropy'
end if
if( lpaw ) then
  message = message(1:len_trim(message)) //'       Onsite E.'
else
  message = message(1:len_trim(message)) //'       ---------'
end if
if( jhybrid /= 0 ) then
  message = message(1:len_trim(message)) //'     SR-HF E.'
else
  message = message(1:len_trim(message)) //'     --------'
end if
if( lplusU ) then
  message = message(1:len_trim(message)) //'      DFT+U E.'
else
  message = message(1:len_trim(message)) //'      --------'
end if
if( lefield ) then
  message = message(1:len_trim(message)) //'     Efield on el.'
else
  message = message(1:len_trim(message)) //'     -------------'
end if
if( lrela_full .or. lsoc_full ) then
  message = message(1:len_trim(message)) //'  SO E.'
else
  message = message(1:len_trim(message)) //'  -----'
end if
  message = message(1:len_trim(message)) //'         Ewald E.'
if( ldftd ) then
  message = message(1:len_trim(message)) //'      DFT-D'
else
  message = message(1:len_trim(message)) //'      -----'
end if
if( lefield ) then
  message = message(1:len_trim(message)) //'        Efield on ion'
else
  message = message(1:len_trim(message)) //'        -------------'
end if
write(nfile(5),'(300a1)') (message(i:i), i = 1, len_trim(message))

!--- Residuals
    write(nfile(29),'(a)') '# Difference of tot E/el. E(HF)-E(KS) & maximum & average residuals'
    write(nfile(29),'(a,a)') '#              difene      difene2     zansa1      zansa2',  &
&                          '      bfzansa1    bfzansa2'

!--- Momentum matrix to calculate optical conductivity
    if( lacconduct ) then
        message ='# Momentum matrix <m|nabra|n> where n/m are occupied/unoccupied states'
        write(nfile(26),'(300a1)') (message(i:i), i = 1, len_trim(message))
    end if

!--- stress ---
    if( lstress ) then
        write(nfile(14),2100) 'Stress in [GPa]               '
        write(nfile(14),2110) '       Pxx          ',  &
&      '   Pyy             P','zz             Pyz  ',  &
&      '           Pzx      ','       Pxy          '
    endif

!--- Mulliken analysis
    if( lmulken ) then
        write(nfile(10),2100) 'Mulliken analysis : s, p, & d ',  &
&                             'population for each atom      '
        write(nfile(11),2100) 'Mulliken analysis : s, p, & d ',  &
&                             'contribution to each band     '
        if( lspdatom ) then
            write(nfile(18),2100) 'Mulliken analysis : s, p, & d ',  &
&                                 'contributions to each band ass',  &
&                                 'ociated with each atom        '
        end if
        if( lmdecmpdos ) then
            write(nfile(34),2100) 'Mulliken analysis : s, pz, px,',  &
&                                 ' py, dz2, dzx, dx2-y2, dyz & d',  &
&                                 'xy contributions to each band '    
        end if
        if( lpmchgat ) then
            write(nfile(35),2100) 'Mulliken analysis : s, pz, px,',  &
&                                 ' py, dz2, dzx, dx2-y2, dyz & d',  &
&                                 'xy population for each atom   '    
        end if
        write(nfile(12),2100) 'Mulliken analysis : overlap po',  &
&                             'pulation between atoms        '
        if( lnoncollinear ) then
            write(nfile(32),2100) 'Mulliken analysis : s, p, & d ',  &
&                                 'contributions to magnetic mome',  &
&                                 'nts of each atom              '
        end if
    end if

!--- Energy Density Analysis (EDA)
    if( leda ) then
        write(nfile(19),2100) 'Energy Density Analysis       '
        write(nfile(19),2110) '         Total E.   ',  &
&      'Kinetic E.  local E.','(1)   Hartree E.    ',  &
&      '  E_xc.      Onsite ','E.  local E.(2)   Co',  &
&      'ulomb E.            '
    end if

!--- Wannier function calculation
    if( lwannier ) then
        write(nfile(13),2100) 'Wannier function  : center of ',  &
&                             'Wannier functions             '
        if( loutuni ) then
        write(nfile(15),2100) 'Wannier function  : unitary ma',  &
&                             'trix                          '
        end if
        write(nfile(17),2100) 'Wannier function  : spread of ',  &
&                             'Wannier functions             '
    end if

    if( ltddft .or. ltddft_fssh .and. lfssh_gsscf ) then
        if( ltddft_fssh .and. lfssh_gsscf ) then
            write(nfile(20),2100) 'Eigenvalues of GS & occupation', &
&                                 's of Excited States           '
        else
            write(nfile(20),2100) 'Diagonal elements of Hamiltoni', &
&                                 'an matrix by wave packets     '
        end if
    end if
    if( lrtddft ) then
        if( ltddft_fssh ) then !.and. lfssh_gsscf ) then
            write(nfile(28),'(a)') '# Excitation energies by LR-TDDFT & occupations'
        else
            write(nfile(28),'(a)') '# Excitation energies by LR-TDDFT'
        end if
    end if
    if( ltddft ) then
        write(nfile(21),2100) 'Maximum of off diagonal elemen', &
                              'ts of Hamiltonian matrix by wa', &
&                             've packets                    '
        write(nfile(22),2100) 'Coefficients of eigenvector ex', &
                              'pansion of wave packets       '
        write(nfile(23),2100) 'Average, max, and min norms of', &
                              ' occupied expanded wave packet', &
                              's                             '
    end if
    if( lrtddft ) then
        write(nfile(25),'(a)') '# spectroscopic oscillator strength'
        write(nfile(25),'(a)')  &
& '# Excitaton E.(Ry)      f_x        f_y        f_z        ave.'
    end if

    if( lnoncollinear ) then
        write(nfile(30),'(a)') '# total magnetic moments (mx, my, mz) and |m|'
        write(nfile(31),'(a)') '# total magnetic moment (mx, my, mz) and |m|'
    end if

    if( lefield ) then
        write(nfile(33),'(a)') '# Polarizations for the whole system'
        if( lsawtooth ) then
            write(nfile(33),'(a)')  &
& '#             electronic polarization                           ionic polarization'
        else
            if( .not.lconstraintD ) then
                write(nfile(33),'(a)')  &
& '#             electronic polarization                           ionic polarization'//  &
& '                               external field D=E+4*pi*P in [Hartree]'
            else
                write(nfile(33),'(a)')  &
& '#             electronic polarization                           ionic polarization'//  &
& '                               electric field E=D-4*pi*P in [Hartree]'
            end if
        end if
    end if
!-----------------------------------------------------------------------
end if myidif

2100 format('#',2x,3a30)
2110 format('#',10a20)
2111 format('#',10x,10a20)
2120 format('#',a20,i2)


return
end subroutine




subroutine remd_open_files( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!    Allocate I/0 files for MD nodes
!
!     nfile(17) is reserved for the calculation of density fluctuation.
!-----------------------------------------------------------------------
use remd_param, only:  &
& iogpsz, ifmd, lmdout, locoor, lovelo, loforc,  &
& lmdout_qm, locoor_qm, lovelo_qm, loforc_qm,  &
& lmdout_cl, locoor_cl, lovelo_cl, loforc_cl,  &
& ncbonds, lstress, nskip_stress_out,  &
& latomic, lheat, ltotmom, lvmd
use remd_variables, only: fname_toc, fname_mts
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
character(50) :: dname = 'data/'
character(50) :: fname
integer :: digit


!-----set data-file names for the next run
fname_toc = dname(1:len_trim(dname)) // 'MD_Tocontinue'
fname_mts = dname(1:len_trim(dname)) // 'MD_mts'


ierror = 0
!-----------------------------------------------------------------------
myidif: if( myid == 0 ) then

    !-----open file for energies
    fname = dname(1:len_trim(dname)) // 'md_eng.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(3), ierror )

    !-----open file for MD cell vectors
    fname = dname(1:len_trim(dname)) // 'md_box.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(4), ierror )

    !-----open file for MD cell vectors
    fname = dname(1:len_trim(dname)) // 'md_cel.d'
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(9), ierror )

!    if( lstress ) then
    if( lstress .and. nskip_stress_out > 0 ) then
        !-----open file for internal sterss tensor
        fname = dname(1:len_trim(dname)) // 'md_str.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(16), ierror )

        !-----open file for diagonarized internal sterss tensor
        fname = dname(1:len_trim(dname)) // 'md_str_diag.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(22), ierror )
    end if

    if( ifmd >= 1 .and. ncbonds > 0 ) then
        fname = dname(1:len_trim(dname)) // 'md_lag.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(18), ierror )
    end if

    if( ifmd == 10 ) then
        fname = dname(1:len_trim(dname)) // 'md_hug.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(23), ierror )
    end if

    if( lvmd ) then
        fname = dname(1:len_trim(dname)) // 'md_ref.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(24), ierror )
    end if

!    if( lheat ) then
!!        fname = dname(1:len_trim(dname)) // 'md_heat.d'
!        fname = dname(1:len_trim(dname)) // 'md_heat'
!        !---Note that a binary file is opened.
!        if( ierror == 0 )  &
!!&       call openf( nfile, myid, nodes, fname, nfile(25), ierror )
!&       call open_a_binary_file( nfile, myid, nodes,  &
!& fname, nfile(25), ierror, .true. )
!    end if

    if( ltotmom ) then
        fname = dname(1:len_trim(dname)) // 'md_mom.d'
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(26), ierror )
    end if

end if myidif


digit = log(dble(nodes)+0.5d0)/log(10.d0) + 1.d0
if( lmdout ) then

ioif: if( mod(myid,iogpsz) == 0 ) then

    !-----open file for atomic species
    fname = dname(1:len_trim(dname)) // 'md_spc.d'
    if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(5), ierror )

    !-----open file for atomic scaled coordinates
    if( locoor ) then
        fname = dname(1:len_trim(dname)) // 'md_ion.d'
        if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(6), ierror )
    end if

    if( ifmd >= 2 ) then
        !-----open file for atomic scaled velocities
        if( lovelo ) then
         fname = dname(1:len_trim(dname)) // 'md_vel.d'
         if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
         if( ierror == 0 )  &
&        call openf( nfile, myid, nodes, fname, nfile(7), ierror )
        end if
    end if

    !-----open file for atomic scaled forces
    if( loforc ) then
        fname = dname(1:len_trim(dname)) // 'md_frc.d'
        if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(8), ierror )
    end if

end if ioif

end if

!-----open file for atomic energies
if( latomic ) then

ioif2: if( mod(myid,iogpsz) == 0 ) then

    fname = dname(1:len_trim(dname)) // 'md_atme.d'
    if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
    if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(19), ierror )

    if( lstress ) then
        fname = dname(1:len_trim(dname)) // 'md_atmstrp.d'
        if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(20), ierror )

        fname = dname(1:len_trim(dname)) // 'md_atmstrk.d'
        if( nodes > 1 ) call get_dumpfname( fname, myid, digit )
        if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(21), ierror )
    end if

end if ioif2

end if


myidif2: if( myid == 0 ) then

    !-----open file for QM/MD run : for QM atoms
    if( lmdout_qm ) then

        !-----open file for scaled coordinates of QM atoms
        if( locoor_qm ) then
            fname = dname(1:len_trim(dname)) // 'md_qm_ion.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(10), ierror )
        end if

        if( ifmd >= 2 ) then
            !-----open file for scaled velocities of QM atoms
            if( lovelo_qm ) then
                fname = dname(1:len_trim(dname)) // 'md_qm_vel.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(11), ierror )
            end if
        end if

        !-----open file for scaled forces of QM atoms
        if( loforc_qm ) then
            fname = dname(1:len_trim(dname)) // 'md_qm_frc.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(12), ierror )
        end if

    end if

    !-----open file for QM/MD run : for MD-cluster atoms
    if( lmdout_cl ) then

        !-----open file for scaled coordinates of MD-cluster atoms
        if( locoor_cl ) then
            fname = dname(1:len_trim(dname)) // 'md_cl_ion.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(13), ierror )
        end if

        if( ifmd >= 2 ) then
            !-----open file for scaled velocities of MD-cluster atoms
            if( lovelo_cl ) then
                fname = dname(1:len_trim(dname)) // 'md_cl_vel.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(14), ierror )
            end if
        end if

        !-----open file for scaled forces of MD-cluster atoms
        if( loforc_cl ) then
            fname = dname(1:len_trim(dname)) // 'md_cl_frc.d'
            if( ierror == 0 )  &
&       call openf( nfile, myid, nodes, fname, nfile(15), ierror )
        end if

    end if

end if myidif2
!-----------------------------------------------------------------------


!--- global maximum
call gimax(ierror)


return
end subroutine




subroutine remd_close_files( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!    Close I/0 files for MD nodes
!-----------------------------------------------------------------------
use remd_param, only: iogpsz, lheat
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


myidif: if( myid == 0 ) then

!    if( lheat ) call close_a_binary_file( nfile, myid, nodes, nfile(25), ierror )

end if myidif


return
end subroutine




subroutine remd_title_in_files( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     title in files for MD nodes
!-----------------------------------------------------------------------
use remd_param, only:  &
& iogpsz, ifmd, ioptmze_cell, ntmd, lmdout, locoor, lovelo, loforc,  &
& lmdout_qm, locoor_qm, lovelo_qm, loforc_qm,  &
& lmdout_cl, locoor_cl, lovelo_cl, loforc_cl,  &
& lQMMDrun, ncbonds, lstress, latomic, ltotmom,  &
& lvmd, dlambda_vmd, treq, h, ntmd, lheat, dtmd, nskip_heat, volume
use remd_param_atom, only: zatom_md, nvmd, ntot, watom_md, omega_vmd
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
character(200) :: header
integer :: nheader, i, it, nfiles


!-----------------------------------------------------------------------
myidif: if( myid == 0 ) then

   !--- energies ---
   if( ifmd <= 1 ) then
       if( ioptmze_cell < 0 ) then
           header = '# step   P.E. [hartree]         '
           nheader = len_trim(header) + 3
         else
           header = '# step  Enthalpy [hr.]    P.E. [hr.]        PV [hr.]'
           nheader = len_trim(header) + 9
       end if
   else if( ifmd == 2 ) then
       header = (( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]                           ' )
       nheader = len_trim(header) + 4
   else if( ifmd == 3 ) then
       header = ((( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]    Total E.[hr.]   bathPE ' )//  &
&                   '[hr.]    bathKE [hr.]' )
       nheader = len_trim(header) + 4
   else if( ifmd == 4 ) then
       header =(((( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]    Total E.[hr.]   bathPE ' )//  &
&                   '[hr.]    bathKE [hr.]    cellPE')//  &
&                   ' [hr.]    cellKE [hr.]' )
       nheader = len_trim(header) + 4
   else if( ifmd == 5 ) then
       header =(((( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]    Total E.[hr.]   bathPE ' )//  &
&                   '[hr.]    bathKE [hr.]    bathPE')//  &
&                   '_atm [hr] bathKE_atm [hr]' )
       nheader = len_trim(header) + 4
   else if( ifmd == 6 ) then
       header =(((( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]    Total E.[hr.]   bathPE ' )//  &
&                   '[hr.]    bathKE [hr.]    ekindif')//  &
&                   ' [hr]    ' )
       nheader = len_trim(header) + 4
   else if( ifmd == 10 ) then
       header =((( '# step    H [hartree]     P.E. [' //  &
&                   'hartree]  K.E. [hartree]     T' )//  &
&                   '[K]    Total E.[hr.]   cellPE')//  &
&                   ' [hr.]    cellKE [hr.]' )
       nheader = len_trim(header) + 4
   end if
   if( lQMMDrun ) then
       header = ( header(1:nheader) //  &
&             'P.E.(MD) [hr.]  P.E.(QM) [hr.]  P.E.(MD cluster)' )
   end if
   write(nfile(3),'(200a1)')  &
&                 ( header(i:i), i = 1, len_trim(header) )


    !--- MD cell vectors ---
    write(nfile(4),2100) 'MD cell vectors               '
    write(nfile(4),2110) '        L_1         ',  &
&      '  L_2           L_3 ','         angle(2-3) ',  &
&      'angle(3-1) angle(1-2',')                   '

    write(nfile(9),'(a)') '# MD cell vectors'
    write(nfile(9),'(a)') '#        L_(1:3,1:3)'

    if( lstress ) then
        !--- stress ---
        write(nfile(16),2100) 'Stress in [GPa]               '
        write(nfile(16),2110) '        Pxx         ',  &
&      '    Pyy             ','Pzz             Pyz ',  &
&      '            Pzx     ','        Pxy         '

        write(nfile(22),2100) 'Diagonalized stress in [GPa]  '
        write(nfile(22),2110) '        Pxx         ',  &
&      '    Pyy             ','Pzz             Eige',  &
&      'nvectors            '
    end if

    if( ifmd >= 1 .and. ncbonds > 0 ) then
        write(nfile(18),2100) 'Lagrange multiplier for bond-l', &
&                             'ength constraint              '
!             write(nfile(18),'(i5)') ncbonds
    end if

    if( ifmd == 10 ) then
        write(nfile(23),'(a)') '# Hugoniot relations'
        write(nfile(23),'(a)')  &
& '# step Particle velocity   Pressure           Energy '
        write(nfile(23),'(a)')  &
& '#       Vs(1-rho0/rho)  rho0 Vs^2(1-rho0/rho) (P0+P)(1-rho0/rho)/rho0/2'
        write(nfile(23),'(a)')  &
& '#        [m/sec]             [GPa]            [a.u.]'
    end if

    if( lvmd ) then
!        call out_vmd( nfile(24) )
        write(nfile(24),'(a)') '# Potential energy of the reference ideal system'
        write(nfile(24),'(a,3es16.8)') '# (lambda) :', dlambda_vmd
        write(nfile(24),'(a,3es16.8)') '#  (k_B*T) :', treq
        write(nfile(24),'(a,3es16.8)') '#    (L_1) :', h(1:3,1,0)
        write(nfile(24),'(a,3es16.8)') '#    (L_2) :', h(1:3,2,0)
        write(nfile(24),'(a,3es16.8)') '#    (L_3) :', h(1:3,3,0)
        write(nfile(24),'(a,i3)') '# (type #) :', ntmd
        write(nfile(24),'(a)') '#          : it nvmd  ntot    watom_md        omega_vmd'
        do it = 1, ntmd
           write(nfile(24),'(a,i3,i3,i8,3es16.8)') '#          :', it,  &
& nvmd(it), ntot(it), watom_md(it), omega_vmd(it)
        end do
        write(nfile(24),'(a)') '# step  P.E. [hartree]  P.E.(ref) [hartree]'
    end if

!    if( lheat ) then
!!        write(nfile(25),'(a)') '# Heat flux vectors'
!!        write(nfile(25),'(a,a,a,a)') '#        ',  &
!!&                              'Kinetic energy transport (molecular motion)     ',  &
!!&                              'Potential energy transport (molecular motion)   ',  &
!!&                              'Energy transfer dut to molecular interaction'
!        write(nfile(25)) dtmd, nskip_heat, volume
!    end if

    if( ltotmom ) then
        write(nfile(26),'(a)') '# Reduced total momentum'
    end if

end if myidif


if( lmdout ) then

  ioif: if( mod(myid,iogpsz) == 0 ) then
     nfiles=min(iogpsz,nodes-myid)

    !--- atomic species ---
    write(nfile(5),'(a,i7)') '# Atomic species              ', nfiles
    write(nfile(5),'(i7,2x,100i4)')  &
& ntmd, ( nint(zatom_md(it)), it = 1, ntmd )

    !--- atomic configuration ---
    if( locoor ) then
        write(nfile(6),'(a,i7)') '# Atomic scaled coordinates   ', nfiles
    end if

    !--- atomic velocities ---
    if( ifmd >= 2 ) then
        if( lovelo ) then
            write(nfile(7),'(a,i7)') '# Atomic scaled velocities    ', nfiles
        end if
    end if

    !--- atomic forces ---
    if( loforc ) then
        write(nfile(8),'(a,i7)') '# Atomic scaled forces        ', nfiles
    end if

  end if ioif

end if


if( latomic ) then

  ioif2: if( mod(myid,iogpsz) == 0 ) then
     nfiles=min(iogpsz,nodes-myid)

    write(nfile(19),'(a,i7)') 'Atomic energies in [hartree] ', nfiles

    if( lstress ) then
        write(nfile(20),'(a,i7)') 'Atomic stresses from potential energy in [a.u.]  ', nfiles
        write(nfile(21),'(a,i7)') 'Atomic stresses from kinetic energy in [a.u.]    ', nfiles
    end if

  end if ioif2

end if


myidif2: if( myid == 0 ) then

    if( lmdout_qm ) then

        !--- atomic configuration for QM atoms ---
        if( locoor_qm ) then
           write(nfile(10),2100) 'Atomic scaled coordinates of Q',  &
&                                'M atoms                       '
        end if

        !--- atomic velocities for QM atoms ---
        if( ifmd >= 2 ) then
            if( lovelo_qm ) then
           write(nfile(11),2100) 'Atomic scaled velocities of QM',  &
&                                ' atoms                        '
            end if
        end if

        !--- atomic forces for QM atoms ---
        if( loforc_qm ) then
           write(nfile(12),2100) 'Atomic scaled forces of QM ato',  &
&                                'ms                            '
        end if

    end if

    if( lmdout_cl ) then

        !--- atomic configuration for MD-cluster atoms ---
        if( locoor_cl ) then
           write(nfile(13),2100) 'Atomic scaled coordinates of M',  &
&                                'D-cluster atoms               '
        end if

        !--- atomic velocities for MD-cluster atoms ---
        if( ifmd >= 2 ) then
            if( lovelo_cl ) then
           write(nfile(14),2100) 'Atomic scaled velocities of MD',  &
&                                '-cluster atoms                '
            end if
        end if

        !--- atomic forces for MD-cluster atoms ---
        if( loforc_cl ) then
           write(nfile(15),2100) 'Atomic scaled forces of MD-clu',  &
&                                'ster atoms                    '
        end if

    end if

end if myidif2
!-----------------------------------------------------------------------

!--- global maximum
call gimax(ierror)

2100 format('#',1x,10a30)
2110 format('#',10a20)
2111 format('#',10x,10a20)
2120 format('#',a20,i2)


return
end subroutine




module openpp_ext
!-----------------------------------------------------------------------
! type declaration and initialization of variables for atoms
!-----------------------------------------------------------------------
implicit none


logical :: ldirname_ext = .false.
character(1) :: dirname_ext

save

end module




subroutine openpp( nfile, myid, nodes,  &
& ifpp, jpl, zatom, fext, jgga, aname, ierror )
!-----------------------------------------------------------------------
!    file open for pseudopotentials
!
! (input)
!   ifpp = 1  : norm conserving pseudopotential
!        = 2  : ultrasoft pseudopotential
!        = 3  : data sets for the peojector-augmented wave (PAW) method
!        = 4  : relativistic norm conserving pseudopotential
!        = 5  : relativistic ultrasoft pseudopotential
!        = 6  : relativistic data sets for the PAW method
!       else  : not supported yet
!
!-----------------------------------------------------------------------
use outfile
use openpp_ext
implicit none
integer :: nfile(*), myid, nodes
integer :: ifpp
integer :: jpl
real*8  :: zatom
character(*) :: fext
integer :: jgga
character(2), dimension(103) :: aname
integer :: ierror

!-----declare local variables
integer :: nzn, nacr
character(50) :: fname, fnamebase
character(14), dimension(6) :: vandir =  &
&   (/  'control/NCPP/ ', 'control/USPP/ ', 'control/PAW/  ',  &
&       'control/RNCPP/', 'control/RUSPP/', 'control/RPAW/ ' /)
save vandir


!-----error trap
if( ifpp < 1 .or. ifpp > 6 ) then
    ierror = 100
    if(loutfile(1)) write(nfile(1),*) 'error in openpp: ifpp =', ifpp
    if(loutfile(2)) write(nfile(2),*) 'error in openpp: ifpp =', ifpp
    return
end if


nzn = nint(zatom)
nacr = 2
if( aname(nzn)(2:2) .eq. ' ' ) nacr = 1

fnamebase = trim(vandir(ifpp)) // aname(nzn)(1:nacr)
if( ldirname_ext ) fnamebase = trim(fnamebase) // dirname_ext
if( jgga == 1 ) then
    fname = trim(fnamebase) // '/'  &
&             // aname(nzn)(1:nacr) // '_' // trim(fext)
  else if( jgga == 2 ) then
    fname = trim(fnamebase) // '.PBE/'  &
&             // aname(nzn)(1:nacr) // '_' // trim(fext)
  else
    fname = trim(fnamebase) // '.RPBE/'  &
&             // aname(nzn)(1:nacr) // '_' // trim(fext)
endif
if(loutfile(1)) write(nfile(1),*) 'open file: ', trim(fname)
if(loutfile(2)) write(nfile(2),*) 'open file: ', trim(fname)

open( jpl, file=fname, status='old', action='read', iostat=ierror )

!--- global maximum
call gimax(ierror)

if( ierror /= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'cannot open file:',trim(fname)
    if(loutfile(2)) write(nfile(2),*) 'cannot open file:',trim(fname)
end if


return
end




subroutine existpp( nfile, myid, nodes,  &
& ifpp, zatom, fext, jgga, aname, lexist, ierror )
!-----------------------------------------------------------------------
!    inquire the existence of pseudopotential files
!
! (input)
!   ifpp = 1  : norm conserving pseudopotential
!        = 2  : ultrasoft pseudopotential
!        = 3  : data sets for the peojector-augmented wave (PAW) method
!       else  : not supported yet
!
! (output)
!   lexist = .true.  : exist
!          = .false. : not exist
!
!-----------------------------------------------------------------------
use outfile
use openpp_ext
implicit none
integer :: nfile(*), myid, nodes
integer :: ifpp
real*8  :: zatom
character(*) :: fext
integer :: jgga
character(2), dimension(103) :: aname
logical :: lexist
integer :: ierror

!-----declare local variables
integer :: nzn, nacr, iexist
character(50) :: fname, fnamebase
character(20), dimension(3) :: vandir =  &
&   (/ 'control/NCPP/', 'control/USPP/', 'control/PAW/ ' /)
save vandir


!-----error trap
if( ifpp < 1 .or. ifpp > 3 ) then
    ierror = 100
    if(loutfile(1)) write(nfile(1),*) 'error in existpp: ifpp =', ifpp
    if(loutfile(2)) write(nfile(2),*) 'error in existpp: ifpp =', ifpp
    return
end if


nzn = nint(zatom)
nacr = 2
if( aname(nzn)(2:2) .eq. ' ' ) nacr = 1

fnamebase = vandir(ifpp)(1:len_trim(vandir(ifpp))) // aname(nzn)(1:nacr)
if( jgga == 1 ) then
    fname = fnamebase(1:len_trim(fnamebase)) // '/'  &
&             // aname(nzn)(1:nacr) // '_' // fext(1:len_trim(fext))
  else if( jgga == 2 ) then
    fname = fnamebase(1:len_trim(fnamebase)) // '.PBE/'  &
&             // aname(nzn)(1:nacr) // '_' // fext(1:len_trim(fext))
  else
    fname = fnamebase(1:len_trim(fnamebase)) // '.RPBE/'  &
&             // aname(nzn)(1:nacr) // '_' // fext(1:len_trim(fext))
endif
if(loutfile(1)) write(nfile(1),*) 'inquire the existence of file: ',  &
&                     fname(1:len_trim(fname))
if(loutfile(2)) write(nfile(2),*) 'inquire the existence of file: ',  &
&                     fname(1:len_trim(fname))

inquire( file=fname, exist=lexist )

if( lexist ) then
    iexist = 0
  else
    iexist = 1
end if

!--- global maximum
call gimax(iexist)

lexist = iexist == 0

if( .not.lexist ) then
    if(loutfile(1)) write(nfile(1),*) 'not exist file:',fname(1:len_trim(fname))
    if(loutfile(2)) write(nfile(2),*) 'not exist file:',fname(1:len_trim(fname))
end if


return
end




module unit_numbers
!-----------------------------------------------------------------------
! type declaration of variables for unit numbers
!-----------------------------------------------------------------------
implicit none

integer, parameter :: unitmin = 10, unitmax = 9999     ! range of unit numbers

logical, dimension(unitmin:unitmax) :: unit_allocated = .false.

save

end module




subroutine allocate_unit_number( iunit )
!-----------------------------------------------------------------------
!  allocate unit number
!
!  return 0, if no unit number is available.
!-----------------------------------------------------------------------
use unit_numbers
implicit none
integer :: iunit

!-----declare local variables
integer :: i


iunit = 0
do i = unitmin, unitmax
   if( unit_allocated(i) ) cycle
   iunit = i
   unit_allocated(iunit) = .true.
   exit
end do


return
end subroutine




subroutine deallocate_unit_number( iunit )
!-----------------------------------------------------------------------
!  deallocate unit number
!-----------------------------------------------------------------------
use unit_numbers
implicit none
integer :: iunit


if( iunit >= unitmin .and. iunit <= unitmax ) then
    unit_allocated(iunit) = .false.
end if


return
end subroutine
