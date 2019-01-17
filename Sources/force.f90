



subroutine force( nfile, myid, nodes,  &
& t_comm, lcstress, ltimecnt, enclpv, elcl )
!-----------------------------------------------------------------------
!     force and internal stress tensor
!
!   floc : force by local pseudopotential
!   fnlc : force by nonlocal pseudopotential
!   fclm : force by direct Coulomb interaction
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: t_comm
logical :: lcstress, ltimecnt
real*8  :: enclpv, elcl

!-----declare local variables
!logical :: lgamma
!logical :: ltddft, ltddft_fssh
!logical :: lfssh_gsscf, ltddft_nscforce
real*8  :: elclv
real*8  :: bufst(6), bufst_kbppr(6), bufstr(6)
integer :: nspin, kvec, i
logical :: leig_end
logical :: lmkpdo_exit
real*8  :: ct0, ct, timecnt, tmrecip, tmreal
real*8  :: tmusrecip, tmusreal, tmcalsld, tmenlpstr, csum1, csum2


!call get_lgamma( lgamma )
!call get_ltddft( ltddft, ltddft_fssh )
!call get_lfssh_gsscf( lfssh_gsscf, ltddft_nscforce )

floc(:,:) = 0.d0
fnlc(:,:) = 0.d0
fclm(:,:) = 0.d0
fpcc(:,:) = 0.d0

strloc(:,:) = 0.d0
strnlc(:,:) = 0.d0
strclm(:,:) = 0.d0
strpcc(:,:) = 0.d0


if( ltimecnt ) then
    ct0 = timecnt()
end if

call fcolmb( nfile, myid, nodes,                                   &
& ecolmb, atom_ecolmb, fclm, strclm,  &
& lclust, ntype, nhk1, nhk2, zv, nel, ratm,                         &
& iatoit, iatmpt_nod, nion_nod, natom, mx1loc, hcell, gamma, rccc2, &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& bufatm, bufatmr, .false., lcstress,  &
& ldouble_grid_recip, idouble_grid_method, lvacuum )

if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '                  direct Coulomb ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '                  direct Coulomb ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if


!-----force by local pseudopotential in real space
if( llclpp_r ) then
    call flocal( nfile, myid, nodes,                               &
& floc, strloc, lclust, nd1v, ntype, nhk1, nhk2,  &
& nhk1_nod, nhk2_nod, zv, nel, ratm, iatoit, iatmpt_nod, nion_nod,  &
& natom, llking, mx1loc, tbflc, tbflca, dltflc, rmxflc,  &
& hcell, lorthrhmbc, lhfull, rdel, rdelv, rho, noddatx, .false.,  &
& lcstress,  &
& mshxyz(mulpit(1)),  &
& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),  &
& tmpk, tmpl, tmpm, tmpn, bufatm, bufatmr )
end if

!-----force by local pseudopotential in reciprocal space
if( llclpp_g ) then
    if( lcstress )  &
&       call calc_vext( nfile, myid, nodes,  &
& elclv, elcl, rho, mshnod(1), rdelv )
    call flocal_g( nfile, myid, nodes,  &
& floc, strloc, lcstress, ntype, nhk1_nod, nhk2_nod, iatoit,  &
& iatmpt_nod, natom, nion_nod, llking, nplw5, nplw5ex, nplw,  &
& nga, ngb, ngc, gx, gy, gz, rhgr, elclv, & !ycos, ysin, nplwcs,  &
& kmax1cs, kmax2cs, kmax3cs, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& .false. )
end if

!-----for the double-grid method
!if( ldouble_grid_recip ) then
!    call flocal_dg( nfile, myid, nodes,  &
!& floc, ntype, nhk1_nod, nhk2_nod, iatoit,  &
!& iatmpt_nod, natom, nion_nod, zv, .false. )
!end if

!-----forces from local pseudopotential
if( llclpp_r .or. llclpp_g )  &
& call gdsum(floc,3*natom,bufatmr)

if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '           local pseudopotential ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '           local pseudopotential ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if


!tddftif: if( .not.ltddft ) then
!if( lvand_r .or. lvand_g ) then
!if( lfssh_gsscf .and. ltddft_nscforce ) then
!    hdiag(1:mshnod(1)) = vhar_out(1:mshnod(1)) - vhar(1:mshnod(1)) + xexc(1:mshnod(1))
!    call unifylc( nfile, myid, nodes,  &
!& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!
!    call calhartdij( nfile, myid, nodes,  &
!& nplw3, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& glocal, eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm )
!
!    if( ltimecnt ) then
!        ct = timecnt()
!        if(loutfile(1)) write(nfile(1),*) '                      calhartdij ',  &
!&                             ':          :', ct - ct0
!        if(loutfile(2)) write(nfile(2),*) '                      calhartdij ',  &
!&                             ':          :', ct - ct0
!        ct0 = ct
!    end if
!end if
!end if


      bufst(1:6) = 0.d0
bufst_kbppr(1:6) = 0.d0
lmkpdo_exit = .true.

!if( lvand_r .or. lvand_g ) then
!    call clear_fbbx( nfile, myid, nodes, natom )
!end if

!if( lkbpp_g .or. lvand_r .or. lvand_g ) then
spindo: do nspin = 1, nspnmx
   if( lspin ) then
       call ldocck( nspin, occ, wegud, nband, nkpnt )
       if( .not.ltddft .and. lgamma ) then
           call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand_r .or. lvand_g ) then
!               call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!           end if
       end if
!       if( lvand_r .or. lvand_g ) then
!           call ldvloc( nspin, vlocud, vexc, mshnod(1) )
!           call ldeig_k( nspin, eig, egvlud, nband, nkpnt )
!           call lddij( nfile, myid, nodes,  &
!& nspin, ntype, nhk1, nhk2, natom, lvandi, rvol )
!       end if
   end if

    !----to keep the original force in scissors corrections
!    call restore_original_eig( nfile, myid, nodes,  &
!& eig, nspin, nband )

    !-----norm-conserving pseudopotentials in KB form
    if( lkbpp_r .or. lkbpp_g ) then

        tmrecip = 0.d0
        tmreal  = 0.d0
        kvecdo: do kvec = 1, nknod
!           if( .not.lgamma ) then
!               call set_kvec( kvec )
!               call get_wfk( rhcr, nspin )
!               call get_plwk( nplw, nplwex, nspnod )
!           end if
!           call loadsumr_g_org( nfile, myid, nodes,  &
!& nbnod, ntype, nhk1, nhk2, natom, lkbppi, lking, lspin, nspin )
           call loadsumr_g( nfile, myid, nodes,  &
& ntype, nhk1, nhk2, lkbppi, lking, lspin, nspin,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& iods, iodsg, iod, iodg, ioag, idstnd )

           leig_end   = kvec == nknod .and. nspin == nspnmx
           call set_all2all_out( leig_end )

           !-----norm-conserving pp. in real space
           if( lkbpp_r ) then

!               if( lgamma ) then
                   call fnnlcl( nfile, myid, nodes,  &
& fnlc, bufst_kbppr, .false., lcstress, t_comm,  &
& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
& lkbppi, lking, rdel_c, rdelg_c, rdelv_c, hcell, hci, lorthrhmbc,  &
& mshglb_c, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
& iods, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkg, nplw2,  &
& nspin, bufatm, bufatmr, natom, mxl, eigr, eigi )
!               else
!                   call set_eikr( nfile, myid, nodes,  &
!& mshglb_c, kfft1, kfft2, kfft3 )
!                   call set_eikl( nfile, myid, nodes,  &
!& ratm, natom, lkbpp_r, .false. )
!                   call fnnlcl_k( nfile, myid, nodes,  &
!& fnlc, bufst_kbppr, .false., lcstress, t_comm,  &
!& ntype, nhk1, nhk2, lmax, lclno, lchk, iatoit,  &
!& lkbppi, lking, rdel_c, rdelg_c, rdelv_c, hcell, hci, lorthrhmbc,  &
!& mshglb_c, kfft1, kfft2, kfft3, mshgnx, mshgny, mshgnz,  &
!& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ,  &
!& iods, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2,  &
!& nspin, bufatm, bufatmr, natom, mxl, nylmmx, mx1, eigr, eigi )
!               end if

               if( ltimecnt ) then
                   ct = timecnt()
                   tmreal = tmreal + ct - ct0
!                   if( myid.eq.0 ) then
!                       write(nfile(1),*) '              NCPP in real space ',  &
!&                                        ':          :', ct - ct0
!                       write(nfile(2),*) '              NCPP in real space ',  &
!&                                        ':          :', ct - ct0
!                   end if
                   ct0 = ct
               end if

           end if

           !-----norm-conserving pp. in reciprocal space
           if( lkbpp_g ) then

!               if( lgamma ) then
                   call kbnlpp_frc( nfile, myid, nodes,  &
& fnlc, bufst, lcstress, rhcr, nplwex, nplw, gx, gy, gz,  &
& nband, nbnod1, nbnod2, nbnod, occ, iods,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
!               else
!                   call set_ycosysin( nfile, myid, nodes,  &
!& natom, ycos, ysin, nplwcs )
!                   call kbnlpp_frc_k( nfile, myid, nodes,  &
!& fnlc, bufst, lcstress, rhcr, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods,  &
!& ntype, nhk1, nhk2, natom, lkbppi, lking )
!               end if

               if( ltimecnt ) then
                   ct = timecnt()
                   tmrecip = tmrecip + ct - ct0
!                   if( myid.eq.0 ) then
!                       write(nfile(1),*) '        NCPP in reciprocal space ',  &
!&                                        ':          :', ct - ct0
!                       write(nfile(2),*) '        NCPP in reciprocal space ',  &
!&                                        ':          :', ct - ct0
!                   end if
                   ct0 = ct
               end if

           end if

        end do kvecdo

        if( ltimecnt ) then
            if( lkbpp_g ) then
                if(loutfile(1)) write(nfile(1),*) '        NCPP in reciprocal space ',  &
&                                 ':          :', tmrecip
                if(loutfile(2)) write(nfile(2),*) '        NCPP in reciprocal space ',  &
&                                 ':          :', tmrecip
            end if

            if( lkbpp_r ) then
                if(loutfile(1)) write(nfile(1),*) '              NCPP in real space ',  &
&                                 ':          :', tmreal
                if(loutfile(2)) write(nfile(2),*) '              NCPP in real space ',  &
&                                 ':          :', tmreal
            end if
        end if

    end if


    !-----ultrasoft pseudopotentials
!    if( lvand_r .or. lvand_g ) then

        !--- unify local potential
!        hdiag(1:mshnod(1)) = vexc(1:mshnod(1)) + vhar(1:mshnod(1)) + vhshdp(1:mshnod(1)) + vext(1:mshnod(1))
!        call unifylc( nfile, myid, nodes,  &
!& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!        call pwpunifylc( nfile, myid, nodes,  &
!& hdiag, vexc, vhar, vhshdp, vext, mshnod(1), glocal, vlocud,  &
!& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d )
!
!        if( ltimecnt ) then
!            ct = timecnt()
!            if(loutfile(1)) write(nfile(1),*) '                         unifylc ',  &
!&                                 ':          :', ct - ct0
!            if(loutfile(2)) write(nfile(2),*) '                         unifylc ',  &
!&                                 ':          :', ct - ct0
!            ct0 = ct
!        end if

!        lmkpdo_exit = .true.
!        if( .not.lgamma ) then
!            call symk_frc_clear( nfile, myid, nodes, lmkpdo_exit )
!        end if

        !-----contribution from Q functions in ultrasoft pp
!        if( lmkpdo_exit ) then
!            call caldij( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, nplw3,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .true., .true. )
!            call calfrcdij( nfile, myid, nodes,  &
!& lspin, nspin, nplw3, gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .true. )
!
!            if( lfssh_gsscf .and. ltddft_nscforce ) then
!                !--- unify local potential
!                hdiag(1:mshnod(1)) = vhar_out(1:mshnod(1)) - vhar(1:mshnod(1)) + xexc(1:mshnod(1))
!                call unifylc( nfile, myid, nodes,  &
!& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!
!                if( lspin ) then
!                    call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                    call ldocck( nspin, occ, wegud, nband, nkpnt )
!                else
!                    call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                end if
!
!                call caldij( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, nplw3,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .true., .true. )
!                call calfrcdij( nfile, myid, nodes,  &
!& lspin, nspin, nplw3, gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .true. )
!
!                !--- unify local potential
!                hdiag(1:mshnod(1)) = vexc(1:mshnod(1)) + vhar(1:mshnod(1))  &
!&                                + vhshdp(1:mshnod(1)) + vext(1:mshnod(1))
!                call unifylc( nfile, myid, nodes,  &
!& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!
!                if( lspin ) then
!                    call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                    call ldocck( nspin, occ, wegud, nband, nkpnt )
!                else
!                    call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                end if
!            end if
!
!            if( ltimecnt ) then
!                ct = timecnt()
!                if(loutfile(1)) write(nfile(1),*) '                       calscrdij ',  &
!&                                 ':          :', ct - ct0
!                if(loutfile(2)) write(nfile(2),*) '                       calscrdij ',  &
!&                                 ':          :', ct - ct0
!                ct0 = ct
!            end if
!        end if
!
!
!     tmusrecip = 0.d0
!     tmusreal  = 0.d0
!     tmcalsld  = 0.d0
!     tmenlpstr = 0.d0
!     csum1     = 0.d0
!     csum2     = 0.d0
!     kvecdo2: do kvec = 1, nknod
!        if( .not.lgamma ) then
!            call set_kvec( kvec )
!            call get_wfk( rhcr, nspin )
!            call get_plwk( nplw, nplwex, nspnod )
!            call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!        end if

!        lmkpdo_exit = .true.
!        if( .not.lgamma ) then
!        mkpif: if( .not.lmkpdo_exit ) then
!            call symk_frc( nfile, myid, nodes,  &
!& fnlc, bufst, lcstress, ltimecnt, ct0,  &
!& tmusrecip, tmusreal, tmcalsld, tmenlpstr, csum1, csum2,  &
!& ycos, ysin, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
!& rhcr, nspin, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods, eig, glocal )

!        else mkpif    !: if( lmkpdo_exit ) then

!        leig_end   = kvec == nknod .and. nspin == nspnmx
!        call set_all2all_out( leig_end )
!        !--- convert atom decomposition to band decomposition 
!        call slm_gdtobd( nfile, myid, nodes, ct0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iods, iodsg, ioag, idstnd, ltimecnt )
!
!        if( lvand_g ) then
!
!            if( .not.ltddft ) then
!
!                if( lgamma ) then
!                    call vdnlpp_frc( nfile, myid, nodes,  &
!& fnlc, nspin, rhcr, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs,  &
!& fft3x, fft3y )
!
!                    if( lfssh_gsscf .and. ltddft_nscforce ) then
!                        if( lspin ) then
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                            call ldocck( nspin, occ, wegud, nband, nkpnt )
!                        else
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                        end if
!
!                        call vdnlpp_nscfrc( nfile, myid, nodes,  &
!& fnlc, nspin, rhcr, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs,  &
!& fft3x )
!
!                        if( lspin ) then
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                        else
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                        end if
!                    end if
!
!                else
!                    call set_ycosysin( nfile, myid, nodes,  &
!& natom, ycos, ysin, nplwcs )
!                    call vdnlpp_frc_k( nfile, myid, nodes,  &
!& fnlc, nspin, rhcr, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& fft3x, fft3y )
!                end if
!
!              else
!
!                    call vdnlpp_frc_in_tddft( nfile, myid, nodes,  &
!& fnlc, nspin, nplwex, nplw, gx, gy, gz,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs, fft3x, fft3y, &
!& bijr, biji, bm1r, bm1i )
!
!            end if
!
!            if( ltimecnt ) then
!                ct = timecnt()
!                tmusrecip = tmusrecip + ct - ct0
!                      if( myid.eq.0 ) then
!                  write(nfile(1),*) '        USPP in reciprocal space ',
!     &                          ':          :', ct - ct0
!                  write(nfile(2),*) '        USPP in reciprocal space ',
!     &                          ':          :', ct - ct0
!                      end if
!                ct0 = ct
!            end if
!
!        end if
!
!
!        if( lvand_r ) then
!
!            if( .not.ltddft ) then
!
!                if( lgamma ) then
!                    call vdnlpp_frc_r( nfile, myid, nodes,  &
!& fnlc, lcstress, rhcr, nplw, nplwex,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, lmax, mxl,  &
!& rdel_c, rdelg_c, hcell, hci, mshglb_c, kfft1, kfft2, kfft3,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, mshgnx, mshgny, mshgnz,  &
!& .false. )
!
!                    if( lfssh_gsscf .and. ltddft_nscforce ) then
!                        if( lspin ) then
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                            call ldocck( nspin, occ, wegud, nband, nkpnt )
!                        else
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                        end if
!
!                        call vdnlpp_nscfrc_r( nfile, myid, nodes,  &
!& fnlc, rhcr, nplw, nplwex,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, lmax, mxl,  &
!& rdel_c, rdelg_c, hcell, hci, mshglb_c, kfft1, kfft2, kfft3,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr,  &
!& .false. )
!
!                        if( lspin ) then
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!                        else
!                            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!                        end if
!                    end if
!
!                  else
!                    call set_eikr( nfile, myid, nodes,  &
!& mshglb_c, kfft1, kfft2, kfft3 )
!                    call set_eikl( nfile, myid, nodes,  &
!& ratm, natom, .false., lvand_r )
!                    call vdnlpp_frc_r_k( nfile, myid, nodes,  &
!& fnlc, lcstress, rhcr, nplw, nplwex,  &
!& nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, lmax, mxl,  &
!& rdel_c, rdelg_c, hcell, hci, mshglb_c, kfft1, kfft2, kfft3,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2,  &
!& eigr, eigi, mshgnx, mshgny, mshgnz,  &
!& .false. )
!                end if
!
!              else
!
!                    call vdnlpp_frc_r_in_tddft( nfile, myid, nodes,  &
!& fnlc, nspin, lcstress, nplw, nplwex,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, lmax, mxl,  &
!& rdel_c, rdelg_c, hcell, hci, mshglb_c, kfft1, kfft2, kfft3,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2,  &
!& eigr, eigi, mshgnx, mshgny, mshgnz,  &
!& bijr, biji, bm1r, bm1i, &
!& .false. )
!
!            end if
!
!            if( ltimecnt ) then
!                ct = timecnt()
!                tmusreal  = tmusreal + ct - ct0
!!                      if( myid.eq.0 ) then
!!                  write(nfile(1),*) '              USPP in real space ',
!!     &                          ':          :', ct - ct0
!!                  write(nfile(2),*) '              USPP in real space ',
!!     &                          ':          :', ct - ct0
!!                      end if
!                ct0 = ct
!            end if
!
!        end if
!
!
!        !--- stress calculation
!        if( lcstress ) then
!
!            if( lvand_g ) then
!
!                if( lgamma ) then
!                    call calsld( nfile, myid, nodes,  &
!& rhcr, nplw, nplwex, gx, gy, gz, nband, nbnod1, nbnod2, nbnod,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )
!                  else
!                    call calsld_k( nfile, myid, nodes,  &
!& rhcr, nplw, nplwex, gx, gy, gz, nband, nbnod1, nbnod2, nbnod,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking )
!                end if
!
!                if( ltimecnt ) then
!                    ct = timecnt()
!                    tmcalsld  = tmcalsld  + ct - ct0
!!                          if( myid.eq.0 ) then
!!                  write(nfile(1),*) '                          calsld ',
!!     &                          ':          :', ct - ct0
!!                  write(nfile(2),*) '                          calsld ',
!!     &                          ':          :', ct - ct0
!!                          end if
!                    ct0 = ct
!                end if
!
!            end if
!
!            if( lgamma ) then
!                call enlpstr( nfile, myid, nodes,  &
!& nspin, lpaw, bufst, nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi )
!              else
!                call enlpstr_k( nfile, myid, nodes,  &
!& nspin, lpaw, bufst, nband, nbnod1, nbnod2, nbnod, occ, iods, eig,  &
!& ntype, nhk1, nhk2, natom, lvandi )
!            end if
!
!            if( ltimecnt ) then
!                ct = timecnt()
!                tmenlpstr = tmenlpstr + ct - ct0
!!                      if( myid.eq.0 ) then
!!                  write(nfile(1),*) '                         enlpstr ',
!!     &                          ':          :', ct - ct0
!!                  write(nfile(2),*) '                         enlpstr ',
!!     &                          ':          :', ct - ct0
!!                      end if
!                ct0 = ct
!            end if
!
!            if( lgamma ) then
!                call chgfrdq( nfile, myid, nodes,  &
!& bufst, nplw3, gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, idstnd,  &
!& occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, iodsg, ioa, ioag, lvandi, ltimecnt, ct0, csum1, csum2 )
!            else
!                call chgfrdq_k( nfile, myid, nodes,  &
!& bufst, nplw3, gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, idstnd,  &
!& occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, iodsg, ioa, ioag, lvandi, ltimecnt, ct0, csum1, csum2 )
!            end if
!
!        end if
!
!
!        end if mkpif
!
!        if( .not.lgamma ) then
!            call svslmi_k( nfile, myid, nodes, lvandi, nspin, .true. )
!        end if
!
!     end do kvecdo2
!
!        !-----contribution from Q functions in ultrasoft pp
!        if( .not.lmkpdo_exit ) then
!
!            call symk_frc_gsum_fbbfr( nfile, myid, nodes,  &
!& ntype, nhk1, nhk2, natom, iatoit, ioa, lvandi )

!            call caldij( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, nplw3,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .true., .false. )
!            call calfrcdij( nfile, myid, nodes,  &
!& lspin, nspin, nplw3, gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .false. )
!
!            if( ltimecnt ) then
!                ct = timecnt()
!                if(loutfile(1)) write(nfile(1),*) '                       calscrdij ',  &
!&                                 ':          :', ct - ct0
!                if(loutfile(2)) write(nfile(2),*) '                       calscrdij ',  &
!&                                 ':          :', ct - ct0
!                ct0 = ct
!            end if
!        end if

!     if( lspin ) then
!         if( lgamma ) then
!             call svslmi( nfile, myid, nodes, nspin, lvandi, .true. )
!         end if
!     end if

!     if( ltimecnt ) then
!     do i = 1, 2
!     if( loutfile(i) ) then
!        if( lvand_g ) then
!            write(nfile(i),*) '        USPP in reciprocal space ',  &
!&                         ':          :', tmusrecip
!        end if
!
!        if( lvand_r ) then
!            write(nfile(i),*) '              USPP in real space ',  &
!&                         ':          :', tmusreal
!        end if
!
!        if( lcstress ) then
!            if( lvand_g ) then
!            write(nfile(i),*) '                          calsld ',  &
!&                         ':          :', tmcalsld
!            end if
!
!            write(nfile(i),*) '                         enlpstr ',  &
!&                         ':          :', tmenlpstr
!            write(nfile(i),*) '         slmidr_bd2ad in chgfrdq ',  &
!&                         ':          :', csum1
!            write(nfile(i),*) '                         chgfrdq ',  &
!&                         ':          :', csum2
!        end if
!     end if
!     end do
!     end if
!
!    end if

    !---scissors corrections to eigenvalues at donor/acceptor interface
!    call reset_scissors_eig( nfile, myid, nodes,  &
!& eig, nspin, nband )

end do spindo
!end if

!if( lvand_r .or. lvand_g ) then
!    !-----add contribution from Q functions to fnlc
!    call addfbbx( nfile, myid, nodes, natom, fnlc )
!end if


!-----forces from nonlocal pseudopotential
call gdsum(fnlc,3*natom,bufatmr)
call unify_sumn(fnlc,3*natom,bufatmr)


!---symmetry operations for nonlocal forces and nonlocal stress
!if( lmkpdo_exit ) then
!    call symk_symfrc( nfile, myid, nodes,  &
!& natom, ratm, fnlc, hcell )
!
!    call strkine_k2_sym( nfile, myid, nodes, bufst_kbppr )
!    call strkine_k2_sym( nfile, myid, nodes, bufst )
!end if

!-----stress from nonlocal pseudopotential
if( lcstress ) then
    call gdsum(bufst_kbppr,6,bufstr)
    call unify_sumn(bufst_kbppr,6,bufstr)
    call gdsum(bufst,6,bufstr)
    call unify_sumn(bufst,6,bufstr)
    strnlc(1,1) = ( bufst_kbppr(1) - bufst(1) )/rvol - 2.d0*enclpv
    strnlc(2,2) = ( bufst_kbppr(2) - bufst(2) )/rvol - 2.d0*enclpv
    strnlc(3,3) = ( bufst_kbppr(3) - bufst(3) )/rvol - 2.d0*enclpv
    strnlc(2,3) = ( bufst_kbppr(4) - bufst(4) )/rvol
    strnlc(3,1) = ( bufst_kbppr(5) - bufst(5) )/rvol
    strnlc(1,2) = ( bufst_kbppr(6) - bufst(6) )/rvol
    strnlc(3,2) = strnlc(2,3)
    strnlc(1,3) = strnlc(3,1)
    strnlc(2,1) = strnlc(1,2)
end if


if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '     gdsum for nonlocal pp. etc. ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '     gdsum for nonlocal pp. etc. ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if

!end if tddftif


!-----forces from the partial core correction
if( lpcc ) then
    if( lspin .or. lnoncollinear ) then
        do i = 1, mshnod(1)
           vexc(i) = ( vlocud(i) + vlocud(i+mshnod(1)) )*0.5d0
        end do
    end if

    hdiag(1:mshnod(1)) = vexc(1:mshnod(1))
    if( lfssh_gsscf .and. ltddft_nscforce ) then
        hdiag(1:mshnod(1)) = hdiag(1:mshnod(1)) + xexc(1:mshnod(1))
    end if

    !-----in real space
    if( lpcc_r ) then
        call pccfrc( nfile, myid, nodes,  &
& hdiag, fpcc, lclust, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& lvacuum,  &
& ntype, natom, nhk1, nhk2, ratm, ltimecnt, bufatm, bufatmr,  &
& lpcci, rpcc, lpking,  &
& mshxyz(mulpit(1)),  &
& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
    end if

    !-----in reciprocal space
    if( lpcc_g ) then
        call unifylc( nfile, myid, nodes,  &
& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
        call pccfrc_g( nfile, myid, nodes,  &
& glocal, fpcc, strpcc, lcstress, ntype, nhk1_nod, nhk2_nod, iatoit,  &
& iatmpt_nod, natom, nion_nod, lpcci, lpking, nplw5, nplw5ex, nplw,  &
& nga, ngb, ngc, gx, gy, gz, thrhgr,  & !ycos, ysin, nplwcs,  &
& kmax1cs, kmax2cs, kmax3cs, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& nd1vks(1), nd1vks(2), nd1vks(3), kfft0d, mfd2ft, ntotfd, ijkgd, nplw5,  &
& fft3x, fft3y, fftwork, .false. )
    end if

    call gdsum(fpcc,3*natom,bufatmr)

    !---symmetry operations on PCC force
!    call symk_symfrc( nfile, myid, nodes,  &
!& natom, ratm, fpcc, hcell )
end if


!-----forces from an empirical correction for vdW : DFT-D
if( ldftd ) call force_dftd( nfile, myid, nodes, &
& lclust, ntype, nhk1, nhk2, ratm, &
& iatoit, iatmpt_nod, nion_nod, natom, hcell, &
& bufatm, bufatmr, ltimecnt, lcstress,  &
& ldouble_grid_recip, idouble_grid_method, lvacuum )


do i = 1, nhk2(ntype)
   frc(1:3,i) = floc(1:3,i) + fnlc(1:3,i) + fclm(1:3,i) + fpcc(1:3,i)
end do

!-----forces from an empirical correction for vdW : DFT-D
if( ldftd ) call add_force_dftd( nfile, myid, nodes, frc, natom )

!-----force from uniform electric field
!if( lefield ) then
!    call force_efield( nfile, myid, nodes, &
!& efield, lclust, ntype, nhk1, nhk2, ratm, iatoit, natom, zv, hcell,  &
!& ltimecnt, lcstress, lefield_islts, lsawtooth_shape, lsawtooth_xyz )
!
!    if( .not.lsawtooth_shape ) then
!        call force_efield_enstr( nfile, myid, nodes, &
!& efield, hcell, ltimecnt, lcstress, lvshape, lefield_islts, lconstraintD )
!    end if
!end if


return
end




subroutine corfrc( ntype, nhk1, nhk2, frc, nt )
!-----------------------------------------------------------------------
!     correction
!-----------------------------------------------------------------------
implicit none

integer :: ntype
integer :: nt
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,nt)   :: frc
real*8,  dimension(3)      :: sumf

integer :: ix, i

do ix = 1, 3
   sumf(ix) = 0.d0
   do i = 1, nhk2(ntype)
      sumf(ix) = sumf(ix) + frc(ix,i)
   end do
   sumf(ix) = sumf(ix)/dble(nhk2(ntype))
   do i = 1, nhk2(ntype)
      frc(ix,i) = frc(ix,i) - sumf(ix)
   end do
end do

return
end




subroutine fcolmb( nfile, myid, nodes,                             &
& ecolmb, atom_ecolmb, fclm, strclm,  &
& lclust, ntype, nhk1, nhk2, zv, nel, ratm,  &
& iatoit, iatmpt, nion, natom, mx1, hcell, gamma, rccc2,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& dbuf, dbufr, ltimecnt, lcstress,  &
& ldouble_grid_recip, idouble_grid_method, lvacuum )
!-----------------------------------------------------------------------
!   fclm : energy and force by direct Coulomb interaction
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: ecolmb
real*8,  dimension(natom)   :: atom_ecolmb
real*8,  dimension(3,natom) :: fclm
real*8,  dimension(3,3)  :: strclm
logical :: lclust
integer :: ntype
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom)   :: ratm
integer, dimension(natom)     :: iatoit
integer, dimension(nion)  :: iatmpt
integer :: nion
real*8,  dimension(3,3) :: hcell
real*8,  dimension(0:mx1) :: tberf, tberfa
real*8 :: drcut
real*8,  dimension(0:mx1) :: tbeff, tbeffa
real*8 :: drcutf

real*8,  dimension(6) :: bufst, bufstr
real*8,  dimension(3*natom) :: dbuf, dbufr
logical   ltimecnt, lcstress, ldouble_grid_recip
integer :: idouble_grid_method
logical, dimension(3) :: lvacuum


if( ltimecnt ) then
    ct0 = timecnt()
end if
ecolmb = 0.d0
do i = 1, natom
   atom_ecolmb(i) = 0.d0
end do
if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
do i_node  = 1, nion
   i   = iatmpt(i_node)
   it1 = iatoit(i)
   do j  = 1, i - 1
      it2 = iatoit(j)
      zv2 = zv(it1)*zv(it2)*2.d0
      xx = ratm(1,j) - ratm(1,i)
      yy = ratm(2,j) - ratm(2,i)
      zz = ratm(3,j) - ratm(3,i)
      r2 = xx*xx + yy*yy + zz*zz
      rr = sqrt( r2 )
      r3 = r2*rr

      tecolmb = zv2/rr
      ecolmb = ecolmb + tecolmb
      atom_ecolmb(i) = atom_ecolmb(i) + tecolmb
      atom_ecolmb(j) = atom_ecolmb(j) + tecolmb

      r3 = zv2/r3
      ffx =  - xx*r3
      ffy =  - yy*r3
      ffz =  - zz*r3
      fclm(1,i) = fclm(1,i) + ffx
      fclm(2,i) = fclm(2,i) + ffy
      fclm(3,i) = fclm(3,i) + ffz
      fclm(1,j) = fclm(1,j) - ffx
      fclm(2,j) = fclm(2,j) - ffy
      fclm(3,j) = fclm(3,j) - ffz
   end do
end do
else
!--- for bulk calculation ----------------------------------------------
do i = 1, 6
   bufst(i) = 0.d0
end do
gamma2 = gamma*gamma
gamma3 = gamma2*gamma
if( lcstress ) then
!--- when stress is needed ------------------
do i_node  = 1, nion
   i   = iatmpt(i_node)
   it1 = iatoit(i)
   do j  = 1, i - 1
      it2 = iatoit(j)
      zv2 = zv(it1)*zv(it2)*2.d0
      q1 = ratm(1,j) - ratm(1,i)
      q2 = ratm(2,j) - ratm(2,i)
      q3 = ratm(3,j) - ratm(3,i)
      if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.0d0,q1)
      if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.0d0,q2)
      if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.0d0,q3)
!            xx = hcell(1,1)*q1
!            yy = hcell(2,2)*q2
!            zz = hcell(3,3)*q3
      xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
      yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
      zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
      r2 = xx*xx + yy*yy + zz*zz
      if( r2.le.rccc2 ) then
          r2g = r2*gamma2
          m = r2g/(2.d0*drcut)
          m = 2*m
          d = 0.5d0*( r2g/drcut - dble(m) )
          cerf = d*( (d-1.d0)*tberfa(m) + tberf(m+2)  &
&                        - tberf(m)   ) + tberf(m)
          cerff= d*( (d-1.d0)*tbeffa(m) + tbeff(m+2)  &
&                        - tbeff(m)   ) + tbeff(m)

          r2 = 1.d0/r2
          rr = sqrt(r2)

          tecolmb = zv2*( rr - gamma*cerf )
          ecolmb = ecolmb + tecolmb
          atom_ecolmb(i) = atom_ecolmb(i) + tecolmb
          atom_ecolmb(j) = atom_ecolmb(j) + tecolmb

          r3 = zv2*( r2*rr - gamma3*cerff )
          ffx =  - xx*r3
          ffy =  - yy*r3
          ffz =  - zz*r3
          fclm(1,i) = fclm(1,i) + ffx
          fclm(2,i) = fclm(2,i) + ffy
          fclm(3,i) = fclm(3,i) + ffz
          fclm(1,j) = fclm(1,j) - ffx
          fclm(2,j) = fclm(2,j) - ffy
          fclm(3,j) = fclm(3,j) - ffz

          bufst(1) = bufst(1) + ffx*xx
          bufst(2) = bufst(2) + ffy*yy
          bufst(3) = bufst(3) + ffz*zz
          bufst(4) = bufst(4) + ffy*zz
          bufst(5) = bufst(5) + ffz*xx
          bufst(6) = bufst(6) + ffx*yy
      end if
   end do
end do
call gdsum(bufst,6,bufstr)
strclm(1,1) = bufst(1)
strclm(2,2) = bufst(2)
strclm(3,3) = bufst(3)
strclm(2,3) = bufst(4)
strclm(3,1) = bufst(5)
strclm(1,2) = bufst(6)
else if( idouble_grid_method /= 3 .or. .not.ldouble_grid_recip ) then
!--- when stress is not needed ------------------
do i_node  = 1, nion
   i   = iatmpt(i_node)
   it1 = iatoit(i)
   do j  = 1, i - 1
      it2 = iatoit(j)
      zv2 = zv(it1)*zv(it2)*2.d0
      q1 = ratm(1,j) - ratm(1,i)
      q2 = ratm(2,j) - ratm(2,i)
      q3 = ratm(3,j) - ratm(3,i)
      if( abs(q1).gt.0.5d0 ) q1 = q1 - sign(1.0d0,q1)
      if( abs(q2).gt.0.5d0 ) q2 = q2 - sign(1.0d0,q2)
      if( abs(q3).gt.0.5d0 ) q3 = q3 - sign(1.0d0,q3)
!            xx = hcell(1,1)*q1
!            yy = hcell(2,2)*q2
!            zz = hcell(3,3)*q3
      xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
      yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
      zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
      r2 = xx*xx + yy*yy + zz*zz
      if( r2.le.rccc2 ) then
          r2g = r2*gamma2
          m = r2g/(2.d0*drcut)
          m = 2*m
          d = 0.5d0*( r2g/drcut - dble(m) )
          cerf = d*( (d-1.d0)*tberfa(m) + tberf(m+2)  &
&                        - tberf(m)   ) + tberf(m)
          cerff= d*( (d-1.d0)*tbeffa(m) + tbeff(m+2)  &
&                        - tbeff(m)   ) + tbeff(m)

          r2 = 1.d0/r2
          rr = sqrt(r2)

          tecolmb = zv2*( rr - gamma*cerf )
          ecolmb = ecolmb + tecolmb
          atom_ecolmb(i) = atom_ecolmb(i) + tecolmb
          atom_ecolmb(j) = atom_ecolmb(j) + tecolmb

          r3 = zv2*( r2*rr - gamma3*cerff )
          ffx =  - xx*r3
          ffy =  - yy*r3
          ffz =  - zz*r3
          fclm(1,i) = fclm(1,i) + ffx
          fclm(2,i) = fclm(2,i) + ffy
          fclm(3,i) = fclm(3,i) + ffz
          fclm(1,j) = fclm(1,j) - ffx
          fclm(2,j) = fclm(2,j) - ffy
          fclm(3,j) = fclm(3,j) - ffz

!                bufst(1) = bufst(1) + ffx*xx
!                bufst(2) = bufst(2) + ffy*yy
!                bufst(3) = bufst(3) + ffz*zz
!                bufst(4) = bufst(4) + ffy*zz
!                bufst(5) = bufst(5) + ffz*xx
!                bufst(6) = bufst(6) + ffx*yy
      end if
   end do
end do
else
!--- idouble_grid_method == 3 (cluster geometry )------------------
do i_node  = 1, nion
   i   = iatmpt(i_node)
   it1 = iatoit(i)
   do j  = 1, i - 1

      it2 = iatoit(j)
      zv2 = zv(it1)*zv(it2)*2.d0
      q1 = ratm(1,j) - ratm(1,i)
      q2 = ratm(2,j) - ratm(2,i)
      q3 = ratm(3,j) - ratm(3,i)
      xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
      yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
      zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
      r2 = xx*xx + yy*yy + zz*zz
      rr = sqrt( r2 )
      r3 = r2*rr

      tecolmb = zv2/rr
      ecolmb = ecolmb + tecolmb
      atom_ecolmb(i) = atom_ecolmb(i) + tecolmb
      atom_ecolmb(j) = atom_ecolmb(j) + tecolmb

      r3 = zv2/r3
      ffx =  - xx*r3
      ffy =  - yy*r3
      ffz =  - zz*r3
      fclm(1,i) = fclm(1,i) + ffx
      fclm(2,i) = fclm(2,i) + ffy
      fclm(3,i) = fclm(3,i) + ffz
      fclm(1,j) = fclm(1,j) - ffx
      fclm(2,j) = fclm(2,j) - ffy
      fclm(3,j) = fclm(3,j) - ffz

   end do
end do
end if
!--- end of direct coulomb interation in real space --------------------
end if
call gdsum(ecolmb,1,dbuf1r)
do i = 1, natom
   atom_ecolmb(i) = atom_ecolmb(i) * 0.5d0
end do
if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '    direct Coulomb in real space ',  &
&                         ': cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '    direct Coulomb in real space ',  &
&                         ': cpu-time :', ct - ct0
    ct0 = ct
end if

if( .not.lclust .and. idouble_grid_method /= 3 .or.  &
&   .not.lclust .and. .not.ldouble_grid_recip ) then
!--- reciprocal space contribution
    nelv = 0.d0
    do it = 1, ntype
       nelv = nelv + (nhk2(it) - nhk1(it) + 1)*nint(zv(it))
    end do

    call frcrec( fclm, strclm, ntype, nhk1, nhk2, zv, nelv,  &
& ecolmb, atom_ecolmb, iatoit, iatmpt, natom, nion,  &
& dbuf, dbufr, lcstress )

!    if( ldouble_grid_recip ) then
!        call frcrec_dg( nfile, myid, nodes,  &
!& fclm, ntype, nhk1, nhk2, zv, nelv,  &
!& ecolmb, atom_ecolmb, iatoit, iatmpt, natom, nion, dbuf, dbufr,  &
!& idouble_grid_method, lvacuum, hcell  )
!    end if

    if( ltimecnt ) then
        call gsync
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) '            in reciproccal space ',  &
&                             ':          :', ct - ct0
        if(loutfile(2)) write(nfile(2),*) '            in reciproccal space ',  &
&                             ':          :', ct - ct0
        ct0 = ct
    end if
end if

do i = 1, nhk2(ntype)
   dbuf(3*i-2) = fclm(1,i)
   dbuf(3*i-1) = fclm(2,i)
   dbuf(3*i-0) = fclm(3,i)
end do

call gdsum(dbuf,3*nhk2(ntype),dbufr)

do i = 1, nhk2(ntype)
   fclm(1,i) = dbuf(3*i-2)
   fclm(2,i) = dbuf(3*i-1)
   fclm(3,i) = dbuf(3*i-0)
end do


return
end




subroutine flocal( nfile, myid, nodes,                             &
& floc, strloc, lclust, nd1v, ntype, nhk1, nhk2,  &
& nhk1_nod, nhk2_nod, zv, nel, ratm, iatoit, iatmpt, nion,  &
& natom, llking, mx1, tbflc, tbflca, dltflc, rmxflc,  &
& hcell, lorthrhmbc, lhfull, rdel, rdelv, rho, noddatx,  &
& ltimecnt, lcstress,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz,  &
& tmpk, tmpl, tmpm, tmpn, dbuf, dbufr )
!-----------------------------------------------------------------------
!   floc : force by local pseudopotential
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8,  dimension(3,natom) :: floc
real*8,  dimension(3,3)  :: strloc
logical :: lclust
integer, dimension(3) :: nd1v
real*8,  dimension(3,3) :: hcell
logical :: lorthrhmbc
logical :: lhfull
integer :: ntype
real*8,  dimension(ntype) :: zv
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(ntype) :: nhk1_nod, nhk2_nod
real*8,  dimension(3,natom)   :: ratm
integer, dimension(natom)     :: iatoit
integer, dimension(nion)  :: iatmpt
integer :: nion
logical, dimension(ntype) :: llking
real*8,  dimension(0:mx1,ntype) :: tbflc, tbflca
real*8,  dimension(ntype)       :: dltflc, rmxflc
real*8,  dimension(3*natom) :: dbuf, dbufr

!--- for bulk calculation ----------------------------------------
real*8  :: sbl1(8) = (/ 1.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0 /)
real*8  :: sbl2(8) = (/ 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0 /)
real*8  :: sbl3(8) = (/ 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0 /)
dimension      bufst(6), bufstr(6)
dimension      mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
!-----------------------------------------------------------------

dimension      rho(noddatx)
dimension      rdel(3)
dimension      mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension tmpk(*), tmpl(*), tmpm(*), tmpn(*)
logical   ltimecnt, lcstress
save sbl1, sbl2, sbl3


if( ltimecnt ) then
    ct0 = timecnt()
end if
if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
#ifdef VECTOR
   do it = 1, ntype
      zvit2  = 2.d0*zv(it)
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
#endif
do m1 = 1, mshnod
   ix = mshnx(m1)
   ix = mshx1 + ix - 1
   iy = mshny(m1)
   iy = mshy1 + iy - 1
   iz = mshnz(m1)
   iz = mshz1 + iz - 1
   xm = dble(ix-1)*rdel(1)
   ym = dble(iy-1)*rdel(2)
   zm = dble(iz-1)*rdel(3)
#ifndef VECTOR
   do it = 1, ntype
      zvit2  = 2.d0*zv(it)
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
#endif
      xx = xm - ratm(1,i)
      yy = ym - ratm(2,i)
      zz = zm - ratm(3,i)
      r2 = xx*xx + yy*yy + zz*zz
      if( r2.le.rmxlc2 ) then
          m = r2/(2.d0*dltflc(it))
          m = 2*m
          d = 0.5d0*( r2/dltflc(it) - dble(m) )
          ffloc = d*( (d-1.d0)*tbflca(m,it) + tbflc(m+2,it)  &
&                          - tbflc(m,it) ) + tbflc(m,it)
        else
          r = sqrt(r2)
          ffloc = zvit2/(r2*r)
      end if
      ffloc = ffloc * rho(m1)
      floc(1,i) = floc(1,i) + ffloc*xx
      floc(2,i) = floc(2,i) + ffloc*yy
      floc(3,i) = floc(3,i) + ffloc*zz
   end do
   end do
end do
else
!--- for bulk calculation ----------------------------------------------
do i = 1, 6
   bufst(i) = 0.d0
end do
if( lhfull ) then
    ihroop = 8
    dhalf  = 0.d0
  else
    ihroop = 1
    dhalf  = 0.5d0
end if
if( lcstress ) then
!--- when stress is needed ------------------
#ifdef VECTOR
   typedo: do it = 1, ntype
   typeif: if( llking(it) ) then
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
do m1 = 1, mshnod
   ix = mshnx(m1)
   ix = mshx1 + ix - 1
   iy = mshny(m1)
   iy = mshy1 + iy - 1
   iz = mshnz(m1)
   iz = mshz1 + iz - 1
   xm = dble(ix-1)/dble(nd1v(1))
   ym = dble(iy-1)/dble(nd1v(2))
   zm = dble(iz-1)/dble(nd1v(3))
#ifndef VECTOR
   typedo: do it = 1, ntype
   typeif: if( llking(it) ) then
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
      q1 = xm - ratm(1,i)
      q2 = ym - ratm(2,i)
      q3 = zm - ratm(3,i)
      if( abs(q1).ge.dhalf ) q1 = q1 - sign(sbl1(lloop),q1)
      if( abs(q2).ge.dhalf ) q2 = q2 - sign(sbl2(lloop),q2)
      if( abs(q3).ge.dhalf ) q3 = q3 - sign(sbl3(lloop),q3)
      if( lorthrhmbc ) then
          xx = hcell(1,1)*q1
          yy = hcell(2,2)*q2
          zz = hcell(3,3)*q3
       else
          xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
          yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
          zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
      end if
      r2 = xx*xx + yy*yy + zz*zz
!#ifdef VECTOR
!            tmpk(m1) = xx
!            tmpl(m1) = yy
!            tmpm(m1) = zz
!            tmpn(m1) = r2
!      end do
!      do m1 = 1, mshnod
!            xx = tmpk(m1)
!            yy = tmpl(m1)
!            zz = tmpm(m1)
!            r2 = tmpn(m1)
!#endif
      if( r2.le.rmxlc2 ) then
          m = r2/(2.d0*dltflc(it))
          m = 2*m
          d = 0.5d0*( r2/dltflc(it) - dble(m) )
          ffloc = d*( (d-1.d0)*tbflca(m,it) + tbflc(m+2,it)  &
&                            - tbflc(m,it) ) + tbflc(m,it)
          ffloc = ffloc * rho(m1)
          fflocx = ffloc*xx
          fflocy = ffloc*yy
          fflocz = ffloc*zz
          floc(1,i) = floc(1,i) + fflocx
          floc(2,i) = floc(2,i) + fflocy
          floc(3,i) = floc(3,i) + fflocz

          bufst(1) = bufst(1) + fflocx*xx
          bufst(2) = bufst(2) + fflocy*yy
          bufst(3) = bufst(3) + fflocz*zz
          bufst(4) = bufst(4) + fflocy*zz
          bufst(5) = bufst(5) + fflocz*xx
          bufst(6) = bufst(6) + fflocx*yy
      end if
#ifndef VECTOR
   end do
   end do
   end if typeif
   end do typedo
#endif
end do
#ifdef VECTOR
   end do
   end do
   end if typeif
   end do typedo
#endif
call gdsum(bufst,6,bufstr)
strloc(1,1) = bufst(1) * rdelv
strloc(2,2) = bufst(2) * rdelv
strloc(3,3) = bufst(3) * rdelv
strloc(2,3) = bufst(4) * rdelv
strloc(3,1) = bufst(5) * rdelv
strloc(1,2) = bufst(6) * rdelv
else
!--- when stress is not needed ------------------
#ifdef VECTOR
   typedo2: do it = 1, ntype
   typeif2: if( llking(it) ) then
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
do m1 = 1, mshnod
   ix = mshnx(m1)
   ix = mshx1 + ix - 1
   iy = mshny(m1)
   iy = mshy1 + iy - 1
   iz = mshnz(m1)
   iz = mshz1 + iz - 1
   xm = dble(ix-1)/dble(nd1v(1))
   ym = dble(iy-1)/dble(nd1v(2))
   zm = dble(iz-1)/dble(nd1v(3))
#ifndef VECTOR
   typedo2: do it = 1, ntype
   typeif2: if( llking(it) ) then
      rmxlc2 = rmxflc(it)*rmxflc(it)
   do i = nhk1(it), nhk2(it)
   do lloop = 1, ihroop
#endif
      q1 = xm - ratm(1,i)
      q2 = ym - ratm(2,i)
      q3 = zm - ratm(3,i)
      if( abs(q1).ge.dhalf ) q1 = q1 - sign(sbl1(lloop),q1)
      if( abs(q2).ge.dhalf ) q2 = q2 - sign(sbl2(lloop),q2)
      if( abs(q3).ge.dhalf ) q3 = q3 - sign(sbl3(lloop),q3)
      if( lorthrhmbc ) then
          xx = hcell(1,1)*q1
          yy = hcell(2,2)*q2
          zz = hcell(3,3)*q3
       else
          xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
          yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
          zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
      end if
      r2 = xx*xx + yy*yy + zz*zz
!#ifdef VECTOR
!            tmpk(m1) = xx
!            tmpl(m1) = yy
!            tmpm(m1) = zz
!            tmpn(m1) = r2
!      end do
!      do m1 = 1, mshnod
!            xx = tmpk(m1)
!            yy = tmpl(m1)
!            zz = tmpm(m1)
!            r2 = tmpn(m1)
!#endif
      if( r2.le.rmxlc2 ) then
          m = r2/(2.d0*dltflc(it))
          m = 2*m
          d = 0.5d0*( r2/dltflc(it) - dble(m) )
          ffloc = d*( (d-1.d0)*tbflca(m,it) + tbflc(m+2,it)  &
&                            - tbflc(m,it) ) + tbflc(m,it)
          ffloc = ffloc * rho(m1)
          fflocx = ffloc*xx
          fflocy = ffloc*yy
          fflocz = ffloc*zz
          floc(1,i) = floc(1,i) + fflocx
          floc(2,i) = floc(2,i) + fflocy
          floc(3,i) = floc(3,i) + fflocz

!                bufst(1) = bufst(1) + fflocx*xx
!                bufst(2) = bufst(2) + fflocy*yy
!                bufst(3) = bufst(3) + fflocz*zz
!                bufst(4) = bufst(4) + fflocy*zz
!                bufst(5) = bufst(5) + fflocz*xx
!                bufst(6) = bufst(6) + fflocx*yy
      end if
#ifndef VECTOR
   end do
   end do
   end if typeif2
   end do typedo2
#endif
end do
#ifdef VECTOR
   end do
   end do
   end if typeif2
   end do typedo2
#endif
end if
!--- end of local pseudopotential in real space  -----------------------
end if
if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '          local pp in real space ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '          local pp in real space ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if

if( .not.lclust ) then
!--- reciprocal space contribution
    nelv = 0
    nelvv = 0
    do it = 1, ntype
       nelvv = nelvv + nint(zv(it))*( nhk2(it) - nhk1(it) + 1 )
       if( llking(it) ) then
           nelv = nelv + nint(zv(it))*( nhk2(it) - nhk1(it) + 1 )
       end if
    end do
    call refloc( nfile, myid, nodes,                               &
& floc, strloc, ntype, nhk1_nod, nhk2_nod, zv, nelvv, nelv, iatoit,  &
& iatmpt, natom, nion, llking, rdelv, lcstress, rho, noddatx,  &
& dbuf, dbufr, mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, mshnod )
    if( ltimecnt ) then
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) '            in reciproccal space ',  &
&                             ':          :', ct - ct0
        if(loutfile(2)) write(nfile(2),*) '            in reciproccal space ',  &
&                             ':          :', ct - ct0
        ct0 = ct
    end if
end if

do i = 1, nhk2(ntype)
   floc(1,i) = floc(1,i) * rdelv
   floc(2,i) = floc(2,i) * rdelv
   floc(3,i) = floc(3,i) * rdelv
end do


return
end




module param_dftd
!-----------------------------------------------------------------------
! type declaration and initialization of variables for DFT-D
!
!  DFT-D second version : S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799
!-----------------------------------------------------------------------
implicit none

integer :: init
!c VDW radii are in ang. determined by
!c atomic ROHF/TZV calculations and then taken as the radius of the 0.01 density contour
!c (scaling done below)
!c data obatained using the 0.005 contour and scale=1.0 gave too large differences
!c between the atoms i.e 1.1 1.5 1.42 1.36 1.3 for H,C,N,O,F
real*8, dimension(103) :: vander =                               &
!c H, He
&   (/ 0.91d0,0.92d0, &
!c Li-Ne
&      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
!c Na-Ar
&      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0, &
!c K, Ca
&      1.35d0,1.34d0, &
!c Sc-Zn
&      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
&      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
!c Ga-Kr
&      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0, &
!c Rb, Sr
&      1.48d0,1.46d0, &
!c Y-Cd
&      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
&      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
!c In, Sn, Sb, Te, I, Xe
&      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0, &
&      (1.d+99,init=1,49) /)

!c C6 parameters in J mol^-1 nm^6
real*8, dimension(103) :: c6 =                                 &
& (/ 0.14d0, 0.08d0, 1.61d0, 1.61d0, 3.13d0, 1.75d0, 1.23d0, 0.70d0, 0.75d0, 0.63d0, &
& 5.71d0, 5.71d0,10.79d0, 9.23d0, 7.84d0, 5.57d0, 5.07d0, 4.61d0,10.80d0,10.80d0, &
& 10.80d0,10.80d0,10.80d0,10.80d0,10.80d0,10.80d0,10.80d0,10.80d0,10.80d0,10.80d0, &
& 16.99d0,17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0,24.67d0,24.67d0, &
& 24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,37.32d0,38.71d0, &
& 38.44d0,31.74d0,31.50d0,29.99d0, (1.d-14,init=1,49) /)


real*8  :: alp = 20.0d0
real*8  :: scalerad = 1.10d0   ! vander <- vander * scalerad
real*8  :: scalec6             ! DFT-functional depenent scale factor

real*8,  allocatable, dimension(:,:) :: c6ij, rrij

!-----energy, force, and stress
real*8 :: edisp                                 ! dispersion energy
real*8, allocatable, dimension(:,:) :: fdftd    ! force
real*8, dimension(3,3) :: strdftd = 0.d0        ! stress
real*8, allocatable, dimension(:) :: atom_dftd  ! for EDA

!-----cutoff distance
real*8  :: cutoff

save

end module




subroutine setdftd( nfile, myid, nodes, &
& jgga, ntype, nhk1, nhk2, zatom, evdj, hrdev, audang, avogad )
!-----------------------------------------------------------------------
! setup DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use outfile
use param_dftd
implicit none
integer :: nfile(*), myid, nodes
integer :: jgga
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(ntype) :: zatom
real*8  :: evdj, hrdev, audang, avogad

!-----declare local variables
integer :: it, jt, iz, jz, natom
real*8  :: const
real*8  :: the_mem
integer :: status
integer :: natomx = 0
logical :: licall = .true.
save licall, natomx


if( licall ) then
    !------allocate memory
    allocate( c6ij(ntype,ntype), rrij(ntype,ntype), stat=status )

if( jgga == 2 ) then
!c PBE
    scalec6 = 0.75d0
else if( jgga == 3 ) then
!c RPBE
    scalec6 = 0.75d0  ! Caution!!  This is not optimized value!
else if( jgga == 4 ) then
!c revPBE
    scalec6 = 0.75d0  ! Caution!!  This is not optimized value!
else
    call fstop( nfile, myid, nodes, 'not yet supported in DFT-D' )
end if
const  = scalec6 /evdj/hrdev*(10d0/audang)**6/avogad * 2.d0 ! in [Ryd.]
c6     = c6     * const
vander = vander * scalerad / audang

if(loutfile(1)) write(nfile(1),*) 'Parameters for DFT-D'
if(loutfile(2)) write(nfile(2),*) 'Parameters for DFT-D'
do it =  1, ntype
   iz = nint(zatom(it))
   if(loutfile(1)) write(nfile(1),*) ' it=', it, ' c6,R0=', c6(iz), vander(iz)
   if(loutfile(2)) write(nfile(2),*) ' it=', it, ' c6,R0=', c6(iz), vander(iz)
end do

do it =  1, ntype
   iz = nint(zatom(it))
do jt = it, ntype
   jz = nint(zatom(jt))

   c6ij(it,jt) = sqrt(c6(iz)*c6(jz))
   rrij(it,jt) = vander(iz) + vander(jz)

   c6ij(jt,it) = c6ij(it,jt)
   rrij(jt,it) = rrij(it,jt)

end do
end do

cutoff = 0.d0
do it =  1, ntype
   cutoff = max( cutoff, rrij(it,it) )
end do
cutoff = cutoff * 4.d0
cutoff = cutoff * cutoff

    licall = .false.
end if


natom = nhk2(ntype)
if( natom > natomx ) then
    !-----if already allocated, deallocate arrays
    if( allocated(fdftd) ) then
        the_mem = 8.d0 * ( size(fdftd) + size(atom_dftd) )
        deallocate( fdftd, atom_dftd, stat=status )

        !------error trap
        call check_dealloc_accum( nfile, myid, nodes,  &
& status, the_mem, 'setdftd', .true. )
    end if


    natomx = natom
    !------allocate memory
    allocate( fdftd(3,natomx), atom_dftd(natomx), stat=status )

    the_mem = 8.d0 * ( size(fdftd) + size(atom_dftd) )

    !------error trap
    call check_alloc_accum( nfile, myid, nodes,  &
& status, the_mem, 'setdftd', .true. )
end if


return
end




subroutine force_dftd( nfile, myid, nodes, &
& lclust, ntype, nhk1, nhk2, ratm,  &
& iatoit, iatmpt, nion, natom, hcell,  &
& dbuf, dbufr, ltimecnt, lcstress,  &
& ldouble_grid_recip, idouble_grid_method, lvacuum )
!-----------------------------------------------------------------------
! forces by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use outfile
use param_dftd
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: lclust
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
integer :: natom
real*8,  dimension(3,natom)   :: ratm
integer, dimension(natom)     :: iatoit
integer :: nion
integer, dimension(nion)  :: iatmpt
real*8,  dimension(3,3) :: hcell
real*8,  dimension(3,natom) :: dbuf, dbufr
logical :: ltimecnt, lcstress, ldouble_grid_recip
integer :: idouble_grid_method
logical, dimension(3) :: lvacuum

!-----declare local variables
integer :: i_node, i, j, it1, it2
integer :: ihxx, ihyy, ihzz, ihx, ihy, ihz
real*8  :: xx, yy, zz, r2, r3, r6, rr, expd, rbr, fdmp, eij, ffx, ffy, ffz
real*8  :: q1, q2, q3, qq1, qq2, qq3, dbuf1r, timecnt
real*8,  dimension(6) :: bufst, bufstr
real*8  :: ct0, ct


if( ltimecnt ) then
    ct0 = timecnt()
end if

edisp = 0.d0
fdftd = 0.d0
strdftd = 0.d0
atom_dftd = 0.d0

if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
      do i_node  = 1, nion
         i   = iatmpt(i_node)
         it1 = iatoit(i)
         do j  = 1, i - 1
            it2 = iatoit(j)
            xx = ratm(1,j) - ratm(1,i)
            yy = ratm(2,j) - ratm(2,i)
            zz = ratm(3,j) - ratm(3,i)
            r2 = xx*xx + yy*yy + zz*zz
            if( r2 < cutoff ) then
                rr = sqrt( r2 )
                rbr = rr/rrij(it1,it2)
                expd = exp(-alp*(rbr-1.d0))
                fdmp = 1.d0/(1.d0+expd)
                r3 = r2*rr
                r6 = r3*r3
                eij = -c6ij(it1,it2)/r6 * fdmp

                edisp = edisp + eij
                atom_dftd(i) = atom_dftd(i) + eij
                atom_dftd(j) = atom_dftd(j) + eij

                eij = eij/r2 * (6.d0 - alp*rbr*expd*fdmp)
                ffx =  - xx*eij
                ffy =  - yy*eij
                ffz =  - zz*eij
                fdftd(1,i) = fdftd(1,i) + ffx
                fdftd(2,i) = fdftd(2,i) + ffy
                fdftd(3,i) = fdftd(3,i) + ffz
                fdftd(1,j) = fdftd(1,j) - ffx
                fdftd(2,j) = fdftd(2,j) - ffy
                fdftd(3,j) = fdftd(3,j) - ffz
            end if
         end do
      end do
else
!--- for bulk calculation ----------------------------------------------
      ihxx = 1
      if( lvacuum(1) ) ihxx = 0
      ihyy = 1
      if( lvacuum(2) ) ihyy = 0
      ihzz = 1
      if( lvacuum(3) ) ihzz = 0
      bufst = 0.d0
      do i_node  = 1, nion
         i   = iatmpt(i_node)
         it1 = iatoit(i)
         do j  = 1, i - 1
            it2 = iatoit(j)
            qq1 = ratm(1,j) - ratm(1,i)
            qq2 = ratm(2,j) - ratm(2,i)
            qq3 = ratm(3,j) - ratm(3,i)
            do ihx = 0, ihxx
               q1 = qq1 - sign(dble(ihx),qq1)
            do ihy = 0, ihyy
               q2 = qq2 - sign(dble(ihy),qq2)
            do ihz = 0, ihzz
               q3 = qq3 - sign(dble(ihz),qq3)
!            if( abs(q1) > 0.5d0 ) q1 = q1 - sign(dble(ihx),q1)
!            if( abs(q2) > 0.5d0 ) q2 = q2 - sign(dble(ihy),q2)
!            if( abs(q3) > 0.5d0 ) q3 = q3 - sign(dble(ihz),q3)
!            xx = hcell(1,1)*q1
!            yy = hcell(2,2)*q2
!            zz = hcell(3,3)*q3
            xx = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
            yy = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
            zz = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
            r2 = xx*xx + yy*yy + zz*zz
            if( r2 < cutoff ) then
                rr = sqrt( r2 )
                rbr = rr/rrij(it1,it2)
                expd = exp(-alp*(rbr-1.d0))
                fdmp = 1.d0/(1.d0+expd)
                r3 = r2*rr
                r6 = r3*r3
                eij = -c6ij(it1,it2)/r6 * fdmp

                edisp = edisp + eij
                atom_dftd(i) = atom_dftd(i) + eij
                atom_dftd(j) = atom_dftd(j) + eij

                eij = eij/r2 * (6.d0 - alp*rbr*expd*fdmp)
                ffx =  - xx*eij
                ffy =  - yy*eij
                ffz =  - zz*eij
                fdftd(1,i) = fdftd(1,i) + ffx
                fdftd(2,i) = fdftd(2,i) + ffy
                fdftd(3,i) = fdftd(3,i) + ffz
                fdftd(1,j) = fdftd(1,j) - ffx
                fdftd(2,j) = fdftd(2,j) - ffy
                fdftd(3,j) = fdftd(3,j) - ffz

            if( lcstress ) then
                bufst(1) = bufst(1) + ffx*xx
                bufst(2) = bufst(2) + ffy*yy
                bufst(3) = bufst(3) + ffz*zz
                bufst(4) = bufst(4) + ffy*zz
                bufst(5) = bufst(5) + ffz*xx
                bufst(6) = bufst(6) + ffx*yy
            end if

            end if
            end do
            end do
            end do
         end do
      end do
      if( lcstress ) then
          call gdsum(bufst,6,bufstr)
          strdftd(1,1) = bufst(1)
          strdftd(2,2) = bufst(2)
          strdftd(3,3) = bufst(3)
          strdftd(2,3) = bufst(4)
          strdftd(3,1) = bufst(5)
          strdftd(1,2) = bufst(6)
      end if
!--- end of calculation in real space --------------------
end if
call gdsum(edisp,1,dbuf1r)
do i = 1, natom
   atom_dftd(i) = atom_dftd(i) * 0.5d0
end do
if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '      DFT-D empirical correction ',  &
&                         ': cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '      DFT-D empirical correction ',  &
&                         ': cpu-time :', ct - ct0
    ct0 = ct
end if

dbuf = fdftd
call gdsum(dbuf,3*nhk2(ntype),dbufr)
fdftd = dbuf


return
end




subroutine add_force_dftd( nfile, myid, nodes, frc, natom )
!-----------------------------------------------------------------------
! forces by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use param_dftd
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: natom
real*8,  dimension(3,natom)   :: frc

frc = frc + fdftd

return
end




subroutine add_energy_dftd( nfile, myid, nodes, sume, etot, edisp_ )
!-----------------------------------------------------------------------
! energy by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use param_dftd
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: sume, etot, edisp_

sume = sume + edisp
etot = etot + edisp
edisp_ = edisp

return
end




subroutine out_energy_dftd( nfile, myid, nodes, elvlt, ii )
!-----------------------------------------------------------------------
! energy by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use param_dftd
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: elvlt
integer :: ii

write(nfile(ii),1120) edisp, edisp*elvlt
1120 format('  * DFT-D empirical E.= ',F15.6,' ( Ryd. )',F13.4,' (eV)')

return
end




subroutine out_force_dftd( nfile, myid, nodes, ii, i )
!-----------------------------------------------------------------------
! energy by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use param_dftd
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ii, i
integer :: ix

write(nfile(ii),'(6x,a12,3e19.10)') &
& ' (DFT-D)   :', ( fdftd(ix,i), ix = 1, 3)

return
end



subroutine stress_dftd( nfile, myid, nodes, &
& str, rvol, prau, ltimecnt )
!-----------------------------------------------------------------------
! energy by DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use outfile
use param_dftd
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8,  dimension(3,3) :: str
real*8  :: rvol, prau
logical :: ltimecnt
!-----declare local variables
integer :: i, j, nf


strdftd = strdftd * (-1.d0) / rvol

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*) 'Stress by DFT-D in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strdftd(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if
721   format(5x,3e16.8)

str = str + strdftd


return
end




