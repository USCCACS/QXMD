



subroutine getmultg( lclust, multg, nd1v, nd2v, mulnd2 )
!-----------------------------------------------------------------------
!     get multigrid level
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
logical        lclust
dimension      nd1v(3)
dimension      mulnd2(multg)


if( lclust ) then
!--- for cluster calculation -------------------------------------------
!cc      multg = multgx
do i = 1, 3
   mul = 1
   nd1vi = nd1v(i)
   do
      if( nd1vi <= 2 ) exit
      mul = mul + 1
      nd1vi = ( nd1vi + 1 )/2
   end do
   multg = min( multg, mul )
end do

mulmsh = min( nd1v(1), nd1v(2), nd1v(3) )
do mul = 1, multg
   if( mulmsh.le.2 ) then
       multg = mul - 1
       return
   end if
   if( mulmsh.lt.2*nd2v+1 ) then
       mulnd2(mul) = (mulmsh-1)/2
       if( mul.eq.1 ) nd2v = mulnd2(mul)
     else
       mulnd2(mul) = nd2v
   end if
   mulmsh = ( mulmsh + 1 )/2
end do
else
!--- for bulk calculation ----------------------------------------------
!cc      multg = multgx
do i = 1, 3
   mul = 1
   nd1vi = nd1v(i)
   do
      if( mod(nd1vi,2) /= 0 .or. nd1vi <= 2 ) exit
      mul = mul + 1
      nd1vi = nd1vi/2
   end do
   multg = min( multg, mul )
end do

mulmsh = min( nd1v(1), nd1v(2), nd1v(3) )
do mul = 1, multg
   if( mulmsh.le.2 ) then
       multg = mul - 1
       return
   end if
   if( mulmsh.lt.2*nd2v+1 ) then
       mulnd2(mul) = (mulmsh-1)/2
       if( mul.eq.1 ) nd2v = mulnd2(mul)
     else
       mulnd2(mul) = nd2v
   end if
   mulmsh = mulmsh/2
end do
!-----------------------------------------------------------------------
end if


return
end




subroutine getsmultg( lclust, multg, nd1v, muls, nd1vs,  &
&                     mulnd2, kulnd2, msrhmx, msrhmy, msrhmz )
!-----------------------------------------------------------------------
!     get multigrid level for serial calculation
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
logical :: lclust
integer :: multg
integer :: muls
integer, dimension(3) :: nd1v
integer, dimension(3) :: nd1vs
integer, dimension(multg) :: mulnd2
integer, dimension(multg) :: kulnd2
integer :: msrhmx, msrhmy, msrhmz


mulms1 = nd1v(1)
mulms2 = nd1v(2)
mulms3 = nd1v(3)
muls   = multg + 1
do mul = 1, multg
   if( mulms1.le.msrhmx .and. mulms2.le.msrhmy .and.  &
&      mulms3.le.msrhmz ) then
       nd1vs(1) = mulms1
       nd1vs(2) = mulms2
       nd1vs(3) = mulms3
       muls     = mul
       exit
   end if

   if( lclust ) then
!--- for cluster calculation
       mulms1 = ( mulms1 + 1 )/2
       mulms2 = ( mulms2 + 1 )/2
       mulms3 = ( mulms3 + 1 )/2
     else
!--- for bulk calculation
       mulms1 = mulms1/2
       mulms2 = mulms2/2
       mulms3 = mulms3/2
   end if
end do

do mul = muls, multg
   kulnd2(mul-muls+1) = mulnd2(mul)
end do


return
end




subroutine setmsh( nfile, myid, nodes, mytid, nodess,  &
& myx, myy, myz, nn, myparity, npx, npy, npz, nproc, lsphere,  &
& lclust, multg, muls, nd1v, rdel, rdelv, rmax, cofmas, lvacuum,  &
& mshnod, mshnew, mshxyz, mulpit, mulx1, mulx2, muly1, muly2,  &
& mulz1, mulz2, mshnx, mshny, mshnz, mshx1, mshy1, mshz1,  &
& mshx, mshy, mshz, mulpms, mulpem, ismx, ismy, ismz,  &
& ismshx, ismshy, ismshz, mulpsx, mulpsy, mulpsz,  &
& mulyzs, mulzxs, mulxys, irmx, irmy, irmz, irmshx, irmshy, irmshz,  &
& mulprx, mulpry, mulprz, mulyzr, mulzxr, mulxyr,  &
& ismfpx, ismfpy, ismfpz, irmfpx, irmfpy, irmfpz, ncomct, mulnd2,  &
& iultgx, iulexf, iuldatx, iulnpx, iulnpy, iulnpz,  &
& iulyzx, iulzxx, iulxyx, iulext,  &
& nmm, nmmcnt, nmmdsp, nmmrc1, nmmrc2,  &
& meshx1, meshx, meshy1, meshy, meshz1, meshz,  &
& jsmshx, jsmshy, jsmshz, nodyzh, nodzxh, nodxyh,  &
& idbuf, idbufr, mshndv, nhit1, mulgmx, mulgmy, mulgmz )
!-----------------------------------------------------------------------
!   nmm : total number of mesh points
!   nxyz(ix,iy,iz)          : number of the mesh point
!   ( nx(i), ny(i), nz(i) ) : coordinate of the i-th mesh point
!-----------------------------------------------------------------------
!   ( mshx1, mshy1, mshz1)     : origin of my subdivision
!   ( mshnx(i), mshny(i), mshnz(i) ) : coordinate of the i-th mesh point
!                                      in my subdivision
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )

integer :: myid, nodes, mytid, nodess
integer, dimension(*) :: nfile
integer, dimension(6) :: nn
integer, dimension(3) :: myparity
logical :: lclust
logical :: lsphere

dimension rmax(3)
dimension nd1v(3)
dimension rdel(3)
dimension cofmas(3)
logical   lvacuum(3)

dimension mshnod(iultgx)
dimension mshnew(iultgx)
dimension mshxyz(iulexf)
dimension mulpit(iultgx), mulx1(iultgx), mulx2(iultgx),  &
&                         muly1(iultgx), muly2(iultgx),  &
&                         mulz1(iultgx), mulz2(iultgx)
dimension mshnx(iuldatx), mshny(iuldatx), mshnz(iuldatx)
dimension mshx1(iultgx), mshy1(iultgx), mshz1(iultgx),  &
&         mshx(iultgx), mshy(iultgx), mshz(iultgx)
dimension mulpms(iultgx)
dimension mulpem(iultgx)

dimension ismx(iulnpx), ismy(iulnpy), ismz(iulnpz),  &
&         ismshx(iulyzx), ismshy(iulzxx), ismshz(iulxyx)
dimension mulpsx(iultgx), mulpsy(iultgx), mulpsz(iultgx),  &
&         mulyzs(iultgx), mulzxs(iultgx), mulxys(iultgx)
dimension irmx(iulnpx), irmy(iulnpy), irmz(iulnpz),  &
&         irmshx(iulyzx), irmshy(iulzxx), irmshz(iulxyx)
dimension jsmshx(3,nodyzh,2), jsmshy(3,nodzxh,2),  &
&         jsmshz(3,nodxyh,2)
dimension mulprx(iultgx), mulpry(iultgx), mulprz(iultgx),  &
&         mulyzr(iultgx), mulzxr(iultgx), mulxyr(iultgx)
dimension ismfpx(iulnpx), ismfpy(iulnpy), ismfpz(iulnpz)
dimension irmfpx(iulnpx), irmfpy(iulnpy), irmfpz(iulnpz)
dimension ncomct(12,iultgx)
dimension mulnd2(iultgx)
dimension idbuf(mshndv), idbufr(mshndv)

dimension nmmcnt(nproc), nmmdsp(nproc)
dimension nmmrc1(0:nproc-1), nmmrc2(0:nproc-1)
dimension meshx1(0:npx-1), meshx(0:npx-1),  &
&         meshy1(0:npy-1), meshy(0:npy-1),  &
&         meshz1(0:npz-1), meshz(0:npz-1)
dimension mulgmx(0:npx-1), mulgmy(0:npy-1), mulgmz(0:npz-1)

real*8  :: mshndm, mulpmul, ibuf1


!--- error trap
if( .not.lclust ) then
    if( mod(nd1v(1),2**(multg-1)).ne.0 .or.  &
&       mod(nd1v(2),2**(multg-1)).ne.0 .or.  &
&       mod(nd1v(3),2**(multg-1)).ne.0      ) then
        do i = 1, 2
        if( loutfile(i) ) then
        write(nfile(i),*) 'Multigrid level is inconsistent',  &
&                         ' with No. of mesh points.'
        write(nfile(i),*) nd1v(1), nd1v(2), nd1v(3), multg
        end if
        end do
!       -----Finalize the parallel environment
!        call end_parallel(ierr)
!        stop
        call fstop( nfile, myid, nodes, '*** error in setmsh' )
    end if
end if


mshxx = nd1v(1)
call divnod( mshxx, npx, myx, mshx1(1), mshx2, mshx(1) )

mshyy = nd1v(2)
call divnod( mshyy, npy, myy, mshy1(1), mshy2, mshy(1) )

mshzz = nd1v(3)
call divnod( mshzz, npz, myz, mshz1(1), mshz2, mshz(1) )


if( lclust ) then
    rmaxr1 = 1.d0/rmax(1)
    rmaxr2 = 1.d0/rmax(2)
    rmaxr3 = 1.d0/rmax(3)
end if
do mul = 1, multg
!if(loutfile(1)) write(nfile(1),'(a23,i3,a3,100a1)')  &
!&     ' ---- multigrid level (', mul, ' ) ',('-',i=1,40)
!if(loutfile(2)) write(nfile(2),'(a23,i3,a3,100a1)')  &
!&     ' ---- multigrid level (', mul, ' ) ',('-',i=1,40)
!=== start setting multigrid ===========================================
if( mul.eq.1 ) then
    mulpit(mul) = 1
    mulpms(mul) = 1
    mulpem(mul) = 1
  else
    mulpit(mul) = mulpit(mul-1) + incxyz
    mulpms(mul) = mulpms(mul-1) + mshnod(mul-1)
    mulpem(mul) = mulpem(mul-1) + mextsh
!--- x-grid
    mshx(mul)  = 0
    if( mod(mshx1(mul-1)+1,2).eq.0 ) then
        mshx1(mul) = ( mshx1(mul-1) + 1 )/2
      else
        mshx1(mul) = ( mshx1(mul-1) + 2 )/2
    end if
    do i = 1, mshx(mul-1)
       if( mod(mshx1(mul-1)+i,2).eq.0 ) then
           mshx(mul) = mshx(mul) + 1
       end if
    end do
!--- y-grid
    mshy(mul)  = 0
    if( mod(mshy1(mul-1)+1,2).eq.0 ) then
        mshy1(mul) = ( mshy1(mul-1) + 1 )/2
      else
        mshy1(mul) = ( mshy1(mul-1) + 2 )/2
    end if
    do i = 1, mshy(mul-1)
       if( mod(mshy1(mul-1)+i,2).eq.0 ) then
           mshy(mul) = mshy(mul) + 1
       end if
    end do
!--- z-grid
    mshz(mul)  = 0
    if( mod(mshz1(mul-1)+1,2).eq.0 ) then
        mshz1(mul) = ( mshz1(mul-1) + 1 )/2
      else
        mshz1(mul) = ( mshz1(mul-1) + 2 )/2
    end if
    do i = 1, mshz(mul-1)
       if( mod(mshz1(mul-1)+i,2).eq.0 ) then
           mshz(mul) = mshz(mul) + 1
       end if
    end do
end if

!if(loutfile(1)) write(nfile(1),*) '         order of FD :', mulnd2(mul)
!if(loutfile(2)) write(nfile(2),*) '         order of FD :', mulnd2(mul)

mulx1(mul)  = -mulnd2(mul) + 1
mulx2(mul)  = mshx(mul) + mulnd2(mul)
muly1(mul)  = -mulnd2(mul) + 1
muly2(mul)  = mshy(mul) + mulnd2(mul)
mulz1(mul)  = -mulnd2(mul) + 1
mulz2(mul)  = mshz(mul) + mulnd2(mul)

incrx = mulx2(mul) - mulx1(mul) + 1
incry = muly2(mul) - muly1(mul) + 1
incrz = mulz2(mul) - mulz1(mul) + 1
incrxy = incrx*incry
incxyz = incrx*incry*incrz
!--- error trap
mulmax = mulpit(mul) + incxyz - 1
call gimax(mulmax)
if( mulmax.gt.iulexf ) then
    if(loutfile(1)) write(nfile(1),*) 'mulmax.gt.iulexf', mulmax, iulexf
    if(loutfile(2)) write(nfile(2),*) 'mulmax.gt.iulexf', mulmax, iulexf
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in setmsh' )
end if



do ix = 1, incxyz
   mshxyz(mulpit(mul)+ix-1) = 0
end do


ierms1 = 0
if( lclust .and. lsphere ) then
    !--- set mesh for cluster calculation ----------------------------------
    !->      rmax2 = rmax * rmax * 1.00001d0
    rmax2 = 1.00001d0
    mshnod(mul)  = 0
    icmul1 = 2**(mul+(muls-1)-1)
    icmul2 = 0
    do j = 0, mul+(muls-1) - 2
       icmul2 = icmul2 - 2**j
    end do
    zdo: do iz = 1, mshz(mul)
       izz = mshz1(mul) + iz - 1
       izz = icmul1*izz + icmul2
       z  = ( dble(izz-1)*rdel(3) - cofmas(3) )*rmaxr3
       z2 = z*z
    ydo: do iy = 1, mshy(mul)
       iyy = mshy1(mul) + iy - 1
       iyy = icmul1*iyy + icmul2
       y  = ( dble(iyy-1)*rdel(2) - cofmas(2) )*rmaxr2
       y2 = y*y
    xdo: do ix = 1, mshx(mul)
       ixx = mshx1(mul) + ix - 1
       ixx = icmul1*ixx + icmul2
       x2 = ( dble(ixx-1)*rdel(1) - cofmas(1) )*rmaxr1
       x2 = x2*x2
       r2 = x2 + y2 + z2
       if( r2.lt.rmax2 ) then
           mshnod(mul) = mshnod(mul) + 1
           m2 = ix + mulnd2(mul) + (iy+mulnd2(mul)-1)*incrx  &
&                                + (iz+mulnd2(mul)-1)*incrxy
           mshxyz(mulpit(mul)+m2-1)    = mshnod(mul)
           !--- error trap
           if( mulpms(mul)+mshnod(mul)-1.gt.iuldatx ) then
               ierms1 = 1
!               if(loutfile(1)) 
                               write(nfile(1),*) 'error : iam =', myid,  &
&                   ' mulpms(mul)+mshnod(mul)-1.gt.iuldatx'
               if(loutfile(2)) write(nfile(2),*) 'error : iam =', myid,  &
&                   ' mulpms(mul)+mshnod(mul)-1.gt.iuldatx'
               exit zdo
           end if
           mshnx(mulpms(mul)+mshnod(mul)-1) = ix
           mshny(mulpms(mul)+mshnod(mul)-1) = iy
           mshnz(mulpms(mul)+mshnod(mul)-1) = iz
       end if
    end do xdo
    end do ydo
    end do zdo
else
    !--- set mesh for bulk calculation -------------------------------------
    mshnod(mul)  = 0
    zdo2: do iz = 1, mshz(mul)
    ydo2: do iy = 1, mshy(mul)
    xdo2: do ix = 1, mshx(mul)
           mshnod(mul) = mshnod(mul) + 1
           m2 = ix + mulnd2(mul) + (iy+mulnd2(mul)-1)*incrx  &
&                                + (iz+mulnd2(mul)-1)*incrxy
           mshxyz(mulpit(mul)+m2-1)    = mshnod(mul)
           !--- error trap
           if( mulpms(mul)+mshnod(mul)-1.gt.iuldatx ) then
               ierms1 = 1
!               if(loutfile(1)) 
                               write(nfile(1),*) 'error : iam =', myid,  &
&                   ' mulpms(mul)+mshnod(mul)-1.gt.iuldatx'
               if(loutfile(2)) write(nfile(2),*) 'error : iam =', myid,  &
&                   ' mulpms(mul)+mshnod(mul)-1.gt.iuldatx'
               exit zdo2
           end if
           mshnx(mulpms(mul)+mshnod(mul)-1) = ix
           mshny(mulpms(mul)+mshnod(mul)-1) = iy
           mshnz(mulpms(mul)+mshnod(mul)-1) = iz
    end do xdo2
    end do ydo2
    end do zdo2
    !--- end of setting mesh -----------------------------------------------
end if
!--- error trap
call gimax(ierms1)
if( ierms1.gt.0 ) then
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in setmsh' )
end if

mextsh = mshnod(mul) + 2*mulnd2(mul)*( mshx(mul)*mshy(mul)  &
&              + mshy(mul)*mshz(mul) + mshz(mul)*mshx(mul) )  &
&              + 4*(mshx(mul) + mshy(mul) + mshz(mul) ) + 8
!--- error trap
mulmax = mulpem(mul) + mextsh - 1
call gimax(mulmax)
if( mulmax.gt.iulext ) then
    if(loutfile(1)) write(nfile(1),*) 'mulmax.gt.iulext', mulmax, iulext
    if(loutfile(2)) write(nfile(2),*) 'mulmax.gt.iulext', mulmax, iulext
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in setmsh' )
end if


mshndm  = mshnod(mul)
mulpmul = mulpms(mul) - 1
if( nproc.gt.1 ) then
!    call gisum(mshndm,1,ibuf1)
!    call gisum(mulpmul,1,ibuf1)
    call gdsum(mshndm,1,ibuf1)
    call gdsum(mulpmul,1,ibuf1)
end if
mulpmul = mulpmul + 1d0
!if(loutfile(1)) write(nfile(1),*) '   No. of mesh point :', mshndm,  &
!&            ' (', mulpmul,' -',mulpmul+mshndm-1,' )'
!if(loutfile(2)) write(nfile(2),*) '   No. of mesh point :', mshndm,  &
!&            ' (', mulpmul,' -',mulpmul+mshndm-1,' )'
!if(loutfile(1)) write(nfile(1),'(a,f16.0,a,f16.0,a,f16.0,a)')  &
!& '   No. of mesh point :', mshndm,  &
!&            ' (', mulpmul,' -',mulpmul+mshndm-1.d0,' )'
!if(loutfile(2)) write(nfile(2),'(a,f16.0,a,f16.0,a,f16.0,a)')  &
!& '   No. of mesh point :', mshndm,  &
!&            ' (', mulpmul,' -',mulpmul+mshndm-1.d0,' )'

!--- set mesh points for communication ---------------------------------

if( mul.eq.1 ) then
    mulpsx(mul) = 1
    mulpsy(mul) = 1
    mulpsz(mul) = 1
    mulprx(mul) = 1
    mulpry(mul) = 1
    mulprz(mul) = 1
  else
    mulpsx(mul) = mulpsx(mul-1) + mulyzs(mul-1)*2
    mulpsy(mul) = mulpsy(mul-1) + mulzxs(mul-1)*2
    mulpsz(mul) = mulpsz(mul-1) + mulxys(mul-1)*2
    mulprx(mul) = mulprx(mul-1) + mulyzr(mul-1)*2
    mulpry(mul) = mulpry(mul-1) + mulzxr(mul-1)*2
    mulprz(mul) = mulprz(mul-1) + mulxyr(mul-1)*2
end if
mulyzs(mul) = mshy(mul)*mshz(mul)*mulnd2(mul)
mulzxs(mul) = mshz(mul)*mshx(mul)*mulnd2(mul) + 2*mshz(mul)
mulxys(mul) = mshx(mul)*mshy(mul)*mulnd2(mul)  &
&           + 2*(mshx(mul) + mshy(mul)) + 4
!--- error trap
ierr1 = 0
if( mulpsx(mul) + mulyzs(mul)*2 - 1 .gt. iulyzx .or.  &
&   mulpsy(mul) + mulzxs(mul)*2 - 1 .gt. iulzxx .or.  &
&   mulpsz(mul) + mulxys(mul)*2 - 1 .gt. iulxyx      ) then
    ierr1 = 1
!    if(loutfile(1)) 
                     write(nfile(1),*) 'error : iam =', myid,  &
&                   ' mulpsx(mul)+mulyzs(mul)*2-1.gt.iulyzx etc.',  &
&                     mulpsx(mul)+mulyzs(mul)*2-1,iulyzx,  &
&                     mulpsy(mul)+mulzxs(mul)*2-1,iulzxx,  &
&                     mulpsz(mul)+mulxys(mul)*2-1,iulxyx
    if(loutfile(2)) write(nfile(2),*) 'error : iam =', myid,  &
&                   ' mulpsx(mul)+mulyzs(mul)*2-1.gt.iulyzx etc.',  &
&                     mulpsx(mul)+mulyzs(mul)*2-1,iulyzx,  &
&                     mulpsy(mul)+mulzxs(mul)*2-1,iulzxx,  &
&                     mulpsz(mul)+mulxys(mul)*2-1,iulxyx
end if
call gimax(ierr1)
if( ierr1.gt.0 ) then
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in setmsh' )
end if
mulyzr(mul) = mulyzs(mul)
mulzxr(mul) = mulzxs(mul)
mulxyr(mul) = mulxys(mul)
mul2x = 2*(npx+1)*(mul-1) + 1
mul2y = 2*(npy+1)*(mul-1) + 1
mul2z = 2*(npz+1)*(mul-1) + 1
call bdset1( lvacuum,  &
&   mshnod(mul), mshnew(mul), mshxyz(mulpit(mul)),  &
&   mulx1(mul), mulx2(mul), muly1(mul), muly2(mul), mulz1(mul),  &
&   mulz2(mul), mshx1(mul), mshy1(mul), mshz1(mul), mshx(mul),  &
&   mshy(mul), mshz(mul), mulnd2(mul),  &
&   ismx(mul2x), ismy(mul2y), ismz(mul2z),  &
&   ismshx(mulpsx(mul)), ismshy(mulpsy(mul)), ismshz(mulpsz(mul)),  &
&   mulyzs(mul), mulzxs(mul), mulxys(mul),  &
&   irmx(mul2x), irmy(mul2y), irmz(mul2z),  &
&   irmshx(mulprx(mul)), irmshy(mulpry(mul)), irmshz(mulprz(mul)),  &
&   mulyzr(mul), mulzxr(mul), mulxyr(mul),  &
&   ismfpx(mul2x), ismfpy(mul2y), ismfpz(mul2z),  &
&   irmfpx(mul2x), irmfpy(mul2y), irmfpz(mul2z),  &
&   ncomct(1,mul),  &
&   jsmshx, jsmshy, jsmshz, nodyzh, nodzxh, nodxyh,  &
&   idbuf, idbufr, mshndv, myid, npx, npy, npz, nproc,  &
&   myx, myy, myz, nn, myparity, nhit1, nfile,  &
&   mulgmx, mulgmy, mulgmz )

end do
!=== end of setting multigrid ==========================================
200 continue



!##### parameters for global mesh ######################################
do node_x = 0, npx - 1
   imshxx = nd1v(1)
   call divnod( imshxx, npx, node_x, imshx1, imshx2, imshx )
   meshx1(node_x) = imshx1
   meshx(node_x)  = imshx
end do
do node_y = 0, npy - 1
   imshyy = nd1v(2)
   call divnod( imshyy, npy, node_y, imshy1, imshy2, imshy )
   meshy1(node_y) = imshy1
   meshy(node_y)  = imshy
end do
do node_z = 0, npz - 1
   imshzz = nd1v(3)
   call divnod( imshzz, npz, node_z, imshz1, imshz2, imshz )
   meshz1(node_z) = imshz1
   meshz(node_z)  = imshz
end do

nmm  = 0
node_id = 0
outdo: do
   node_x=node_id/(npy*npz)
   node_y=mod(node_id/npz,npy)
   node_z=mod(node_id,npz)
   imshxx = nd1v(1)
   call divnod( imshxx, npx, node_x, imshx1, imshx2, imshx )
   imshyy = nd1v(2)
   call divnod( imshyy, npy, node_y, imshy1, imshy2, imshy )
   imshzz = nd1v(3)
   call divnod( imshzz, npz, node_z, imshz1, imshz2, imshz )

   mul = 1
   nmmrc1(node_id) = nmm + 1
   if( lclust .and. lsphere ) then
       !--- set mesh for cluster calculation ----------------------------------
       icmul1 = 2**(mul+(muls-1)-1)
       icmul2 = 0
       do j = 0, mul+(muls-1) - 2
          icmul2 = icmul2 - 2**j
       end do
       do iz = 1, imshz
          izz = imshz1 + iz - 1
          izz = icmul1*izz + icmul2
          z  = ( dble(izz-1)*rdel(3) - cofmas(3) )*rmaxr3
          z2 = z*z
       do iy = 1, imshy
          iyy = imshy1 + iy - 1
          iyy = icmul1*iyy + icmul2
          y  = ( dble(iyy-1)*rdel(2) - cofmas(2) )*rmaxr2
          y2 = y*y
       do ix = 1, imshx
          ixx = imshx1 + ix - 1
          ixx = icmul1*ixx + icmul2
          x2 = ( dble(ixx-1)*rdel(1) - cofmas(1) )*rmaxr1
          x2 = x2*x2
          r2 = x2 + y2 + z2
          if( r2.lt.rmax2 ) then
              nmm = nmm + 1
!cc             nxyz(ixx,iyy,izz) = nmm
!cc             nx(nmm) = ixx
!cc             ny(nmm) = iyy
!cc             nz(nmm) = izz
            else
!cc             nxyz(ixx,iyy,izz) = 0
          end if
       end do
       end do
       end do
   else
       !--- set mesh for bulk calculation -------------------------------------
       do iz = 1, imshz
          izz = imshz1 + iz - 1
       do iy = 1, imshy
          iyy = imshy1 + iy - 1
       do ix = 1, imshx
          ixx = imshx1 + ix - 1
              nmm = nmm + 1
!cc             nxyz(ixx,iyy,izz) = nmm
!cc             nx(nmm) = ixx
!cc             ny(nmm) = iyy
!cc             nz(nmm) = izz
       end do
       end do
       end do
       !--- end of setting mesh -----------------------------------------------
   end if
   nmmrc2(node_id) = nmm

   node_id = node_id + 1
   if( node_id > nodess - 1 ) exit
end do outdo


!*** error trap
msherr = nmmrc2(mytid)-nmmrc1(mytid)+1 - mshnod(1)
if( nmmrc2(mytid)-nmmrc1(mytid)+1.ne.mshnod(1) ) then
!    if(loutfile(1)) 
                    write(nfile(1),*) 'error in setmsh : iam =',  &
&               myid, nmmrc1(mytid), nmmrc2(mytid), mshnod(1)
    if(loutfile(2)) write(nfile(2),*) 'error in setmsh : iam =',  &
&               myid, nmmrc1(mytid), nmmrc2(mytid), mshnod(1)
end if
call gimax(msherr)
if( msherr.gt.0 ) then
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in setmsh' )
end if


!--- total No. of mesh points     : nmm
do i = 1, nodess
   node_id = i - 1
   nmmcnt(i) = nmmrc2(node_id) - nmmrc1(node_id) + 1
end do
nmmdsp(1) = 0
do i = 2, nodess
   nmmdsp(i) = nmmdsp(i-1) + nmmcnt(i-1)
end do


nmm1    = nmmrc1(mytid)
nmm2    = nmmrc2(mytid)
!#######################################################################


return
end




subroutine cfftmsh( nfile, myid, nodes,  &
& rdel, rdelg, rdelv, ntotfd, nd1vks, mshglb,  &
& lclust, lsphere, rmax, cofmas, rba, rvol,  &
& kfft1, kfft2, kfft3, kfft0 )
!-----------------------------------------------------------------------
!     set FFT mesh for coarse grid
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8,  dimension(3) :: rdel
real*8,  dimension(3,3) :: rdelg
real*8  :: rdelv
integer :: ntotfd
integer, dimension(3) :: nd1vks
integer :: kfft1, kfft2, kfft3, kfft0
integer, dimension(kfft0) :: mshglb
logical :: lclust, lsphere
real*8,  dimension(3) :: rmax
real*8,  dimension(3) :: cofmas
real*8,  dimension(3,3) :: rba
real*8  :: rvol


nd1vks(1) = kfft1
nd1vks(2) = kfft2
nd1vks(3) = kfft3
do i = 1, 3
   rdel(i) = rba(i,i)/dble(nd1vks(i))
   do j = 1, 3
      rdelg(j,i) = rba(j,i)/dble(nd1vks(i))
   end do
end do

rdelv = rvol/dble(kfft0)


do i = 1, kfft0
   mshglb(i) = 0
end do

if( lclust .and. lsphere ) then
    rmaxr1 = 1.d0/rmax(1)
    rmaxr2 = 1.d0/rmax(2)
    rmaxr3 = 1.d0/rmax(3)
!--- set mesh for cluster calculation ----------------------------------
    rmax2  = 1.00001d0
    ntotfd = 0
    do iz = 1, kfft3
       z  = ( dble(iz-1)*rdel(3) - cofmas(3) )*rmaxr3
       z2 = z*z
    do iy = 1, kfft2
       y  = ( dble(iy-1)*rdel(2) - cofmas(2) )*rmaxr2
       y2 = y*y
    do ix = 1, kfft1
       x2 = ( dble(ix-1)*rdel(1) - cofmas(1) )*rmaxr1
       x2 = x2*x2
       r2 = x2 + y2 + z2
       if( r2.lt.rmax2 ) then
           ntotfd = ntotfd + 1
           m2 = ix + kfft1*( (iy-1) + kfft2*(iz-1) )
           mshglb(m2) = ntotfd
       end if
    end do
    end do
    end do
else
!--- set mesh for bulk calculation -------------------------------------
    ntotfd = kfft0
    do i = 1, kfft0
       mshglb(i) = i
    end do
!--- end of setting mesh -----------------------------------------------
end if


return
end




subroutine bdset1( lvacuum,  &
& mshnod, mshnew, mshxyz,  &
& mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz, nd2v,  &
& ismx, ismy, ismz, ismshx, ismshy, ismshz, mulyzs, mulzxs, mulxys,  &
& irmx, irmy, irmz, irmshx, irmshy, irmshz, mulyzr, mulzxr, mulxyr,  &
& ismfpx, ismfpy, ismfpz, irmfpx, irmfpy, irmfpz,  &
& ncomct, jsmshx, jsmshy, jsmshz, nodyzh, nodzxh, nodxyh,  &
& idbuf, idbufr, mshndv, myid, npx, npy, npz, nproc, myx, myy, myz,  &
& nn, myparity, nhit1, nfile, mulgmx, mulgmy, mulgmz )
!-----------------------------------------------------------------------
!  Set No. of data for exchanges boundary among neighbor nodes.
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nn(6), myparity(3)
dimension nfile(*)
logical   lvacuum(3)

dimension      mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension      ismx(0:npx,2), ismshx(mulyzs,2),  &
&              ismy(0:npy,2), ismshy(mulzxs,2),  &
&              ismz(0:npz,2), ismshz(mulxys,2)
dimension      jsmshx(3,nodyzh,2), jsmshy(3,nodzxh,2),  &
&              jsmshz(3,nodxyh,2)
dimension      irmx(0:npx,2), irmshx(mulyzr,2),  &
&              irmy(0:npy,2), irmshy(mulzxr,2),  &
&              irmz(0:npz,2), irmshz(mulxyr,2)
dimension      ismfpx(0:npx,2), ismfpy(0:npy,2), ismfpz(0:npz,2)
dimension      irmfpx(0:npx,2), irmfpy(0:npy,2), irmfpz(0:npz,2)
dimension      idbuf(mshndv), idbufr(mshndv)
dimension ibuf(2), ibufr(2)
dimension ncomct(6,2)
dimension mulgmx(0:npx-1), mulgmy(0:npy-1), mulgmz(0:npz-1)
integer   root


!--- get No. of communication in each direction ------------------------
node_y = 0
node_z = 0
mshxmul = 0
do node_x = 0, npx - 1
   root = node_x*(npy*npz) + node_y*npz + node_z
   mulgmx(node_x) = mshx
   if( npx.gt.1 ) call ibcast(mulgmx(node_x),1,root)
   mshxmul = mshxmul + mulgmx(node_x)
end do
node_z = 0
node_x = 0
mshymul = 0
do node_y = 0, npy - 1
   root = node_x*(npy*npz) + node_y*npz + node_z
   mulgmy(node_y) = mshy
   if( npy.gt.1 ) call ibcast(mulgmy(node_y),1,root)
   mshymul = mshymul + mulgmy(node_y)
end do
node_x = 0
node_y = 0
mshzmul = 0
do node_z = 0, npz - 1
   root = node_x*(npy*npz) + node_y*npz + node_z
   mulgmz(node_z) = mshz
   if( npz.gt.1 ) call ibcast(mulgmz(node_z),1,root)
   mshzmul = mshzmul + mulgmz(node_z)
end do
!--- in x direction
do i = 1, 2
   if( i.eq.1 ) then
       ifugo = 1
     else
       ifugo = -1
   end if
   ncomct(i,1) = 0
   ncomct(i,2) = 0
   ngetdat   = 0
   do
      ncomct(i,1) = ncomct(i,1) + 1
      node_x = mod(myx+ifugo*ncomct(i,1)+npx,npx)
      ngetdat = ngetdat + mulgmx(node_x)
      if( ncomct(i,2).eq.0 .and. ngetdat.ge.nhit1 )  &
&         ncomct(i,2) = ncomct(i,1)
      if( ngetdat >= nd2v ) exit
   end do
end do
if( lvacuum(1) ) then
    if( myx.le.0 ) then
        ncomct(2,1) = 0
        ncomct(2,2) = 0
    end if
    if( myx.ge.npx-1 ) then
        ncomct(1,1) = 0
        ncomct(1,2) = 0
    end if
end if
!--- in y direction
do i = 3, 4
   if( i.eq.3 ) then
       ifugo = 1
     else
       ifugo = -1
   end if
   ncomct(i,1) = 0
   ncomct(i,2) = 0
   ngetdat   = 0
   do
      ncomct(i,1) = ncomct(i,1) + 1
      node_y = mod(myy+ifugo*ncomct(i,1)+npy,npy)
      ngetdat = ngetdat + mulgmy(node_y)
      if( ncomct(i,2).eq.0 .and. ngetdat.ge.nhit1 )  &
&         ncomct(i,2) = ncomct(i,1)
      if( ngetdat >= nd2v ) exit
   end do
end do
if( lvacuum(2) ) then
    if( myy.le.0 ) then
        ncomct(4,1) = 0
        ncomct(4,2) = 0
    end if
    if( myy.ge.npy-1 ) then
        ncomct(3,1) = 0
        ncomct(3,2) = 0
    end if
end if
!--- in z direction
do i = 5, 6
   if( i.eq.5 ) then
       ifugo = 1
     else
       ifugo = -1
   end if
   ncomct(i,1) = 0
   ncomct(i,2) = 0
   ngetdat   = 0
   do
      ncomct(i,1) = ncomct(i,1) + 1
      node_z = mod(myz+ifugo*ncomct(i,1)+npz,npz)
      ngetdat = ngetdat + mulgmz(node_z)
      if( ncomct(i,2).eq.0 .and. ngetdat.ge.nhit1 )  &
&         ncomct(i,2) = ncomct(i,1)
      if( ngetdat >= nd2v ) exit
   end do
end do
if( lvacuum(3) ) then
    if( myz.le.0 ) then
        ncomct(6,1) = 0
        ncomct(6,2) = 0
    end if
    if( myz.ge.npz-1 ) then
        ncomct(5,1) = 0
        ncomct(5,2) = 0
    end if
end if

do i = 1, 2
do j = 1, 6
   call gimax(ncomct(j,i))
end do
end do
!do i = 1, 2
!if( loutfile(i) ) then
!    write(nfile(i),*) '         No. of grid :',  &
!&                      mshxmul, mshymul, mshzmul
!    write(nfile(i),*) 'No. of grid/node in x:',  &
!&                   ( mulgmx(ii), ii = 0, npx - 1 )
!    write(nfile(i),*) '                 in y:',  &
!&                   ( mulgmy(ii), ii = 0, npy - 1 )
!    write(nfile(i),*) '                 in z:',  &
!&                   ( mulgmz(ii), ii = 0, npz - 1 )
!    write(nfile(i),*) '  data exchange in x :',  &
!&                  ncomct(1,1),ncomct(2,1),ncomct(1,2),ncomct(2,2)
!    write(nfile(i),*) '                in y :',  &
!&                  ncomct(3,1),ncomct(4,1),ncomct(3,2),ncomct(4,2)
!    write(nfile(i),*) '                in z :',  &
!&                  ncomct(5,1),ncomct(6,1),ncomct(5,2),ncomct(6,2)
!end if
!end do


mshnew = mshnod


!-------Loop over the lower & higher x directions ----------------------
kd = 1
do 200 kdd= -1,0

!---------Neighbor node ID 
    ku=2*kd+kdd
    inode=nn(ku)

ismfpx(0,2+kdd) = 0
  ismx(0,2+kdd) = 0
irmfpx(0,2+kdd) = 0
  irmx(0,2+kdd) = 0
ixmin = mshx + 1
ixmax = 0
do 200 ncom = 1, ncomct(ku,1)
   i = 2+kdd
   ismfpx(ncom,i) = ismfpx(ncom-1,i)
     ismx(ncom,i) =   ismx(ncom-1,i)
if( i.eq.1 ) then
!--- to send data to neighbor node in -x direction ---------------------
    node_x = mod(myx+(ncom-1)+npx,npx)
    ixmin = ixmax + 1
    ixmax = ixmax + mulgmx(node_x)
    if( .not.lvacuum(1) .or. myx.gt.0 ) then
!cc              do ix = 1, nhit1
        do ix = ixmin, min( nhit1, ixmax )
        do iy = 1, mshy
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismx(ncom,i) = ismx(ncom,i) + 1
               jsmshx(1,ismx(ncom,i),i) = ix
               jsmshx(2,ismx(ncom,i),i) = iy
               jsmshx(3,ismx(ncom,i),i) = iz
                 ismshx(ismx(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( ixmin .le. min( nhit1, ixmax ) )  &
&                                    ismfpx(ncom,i) = ismx(ncom,i)

!cc              do ix = nhit1+1, nd2v
        do ix = max( nhit1+1, ixmin ), min( nd2v, ixmax )
        do iy = 1, mshy
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismx(ncom,i) = ismx(ncom,i) + 1
               jsmshx(1,ismx(ncom,i),i) = ix
               jsmshx(2,ismx(ncom,i),i) = iy
               jsmshx(3,ismx(ncom,i),i) = iz
                 ismshx(ismx(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
  else
!--- to send data to neighbor node in +x direction ---------------------
    node_x = mod(myx-(ncom-1)+npx,npx)
    ixmax = ixmin - 1
    ixmin = ixmax - mulgmx(node_x) + 1
    if( .not.lvacuum(1) .or. myx.lt.npx-1 ) then
!cc              do ix = mshx-nhit1+1, mshx
        do ix = max( mshx-nhit1+1, ixmin ), min( mshx, ixmax )
        do iy = 1, mshy
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismx(ncom,i) = ismx(ncom,i) + 1
               jsmshx(1,ismx(ncom,i),i) = ix - mshx
               jsmshx(2,ismx(ncom,i),i) = iy
               jsmshx(3,ismx(ncom,i),i) = iz
                 ismshx(ismx(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( max( mshx-nhit1+1, ixmin ) .le. min( mshx, ixmax ) )  &
&                                    ismfpx(ncom,i) = ismx(ncom,i)

        do ix =max( mshx-nd2v+1, ixmin ), min( mshx-nhit1, ixmax )
        do iy = 1, mshy
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismx(ncom,i) = ismx(ncom,i) + 1
               jsmshx(1,ismx(ncom,i),i) = ix - mshx
               jsmshx(2,ismx(ncom,i),i) = iy
               jsmshx(3,ismx(ncom,i),i) = iz
                 ismshx(ismx(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
end if


!---------Send the # of boundary data
    nsd = ismx(ncom,2+kdd) - ismx(ncom-1,2+kdd)
    nsd3= nsd*3

    ibuf(1) =   ismx(ncom,2+kdd) -   ismx(ncom-1,2+kdd)
    ibuf(2) = ismfpx(ncom,2+kdd) - ismfpx(ncom-1,2+kdd)
!---------Even: send & recv
    if (myparity(kd).eq.0) then
      call cisend(100+ku,ibuf,2,inode,0)
      call cirecv(110+ku,ibufr,2,0)
!---------Odd: recv & send
    else if (myparity(kd).eq.1) then
      call cirecv(100+ku,ibufr,2,0)
      call cisend(110+ku,ibuf,2,inode,0)
!---------Exchange information with myself
    else
      ibufr(1) = ibuf(1)
      ibufr(2) = ibuf(2)
    end if

    nrc = ibufr(1)
    nrc3= nrc*3

!---------Message buffering
    ii = 0
    do i = ismx(ncom-1,2+kdd) + 1, ismx(ncom,2+kdd)
    do j = 1, 3
       ii = ii + 1
       idbuf(ii)=jsmshx(j,i,2+kdd)
    end do
    end do

!---------Even: send & recv, if not empty
    if (myparity(kd).eq.0) then
      if (nsd.ne.0) call cisend(120+ku,idbuf,nsd3,inode,0)
      if (nrc.ne.0) call cirecv(130+ku,idbufr,nrc3,0)
!---------Odd: recv & send, if not empty
    else if (myparity(kd).eq.1) then
      if (nrc.ne.0) call cirecv(120+ku,idbufr,nrc3,0)
      if (nsd.ne.0) call cisend(130+ku,idbuf,nsd3,inode,0)
!---------Exchange information with myself
    else
      do i=1,nrc3
        idbufr(i)=idbuf(i)
      end do
    end if

!---------Message storing
      irmx(ncom,2+kdd) =   irmx(ncom-1,2+kdd) + nrc
    irmfpx(ncom,2+kdd) = irmfpx(ncom-1,2+kdd) + ibufr(2)
    if( 2+kdd.eq.1 ) then
        mtasu = mshx
      else
        mtasu = 0
    end if
    ii = 0
    do i = irmx(ncom-1,2+kdd) + 1, irmx(ncom,2+kdd)
       ii = ii + 3
       ix = idbufr(ii-2) + mtasu
       iy = idbufr(ii-1)
       iz = idbufr(ii-0)
       mshnew = mshnew + 1
       mshxyz(ix,iy,iz) = mshnew
       irmshx(i,2+kdd)  = mshnew
    end do

!---------Message passing in one direction done-------------------------

!---------Internode synchronization
    if( nproc.gt.1 ) call gsync

!-------Enddo over the lower & higher x directions ---------------------
200 continue


!-------Loop over the lower & higher y directions ----------------------
kd = 2
do 210 kdd= -1,0

!---------Neighbor node ID 
    ku=2*kd+kdd
    inode=nn(ku)

ismfpy(0,2+kdd) = 0
  ismy(0,2+kdd) = 0
irmfpy(0,2+kdd) = 0
  irmy(0,2+kdd) = 0
iymin = mshy + 1
iymax = 0
do 210 ncom = 1, ncomct(ku,1)
   i = 2+kdd
   ismfpy(ncom,i) = ismfpy(ncom-1,i)
     ismy(ncom,i) =   ismy(ncom-1,i)
if( i.eq.1 ) then
!--- to send data to neighbor node in -y direction ---------------------
    node_y = mod(myy+(ncom-1)+npy,npy)
    iymin = iymax + 1
    iymax = iymax + mulgmy(node_y)
    if( .not.lvacuum(2) .or. myy.gt.0 ) then
!cc              do ix = 1, mshx
        do ix = 0, mshx + 1
!cc              do iy = 1, nhit1
        do iy = iymin, min( nhit1, iymax )
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismy(ncom,i) = ismy(ncom,i) + 1
               jsmshy(1,ismy(ncom,i),i) = ix
               jsmshy(2,ismy(ncom,i),i) = iy
               jsmshy(3,ismy(ncom,i),i) = iz
                 ismshy(ismy(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( iymin .le. min( nhit1, iymax ) )  &
&                                    ismfpy(ncom,i) = ismy(ncom,i)

        do ix = 1, mshx
!cc              do iy = nhit1+1, nd2v
        do iy = max( nhit1+1, iymin ), min( nd2v, iymax )
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismy(ncom,i) = ismy(ncom,i) + 1
               jsmshy(1,ismy(ncom,i),i) = ix
               jsmshy(2,ismy(ncom,i),i) = iy
               jsmshy(3,ismy(ncom,i),i) = iz
                 ismshy(ismy(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
  else
!--- to send data to neighbor node in +y direction ---------------------
    node_y = mod(myy-(ncom-1)+npy,npy)
    iymax = iymin - 1
    iymin = iymax - mulgmy(node_y) + 1
    if( .not.lvacuum(2) .or. myy.lt.npy-1 ) then
!cc              do ix = 1, mshx
        do ix = 0, mshx + 1
!cc              do iy = mshy-nhit1+1, mshy
        do iy = max( mshy-nhit1+1, iymin ), min( mshy, iymax )
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismy(ncom,i) = ismy(ncom,i) + 1
               jsmshy(1,ismy(ncom,i),i) = ix
               jsmshy(2,ismy(ncom,i),i) = iy - mshy
               jsmshy(3,ismy(ncom,i),i) = iz
                 ismshy(ismy(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( max( mshy-nhit1+1, iymin ) .le. min( mshy, iymax ) )  &
&                                    ismfpy(ncom,i) = ismy(ncom,i)

        do ix = 1, mshx
!cc              do iy = mshy-nd2v+1, mshy-nhit1
        do iy =max( mshy-nd2v+1, iymin ), min( mshy-nhit1, iymax )
        do iz = 1, mshz
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismy(ncom,i) = ismy(ncom,i) + 1
               jsmshy(1,ismy(ncom,i),i) = ix
               jsmshy(2,ismy(ncom,i),i) = iy - mshy
               jsmshy(3,ismy(ncom,i),i) = iz
                 ismshy(ismy(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
end if


!---------Send the # of boundary data
    nsd = ismy(ncom,2+kdd) - ismy(ncom-1,2+kdd)
    nsd3= nsd*3

    ibuf(1) =   ismy(ncom,2+kdd) -   ismy(ncom-1,2+kdd)
    ibuf(2) = ismfpy(ncom,2+kdd) - ismfpy(ncom-1,2+kdd)
!---------Even: send & recv
    if (myparity(kd).eq.0) then
      call cisend(100+ku,ibuf,2,inode,0)
      call cirecv(110+ku,ibufr,2,0)
!---------Odd: recv & send
    else if (myparity(kd).eq.1) then
      call cirecv(100+ku,ibufr,2,0)
      call cisend(110+ku,ibuf,2,inode,0)
!---------Exchange information with myself
    else
      ibufr(1) = ibuf(1)
      ibufr(2) = ibuf(2)
    end if

    nrc = ibufr(1)
    nrc3= nrc*3

!---------Message buffering
    ii = 0
    do i = ismy(ncom-1,2+kdd) + 1, ismy(ncom,2+kdd)
    do j = 1, 3
       ii = ii + 1
       idbuf(ii)=jsmshy(j,i,2+kdd)
    end do
    end do

!---------Even: send & recv, if not empty
    if (myparity(kd).eq.0) then
      if (nsd.ne.0) call cisend(120+ku,idbuf,nsd3,inode,0)
      if (nrc.ne.0) call cirecv(130+ku,idbufr,nrc3,0)
!---------Odd: recv & send, if not empty
    else if (myparity(kd).eq.1) then
      if (nrc.ne.0) call cirecv(120+ku,idbufr,nrc3,0)
      if (nsd.ne.0) call cisend(130+ku,idbuf,nsd3,inode,0)
!---------Exchange information with myself
    else
      do i=1,nrc3
        idbufr(i)=idbuf(i)
      end do
    end if

!---------Message storing
      irmy(ncom,2+kdd) =   irmy(ncom-1,2+kdd) + nrc
    irmfpy(ncom,2+kdd) = irmfpy(ncom-1,2+kdd) + ibufr(2)
    if( 2+kdd.eq.1 ) then
        mtasu = mshy
      else
        mtasu = 0
    end if
    ii = 0
    do i = irmy(ncom-1,2+kdd) + 1, irmy(ncom,2+kdd)
       ii = ii + 3
       ix = idbufr(ii-2)
       iy = idbufr(ii-1) + mtasu
       iz = idbufr(ii-0)
       mshnew = mshnew + 1
       mshxyz(ix,iy,iz) = mshnew
       irmshy(i,2+kdd)  = mshnew
    end do

!---------Message passing in one direction done-------------------------

!---------Internode synchronization
    if( nproc.gt.1 ) call gsync

!-------Enddo over the lower & higher y directions ---------------------
210 continue


!-------Loop over the lower & higher z directions ----------------------
kd = 3
do 220 kdd= -1,0

!---------Neighbor node ID 
    ku=2*kd+kdd
    inode=nn(ku)

ismfpz(0,2+kdd) = 0
  ismz(0,2+kdd) = 0
irmfpz(0,2+kdd) = 0
  irmz(0,2+kdd) = 0
izmin = mshz + 1
izmax = 0
do 220 ncom = 1, ncomct(ku,1)
   i = 2+kdd
   ismfpz(ncom,i) = ismfpz(ncom-1,i)
     ismz(ncom,i) =   ismz(ncom-1,i)
if( i.eq.1 ) then
!--- to send data to neighbor node in -z direction ---------------------
    node_z = mod(myz+(ncom-1)+npz,npz)
    izmin = izmax + 1
    izmax = izmax + mulgmz(node_z)
    if( .not.lvacuum(3) .or. myz.gt.0 ) then
!cc              do ix = 1, mshx
!cc              do iy = 1, mshy
        do ix = 0, mshx + 1
        do iy = 0, mshy + 1
!cc              do iz = 1, nhit1
        do iz = izmin, min( nhit1, izmax )
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismz(ncom,i) = ismz(ncom,i) + 1
               jsmshz(1,ismz(ncom,i),i) = ix
               jsmshz(2,ismz(ncom,i),i) = iy
               jsmshz(3,ismz(ncom,i),i) = iz
                 ismshz(ismz(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( izmin .le. min( nhit1, izmax ) )  &
&                                    ismfpz(ncom,i) = ismz(ncom,i)

        do ix = 1, mshx
        do iy = 1, mshy
!cc              do iz = nhit1+1, nd2v
        do iz = max( nhit1+1, izmin ), min( nd2v, izmax )
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismz(ncom,i) = ismz(ncom,i) + 1
               jsmshz(1,ismz(ncom,i),i) = ix
               jsmshz(2,ismz(ncom,i),i) = iy
               jsmshz(3,ismz(ncom,i),i) = iz
                 ismshz(ismz(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
  else
!--- to send data to neighbor node in +z direction ---------------------
    node_z = mod(myz-(ncom-1)+npz,npz)
    izmax = izmin - 1
    izmin = izmax - mulgmz(node_z) + 1
    if( .not.lvacuum(3) .or. myz.lt.npz-1 ) then
!cc              do ix = 1, mshx
!cc              do iy = 1, mshy
        do ix = 0, mshx + 1
        do iy = 0, mshy + 1
!cc              do iz = mshz-nhit1+1, mshz
        do iz = max( mshz-nhit1+1, izmin ), min( mshz, izmax )
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismz(ncom,i) = ismz(ncom,i) + 1
               jsmshz(1,ismz(ncom,i),i) = ix
               jsmshz(2,ismz(ncom,i),i) = iy
               jsmshz(3,ismz(ncom,i),i) = iz - mshz
                 ismshz(ismz(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
        if( max( mshz-nhit1+1, izmin ) .le. min( mshz, izmax ) )  &
&                                    ismfpz(ncom,i) = ismz(ncom,i)

        do ix = 1, mshx
        do iy = 1, mshy
!cc              do iz = mshz-nd2v+1, mshz-nhit1
        do iz =max( mshz-nd2v+1, izmin ), min( mshz-nhit1, izmax )
           if( mshxyz(ix,iy,iz).gt.0 ) then
               ismz(ncom,i) = ismz(ncom,i) + 1
               jsmshz(1,ismz(ncom,i),i) = ix
               jsmshz(2,ismz(ncom,i),i) = iy
               jsmshz(3,ismz(ncom,i),i) = iz - mshz
                 ismshz(ismz(ncom,i),i) = mshxyz(ix,iy,iz)
           end if
        end do
        end do
        end do
    end if
end if


!---------Send the # of boundary data
    nsd = ismz(ncom,2+kdd) - ismz(ncom-1,2+kdd)
    nsd3= nsd*3

    ibuf(1) =   ismz(ncom,2+kdd) -   ismz(ncom-1,2+kdd)
    ibuf(2) = ismfpz(ncom,2+kdd) - ismfpz(ncom-1,2+kdd)
!---------Even: send & recv
    if (myparity(kd).eq.0) then
      call cisend(100+ku,ibuf,2,inode,0)
      call cirecv(110+ku,ibufr,2,0)
!---------Odd: recv & send
    else if (myparity(kd).eq.1) then
      call cirecv(100+ku,ibufr,2,0)
      call cisend(110+ku,ibuf,2,inode,0)
!---------Exchange information with myself
    else
      ibufr(1) = ibuf(1)
      ibufr(2) = ibuf(2)
    end if

    nrc = ibufr(1)
    nrc3= nrc*3

!---------Message buffering
    ii = 0
    do i = ismz(ncom-1,2+kdd) + 1, ismz(ncom,2+kdd)
    do j = 1, 3
       ii = ii + 1
       idbuf(ii)=jsmshz(j,i,2+kdd)
    end do
    end do

!---------Even: send & recv, if not empty
    if (myparity(kd).eq.0) then
      if (nsd.ne.0) call cisend(120+ku,idbuf,nsd3,inode,0)
      if (nrc.ne.0) call cirecv(130+ku,idbufr,nrc3,0)
!---------Odd: recv & send, if not empty
    else if (myparity(kd).eq.1) then
      if (nrc.ne.0) call cirecv(120+ku,idbufr,nrc3,0)
      if (nsd.ne.0) call cisend(130+ku,idbuf,nsd3,inode,0)
!---------Exchange information with myself
    else
      do i=1,nrc3
        idbufr(i)=idbuf(i)
      end do
    end if

!---------Message storing
      irmz(ncom,2+kdd) =   irmz(ncom-1,2+kdd) + nrc
    irmfpz(ncom,2+kdd) = irmfpz(ncom-1,2+kdd) + ibufr(2)
    if( 2+kdd.eq.1 ) then
        mtasu = mshz
      else
        mtasu = 0
    end if
    ii = 0
    do i = irmz(ncom-1,2+kdd) + 1, irmz(ncom,2+kdd)
       ii = ii + 3
       ix = idbufr(ii-2)
       iy = idbufr(ii-1)
       iz = idbufr(ii-0) + mtasu
       mshnew = mshnew + 1
       mshxyz(ix,iy,iz) = mshnew
       irmshz(i,2+kdd)  = mshnew
    end do

!---------Message passing in one direction done-------------------------

!---------Internode synchronization
    if( nproc.gt.1 ) call gsync

!-------Enddo over the lower & higher z directions ---------------------
220 continue


!--- error trap
ierror = 0
kd = 1
kdddo1: do kdd = -1, 0
   ku=2*kd+kdd
   do ncom = ncomct(ku,1), 1, -1
      i = 2+kdd
      if( ismfpx(ncom,i)-ismfpx(ncom-1,i).gt.0 ) then
          if( ncom.gt.ncomct(ku,2) ) then
              ierror = 1
!              if(loutfile(1)) 
                              write(nfile(1),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismx(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpx(j,i), j = 0, ncomct(ku,1) )
              if(loutfile(2)) write(nfile(2),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismx(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpx(j,i), j = 0, ncomct(ku,1) )
              exit kdddo1
          end if
          exit
      end if
   end do
end do kdddo1

if( ierror == 0 ) then
kd = 2
kdddo2: do kdd = -1, 0
   ku=2*kd+kdd
   do ncom = ncomct(ku,1), 1, -1
      i = 2+kdd
      if( ismfpy(ncom,i)-ismfpy(ncom-1,i).gt.0 ) then
          if( ncom.gt.ncomct(ku,2) ) then
              ierror = 1
!              if(loutfile(1)) 
                              write(nfile(1),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismy(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpy(j,i), j = 0, ncomct(ku,1) )
              if(loutfile(2)) write(nfile(2),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismy(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpy(j,i), j = 0, ncomct(ku,1) )
              exit kdddo2
          end if
          exit
      end if
   end do
end do kdddo2
end if

if( ierror == 0 ) then
kd = 3
kdddo3: do kdd = -1, 0
   ku=2*kd+kdd
   do ncom = ncomct(ku,1), 1, -1
      i = 2+kdd
      if( ismfpz(ncom,i)-ismfpz(ncom-1,i).gt.0 ) then
          if( ncom.gt.ncomct(ku,2) ) then
              ierror = 1
!              if(loutfile(1)) 
                              write(nfile(1),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismz(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpz(j,i), j = 0, ncomct(ku,1) )
              if(loutfile(2)) write(nfile(2),*) 'error in data exchange : iam =',  &
&                           myid, ku, ncomct(ku,1), ncomct(ku,2),  &
&                           ( ismz(j,i),   j = 0, ncomct(ku,1) ),  &
&                           ( ismfpz(j,i), j = 0, ncomct(ku,1) )
              exit kdddo3
          end if
          exit
      end if
   end do
end do kdddo3
end if

call gimax(ierror)
if( ierror.gt.0 ) then
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in bdset1' )
end if


return
end




subroutine ntset( myx, myy, myz, nn, myparity, npx, npy, npz )
!----------------------------------------------------------------------c
!  Prepares a neighbor-node-ID table, NN, & a shift-vector table, SV,
!  for internode message passing.  Also prepares the node parity table,
!  MYPARITY.
!----------------------------------------------------------------------c
implicit real*8 ( a-h, o-z )
integer :: npx, npy, npz
integer :: myx, myy, myz
integer, dimension(6) :: nn
integer, dimension(3) :: myparity

!-----Integer vectors to specify the six neighbor nodes
integer, dimension(3,6) :: iv =  &
& reshape((/-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1/), (/3,6/))

do ku = 1, 6
!-------Vector index of neighbor ku
  k1x=mod(myx+iv(1,ku)+npx,npx)
  k1y=mod(myy+iv(2,ku)+npy,npy)
  k1z=mod(myz+iv(3,ku)+npz,npz)
!-------Scalar neighbor ID, NN
  nn(ku)=k1x*(npy*npz)+k1y*npz+k1z
end do

!-----Set up the node parity table, MY_PARITY
call set_parity(myx,npx,iparity)
myparity(1)=iparity
call set_parity(myy,npy,iparity)
myparity(2)=iparity
call set_parity(myz,npz,iparity)
myparity(3)=iparity

return
end


subroutine set_parity(myx,nx,iparity)
!----------------------------------------------------------------------c
implicit real*8(a-h,o-z)
!----------------------------------------------------------------------c
if(nx.eq.1)then
  iparity=2
elseif(mod(myx,2).eq.0)then
  iparity=0
else
  iparity=1
end if

return
end




subroutine baset( nfile, myid, nodes,  &
& lvacuum, mshx1, mshy1, mshz1, mshx, mshy, mshz, lclust,  &
& ntype, nhk1, nhk2, iatmpt, nionnew, rdelg, hcell, ratm, cratm,  &
& lkbppi, lking, rmxnlx, xc, ntall )
!-----------------------------------------------------------------------
!  set No. of boundary atoms.
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension nhk1(*), nhk2(*), iatmpt(*)
dimension ratm(3,*)
dimension cratm(3,*)
dimension rdelg(3,3)
dimension hcell(3,3)
logical   lvacuum(3)
logical   lclust
logical, dimension(ntype) :: lkbppi, lking
real*8,  dimension(3) :: xc

dimension celbnd(2,6)
dimension rmxdrn(3)
real*8,  dimension(3) :: d1vr


nionnew = 0
do it = 1, ntype
   if( lkbppi(it) .and. lking(it) ) then
       do i = nhk1(it), nhk2(it)
          nionnew = nionnew + 1
          cratm(1,nionnew) = ratm(1,i)
          cratm(2,nionnew) = ratm(2,i)
          cratm(3,nionnew) = ratm(3,i)
       end do
   end if
end do
if( lclust ) return


do i = 1, 3
!   d1vr(i) = rdel(i)/hcell(i,i)
   d1vr(i) = vecratio(rdelg(1,i),hcell(1,i))
end do
! --- set boundary
do i = 1, 3
   celbnd(1,2*i-1) = 0.d0
   celbnd(2,2*i-1) = - d1vr(i)
   celbnd(1,2*i  ) = 1.d0
   celbnd(2,2*i  ) = 1.d0
end do

do i = 1, 3
   rmxdrn(i) = xc(i)
end do


! --- check atoms on boundary
xsmall = 1.d-14
do ix = 1, 3
do i  = 1, nionnew
   if( cratm(ix,i) .le. 0.d0 ) cratm(ix,i) = xsmall
   if( cratm(ix,i) .ge. 1.d0 ) cratm(ix,i) = 1.d0 - xsmall
end do
end do


ierror = 0
!=====Main loop over x, y, & z directions starts=======================c
kddo: do kd = 1, 3
   if(      kd.eq.1 ) then
            kd1 = 2
            kd2 = 3
   else if( kd.eq.2 ) then
            kd1 = 3
            kd2 = 1
   else if( kd.eq.3 ) then
            kd1 = 1
            kd2 = 2
   end if

   if( .not.lvacuum(kd) ) then

   !-------For both low & high directions
   kdddo: do kdd= -1,0
      nion1 = 1
      do
         nionold = nionnew
         kul = 2*kd+kdd
         if( kdd.eq.-1 ) then
             fugo = 1.d0
           else
             fugo = -1.d0
         end if
         !---   pre-scan the residents
         lsb = 0
         do i_nod = nion1, nionnew
            if( fugo*(cratm(kd,i_nod)-celbnd(1,kul)).gt.0.d0 .and.  &
&                 abs(cratm(kd,i_nod)-celbnd(2,kul)).le.rmxdrn(kd) ) then
                lsb = lsb + 1
            end if
         end do
         if( lsb.le.0 ) cycle kdddo
         !----------------------------------- break if no atom is transferred.
         !--- error trap
         if( nionnew+lsb .gt. ntall ) then
             ierror = 1
             exit kddo
         end if

         !-------Scan all the residents & copies to find boundary atoms
         lsb = nionnew
         do i_nod = nion1, nionnew
            if( fugo*(cratm(kd,i_nod)-celbnd(1,kul)).gt.0.d0 .and.  &
&                 abs(cratm(kd,i_nod)-celbnd(2,kul)).le.rmxdrn(kd) ) then
               lsb = lsb + 1
               iatmpt(lsb)   = iatmpt(i_nod)
               cratm(kd, lsb) = cratm(kd, i_nod) + fugo
               cratm(kd1,lsb) = cratm(kd1,i_nod)
               cratm(kd2,lsb) = cratm(kd2,i_nod)
            end if
         !-------Enddo resident & copied atoms i 
         end do

         !---------Increase the # of received boundar atoms
         nionnew = lsb

         !-------Enddo over the lower & higher directions
         nion1 = nionold + 1
      end do
   end do kdddo

   end if
!-----Enddo x, y, & z directions
end do kddo

!--- error trap
call gimax(ierror)
if( ierror.ne.0 ) then
    if(loutfile(1)) write(nfile(1),*) 'error in baset: nionnew.gt.ntall'
    if(loutfile(2)) write(nfile(2),*) 'error in baset: nionnew.gt.ntall'
!      -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in baset' )
end if


!--- Error trap: check No. of ions in each node
nionmx = nionnew
call gimax(nionmx)
!if(loutfile(1)) write(nfile(1),*) ' the largest No. of atoms      :',nionmx
!if(loutfile(2)) write(nfile(2),*) ' the largest No. of atoms      :',nionmx
!--- Error trap: check No. of ions exchanged
ierror = 0
if( nionnew.ne.nionmx ) then
    ierror = 1
!    if(loutfile(1)) 
                    write(nfile(1),*) 'error in baset(1) iam=',  &
&                      myid, nionnew, nionmx
    if(loutfile(2)) write(nfile(2),*) 'error in baset(1) iam=',  &
&                      myid, nionnew, nionmx
end if
call gimax(ierror)
if( ierror.ne.0 ) then
!       -----Finalize the parallel environment
!    call end_parallel(ierr)
!    stop
    call fstop( nfile, myid, nodes, '*** error in baset' )
end if

!----------------------------------------------


if( .not.lclust ) then
    do i = 1, nionnew
       q1 = cratm(1,i)
       q2 = cratm(2,i)
       q3 = cratm(3,i)
       do j = 1, 3
          cratm(j,i) = hcell(j,1)*q1 + hcell(j,2)*q2  &
&                                    + hcell(j,3)*q3
       end do
    end do
end if


return
end




