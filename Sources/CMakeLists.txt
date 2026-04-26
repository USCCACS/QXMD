set(sources
    allinc.f90
    modules.f90
    ftmain.f90
    constant.f90
    input.f90
    pwlda.f90
    pwreset.f90
    kpoints.f90
    bulk.f90
    stress.f90
    openfile.f90
    reclat.f90
    ba.f90
    eigen.f90
    schmidt.f90
    fermi.f90
    pcc.f90
    vxc.f90
    chgdns.f90
    funcs.f90
    ppkb.f90
    ppkbr.f90
    ppvd.f90
    nlkbpp.f90
    nlkbppr.f90
    force.f90
    tddft_fssh.f90
    sbalin.f90
    savedd.f90
    engrad.f90
    commwv.f90
    fft3.f90
    datadump.f90
    remd_modules.f90
    remd.f90
    remd_in.f90
    remd_disp.f90)

qxmd_add_library(qxmd_lib ${sources})
target_link_libraries(qxmd_lib PUBLIC qxmd_fft)

#######################################################
# QXMD executable
#######################################################

# This will build the default MPI build.
# For pvm_comm and no_comm are also availabe but not commonly used.
qxmd_add_library(mpi_wrapper mpi_comm.f90)
target_link_libraries(mpi_wrapper PUBLIC qxmd_mpi_fortran)

qxmd_add_executable(qxmd_exe ftmain.f90)
target_link_libraries(qxmd_exe PUBLIC qxmd_lib mpi_wrapper)

#######################################################
# Enable cholesky decomposition 
#######################################################
if(QXMD_ENABLE_CHOLESKY)
    qxmd_add_global_compile_definitions(PCHOLESKY=1)
else()
    qxmd_add_global_compile_definitions(PCHOLESKY=0)
endif()

#######################################################
# QXMD real precision option
#######################################################
qxmd_add_global_compile_definitions(CCREAL4=0)

set_target_properties(qxmd_exe PROPERTIES OUTPUT_NAME qxmd)
