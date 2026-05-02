###########################################################
# QXMD build helper functions
###########################################################

if(NOT TARGET QXMDGlobalCompileDefinitions)
    add_library(QXMDGlobalCompileDefinitions INTERFACE)
endif()

function(qxmd_add_global_compile_definitions DEF)
    if(TARGET QXMDGlobalCompileDefinitions)
        set_property(TARGET QXMDGlobalCompileDefinitions APPEND
                     PROPERTY INTERFACE_COMPILE_DEFINITIONS ${DEF} ${ARGN})
    endif()
endfunction(qxmd_add_global_compile_definitions)

function(qxmd_get_global_compile_definitions OUTVAR)
    if(TARGET QXMDGlobalCompileDefinitions)
        get_target_property(defs QXMDGlobalCompileDefinitions
            INTERFACE_COMPILE_DEFINITIONS)
        set(${OUTVAR} ${defs} PARENT_SCOPE)
    endif()
endfunction(qxmd_get_global_compile_definitions)

function(qxmd_get_fortran_cpp_flag OUTVAR)
    if(DEFINED Fortran_PREPROCESSOR_FLAGS)
        set(${OUTVAR} "${Fortran_PREPROCESSOR_FLAGS}" PARENT_SCOPE)
    else()
        # TODO actual flag check
	if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
	  set (${OUTVAR} "-fpp;-allow;nofpp_comments" PARENT_SCOPE)
        else() 
          set(${OUTVAR} "-cpp" PARENT_SCOPE)
        endif() 
    endif()
endfunction(qxmd_get_fortran_cpp_flag)

function(qxmd_preprocess_source IN OUT)
    qxmd_get_global_compile_definitions(global_defs)
    foreach(DEF ${global_defs})
        list(APPEND global_flags "-D${DEF}")
    endforeach()
    get_filename_component(out_dir ${OUT} DIRECTORY)
    if(NOT EXISTS ${out_dir})
        file(MAKE_DIRECTORY ${out_dir})
    endif()
    add_custom_command(
        OUTPUT ${OUT}
        COMMAND cpp -P ${global_flags} -E ${IN} > ${OUT}
        MAIN_DEPENDENCY ${IN}
        COMMENT "Preprocessing ${IN}"
        VERBATIM)    
endfunction(qxmd_preprocess_source)

function(qxmd_fix_fortran_modules TGT)
    set(targets ${TGT} ${ARGN})
    foreach(tgt IN LISTS targets)
        get_target_property(tgt_type ${tgt} TYPE)
        # All of the following target modifications make
        # sense on non-interfaces only
        if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
            get_target_property(tgt_module_dir ${tgt} Fortran_MODULE_DIRECTORY)
            # set module path to tgt_binary_dir/mod
            get_target_property(tgt_binary_dir ${tgt} BINARY_DIR)
            set_target_properties(${tgt}
                PROPERTIES
                    Fortran_MODULE_DIRECTORY ${tgt_binary_dir}/mod/${TGT})
            # make module directory available for clients of TGT 
            target_include_directories(${tgt}
                PUBLIC
                    $<BUILD_INTERFACE:${tgt_binary_dir}/mod/${TGT}>
                INTERFACE
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qxmd/${TGT}>)
        endif()
    endforeach()
endfunction(qxmd_fix_fortran_modules)

function(qxmd_add_executable EXE)
    add_executable(${EXE} ${ARGN})
    _qxmd_add_target(${EXE} ${ARGN})
endfunction(qxmd_add_executable)

function(qxmd_add_library LIB)
    add_library(${LIB} ${ARGN})
    _qxmd_add_target(${LIB} ${ARGN})
endfunction(qxmd_add_library)

# Only use this one for Fortran targets
function(_qxmd_add_target TGT)
    if(TARGET QXMDGlobalCompileDefinitions)
        target_link_libraries(${TGT} PUBLIC QXMDGlobalCompileDefinitions)
    endif()
    qxmd_fix_fortran_modules(${TGT})
    qxmd_get_fortran_cpp_flag(f_cpp_flag)
    target_compile_options(${TGT} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${f_cpp_flag}>)
endfunction(_qxmd_add_target)

function(qxmd_install_targets TGT)
    set(targets ${TGT} ${ARGN})
    install(TARGETS ${targets}
        EXPORT qxmdTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # Windows needs RUNTIME also for libraries
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qxmd/${TGT})
    # Retrieving non-whitelisted properties leads to an hard
    # error, let's skip the following section for interface
    # targets. See here for details:
    # https://gitlab.kitware.com/cmake/cmake/issues/17640
    foreach(tgt IN LISTS targets)
        get_target_property(tgt_type ${tgt} TYPE)
        if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
            # If the target generates Fortran modules, make sure
            # to install them as well to a proper location
            get_target_property(tgt_module_dir ${tgt} Fortran_MODULE_DIRECTORY)
            if(tgt_module_dir)
                install(DIRECTORY ${tgt_module_dir}/
                    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qxmd/${TGT})
            endif()
        endif()        
    endforeach()
endfunction(qxmd_install_targets)

function(qxmd_ensure_build_type DEFAULT)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to '${DEFAULT}' as none was specified")
        set(CMAKE_BUILD_TYPE "${DEFAULT}"
            CACHE STRING "Choose the type of build." FORCE)
        # Set the possible values of build type for cmake-gui
        set_property(CACHE CMAKE_BUILD_TYPE
            PROPERTY
                STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
    endif()
endfunction(qxmd_ensure_build_type)

if(TARGET QXMDGlobalCompileDefinitions)
    qxmd_install_targets(QXMDGlobalCompileDefinitions)
endif()
