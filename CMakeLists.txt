# (C) Copyright 2022 .
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


cmake_minimum_required( VERSION 3.12)
project( imsproc VERSION 2022.10 LANGUAGES Fortran )

# ecbuild integration
find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( imsproc_compiler_flags )

################################################################################
# Dependencies
################################################################################

find_package(OpenMP COMPONENTS C Fortran)
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran )

################################################################################
# Sources
################################################################################

#add_definitions ( -DCPP=0 )

set( IMSPROC_LINKER_LANGUAGE Fortran )

list ( APPEND sorc_files
./sorc/IMSaggregate_mod.f90
)

ecbuild_add_library( TARGET   imsproc
                     SOURCES  ${sorc_files}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${IMSPROC_LINKER_LANGUAGE}
                    )


target_link_libraries(imsproc PUBLIC NetCDF::NetCDF_Fortran)
target_link_libraries(imsproc PUBLIC MPI::MPI_Fortran)
target_link_libraries(imsproc PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)

#target_compile_options(imsproc PRIVATE  -r8 -init=huge -cpp, -DCCPP=.false. )

# Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}>)

set( IMSPROC_INCLUDE_DIRS  ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( IMSPROC_LIBRARIES imsproc )

# calcfIMS.exe executable
#-------------------------------------------------------------------------------
set ( exe_files ./sorc/driver_fIMS.f90 )

ecbuild_add_executable( TARGET   calcfIMS.exe
                        SOURCES  ${exe_files}
                        LIBS     imsproc
                        LINKER_LANGUAGE ${IMSPROC_LINKER_LANGUAGE}
                       )


################################################################################
# Finalise configuration
################################################################################

ecbuild_install_project( NAME imsproc )
ecbuild_print_summary()
