cmake_minimum_required(VERSION 3.15)

project(land-IMS_proc
        VERSION 1.0.0
        LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel)$")
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} is not supported.")
endif()

# Dependencies
#find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED COMPONENTS Fortran)

# Get list of files
list(APPEND land-IMS_proc_src driver_fIMS.f90  IMSaggregate_mod.f90 )

add_executable(calcfIMS ${land-IMS_proc_src})
target_link_libraries(calcfIMS NetCDF::NetCDF_Fortran)


