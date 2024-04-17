cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(HK-IDLISOAnalysis VERSION 1.0)

set(CMAKE_CXX_STADARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")
endif()

cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

################ ROOT #################

if(NOT DEFINED ENV{ROOTSYS})
        find_program(root_exists "root-config")
       	if (root_exists)
	   cmessage(STATUS "$ROOTSYS is not defined, but root-config exists.")
       	else()
	   cmessage(FATAL_ERROR "$ROOTSYS is not defined, please set up ROOT first.")
	endif()
else()
	cmessage(STATUS "Using ROOT installed at $ENV{ROOTSYS}")
	set(CMAKE_ROOTSYS $ENV{ROOTSYS})
endif()

find_package(ROOT REQUIRED COMPONENTS Tree MathCore Hist RIO Physics)


################ WCSim ################
if(NOT USE_WCSIM)
        set(USE_WCSIM 0)
elseif(USE_WCSIM AND NOT DEFINED ENV{WCSIMDIR})
        set(USE_WCSIM 0)
	cmessage(WARNING "No WCSIMDIR path defined. Will build without WCSim.")
else()
	set(USE_WCSIM 1)
endif() 


############### INSTALL ################

add_subdirectory(app)
add_subdirectory(src)

configure_file(cmake/setup.sh.in 
	"${CMAKE_SOURCE_DIR}/setup.sh" @ONLY)

unset(USE_WCSIM CACHE)
