cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(HK-IDLISOAnalysis VERSION 1.0)

set(CMAKE_CXX_STADARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
cmessage(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")
endif()

message(STATUS "CMAKE_INSTALL_PREFIX: \"${CMAKE_INSTALL_PREFIX}\"")

################ ROOT #################
find_package(ROOT REQUIRED COMPONENTS Tree MathCore Hist RIO)
#include("${ROOT_USE_FILE}")
#separate_arguments(ROOT_EXE_LINKER_FLAGS)

################ WCSim ################
if(NOT USE_WCSIM)
        set(USE_WCSIM 0)
elseif(USE_WCSIM AND NOT DEFINED ENV{WCSIMDIR})
        set(USE_WCSIM 0)
	message(WARNING "No WCSIMDIR path defined. Will build without WCSim.")
else()
	set(USE_WCSIM 1)
endif() 

if(USE_WCSIM)
	add_subdirectory(TreeConverter)
	include_directories(TreeConverter)
endif()

############### opticalAnalysis ############

add_executable(opticalAnalysis analysis.cc)

#target_link_libraries(opticalAnalysis PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})
target_link_libraries(opticalAnalysis PUBLIC ROOT::Tree ROOT::MathCore ROOT::Hist ROOT::RIO)

############### INSTALL ################

install(TARGETS opticalAnalysis DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)