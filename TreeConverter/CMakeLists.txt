if(USE_WCSIM)
	
	if ( NOT DEFINED ENV{WCSIMDIR} )
	  message (FATAL_ERROR "$WCSIMDIR is not defined, please set up WCSim first.")
	else()
  	  message(STATUS "Using WCSim installed at $ENV{WCSIMDIR}")
  	  set(CMAKE_WCSIMDIR $ENV{WCSIMDIR})
	endif()

	find_package(WCSim)

	set(WCSIM_CXX_FLAGS "-I${CMAKE_WCSIMDIR}/src;-I${CMAKE_WCSIMDIR}/include")
	set(WCSIM_LIBDIR "${CMAKE_WCSIMDIR}")

	list(APPEND WCSIM_LIBS
  	  WCSimRoot
	)	  

	link_directories(${WCSIM_LIBDIR})
	
	add_executable(TreeConverter TreeConverter.cc)
	target_link_libraries(TreeConverter ROOT::Tree ${WCSIM_LIBS})
	target_include_directories(TreeConverter PUBLIC ${CMAKE_SOURCE_DIR}/utils)
	target_compile_options(TreeConverter PUBLIC ${WCSIM_CXX_FLAGS})
endif()

