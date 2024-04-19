if(USE_WCSIM)
	
	if ( NOT DEFINED ENV{WCSIMDIR} )
	  cmessage (FATAL_ERROR "$WCSIMDIR is not defined, please set up WCSim first.")
	else()
  	  cmessage(STATUS "Using WCSim installed at $ENV{WCSIMDIR}")
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
	target_include_directories(TreeConverter PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/utils)
	target_compile_options(TreeConverter PUBLIC ${WCSIM_CXX_FLAGS})

	install(TARGETS TreeConverter DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

endif()

add_executable(opticalAnalysis analysis.cc)
target_include_directories(opticalAnalysis PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(opticalAnalysis PUBLIC ROOT::Tree ROOT::MathCore ROOT::Hist ROOT::RIO ROOT::Physics HK-IDLISOAnalysis)

install(TARGETS opticalAnalysis DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
