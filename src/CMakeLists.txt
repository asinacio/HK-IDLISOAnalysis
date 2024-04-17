set(HEADERS 
    DataFilter.hh
    DataPMT.hh
    DataReader.hh
    DataRunInfo.hh
    OpticalModel.hh
)


add_library(HK-IDLISOAnalysis SHARED
    DataFilter.cc  
    DataPMT.cc     
    DataReader.cc  
    DataRunInfo.cc 
    OpticalModel.cc
)


set_target_properties(HK-IDLISOAnalysis PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_link_libraries(HK-IDLISOAnalysis PUBLIC ROOT::Tree ROOT::MathCore ROOT::Hist ROOT::RIO ROOT::Physics)
target_include_directories(HK-IDLISOAnalysis PUBLIC ${CMAKE_SOURCE_DIR}/src)
root_generate_dictionary(G__HKIDLISOAnalysis ${HEADERS} MODULE HK-IDLISOAnalysis LINKDEF LinkDef.h)

install(TARGETS HK-IDLISOAnalysis DESTINATION lib/ PUBLIC_HEADER DESTINATION include/)
install(TARGETS G__HKIDLISOAnalysis DESTINATION lib/)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libHK-IDLISOAnalysis_rdict.pcm
	      ${CMAKE_CURRENT_BINARY_DIR}/libHK-IDLISOAnalysis.rootmap
	      DESTINATION lib/)