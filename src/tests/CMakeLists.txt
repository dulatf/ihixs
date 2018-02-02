

add_executable(ihixs_nlo_exact_test
            ihixs_nlo_exact_test.cpp)
target_link_libraries(ihixs_nlo_exact_test
                        higgs_lib
                        models
                        cppchaplin
                        tools
                        gtest)


#-------------------------------------------------------------------------------

                    
add_executable(
            ihixs_as_series
            ihixs_as_series.cpp
            )
target_link_libraries(
                    ihixs_as_series
                    higgs_lib
                    tools
                    #cuba
                    ${CUBA_LIB_NAMES}
                    gtest
                    )
#-------------------------------------------------------------------------------

add_executable(ihixs_eft ihixs_eft.cpp)

target_link_libraries(ihixs_eft
    higgs_lib
    cppchaplin
    models
    tools
    #cuba
    ${CUBA_LIB_NAMES}
    gtest
)

#-------------------------------------------------------------------------------

add_executable(
        quark_evolution_test
        quark_evolution_test.cpp
)

target_link_libraries(quark_evolution_test
    models
    tools
    gtest)

#-------------------------------------------------------------------------------

add_executable(mt_expansion_test mt_expansion_test.cpp)
target_link_libraries(mt_expansion_test
    higgs_lib
    models
    tools
    cppchaplin
    #cuba
    ${CUBA_LIB_NAMES}
    #LHAPDF
    ${LHAPDF_LIB_NAMES}
    gtest
)




          
#-------------------------------------------------------------------------------

add_executable(eft_q1q2_n3lo_Lf_test eft_q1q2_n3lo_Lf_test.cpp)
target_link_libraries(eft_q1q2_n3lo_Lf_test
    #chaplin
    cppchaplin
    higgs_lib
    gtest)

#-------------------------------------------------------------------------------

              
                    


