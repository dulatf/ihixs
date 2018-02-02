add_library(higgs_lib
    ../core/channel.cpp
    ../core/inclusive_process.cpp
    ../core/luminosity_integrals.cpp
    ../core/sigma_term.cpp
    ../core/input_parameters.cpp
    effective_theory/higgs_eft.cpp
    effective_theory/eft_gg_delta.cpp
    effective_theory/eft_gg_plus.cpp
    effective_theory/eft_gg_reg.cpp
    effective_theory/eft_gg_reg_nnlo.cpp
    effective_theory/eft_gg_reg_n3lo.cpp
    effective_theory/eft_qg_reg.cpp
    effective_theory/eft_qqbar_reg.cpp
    effective_theory/eft_qq_reg.cpp
    effective_theory/eft_q1q2_reg.cpp
    effective_theory/wilson_coefficients.cpp
    effective_theory/inclusive_higgs_eft.cpp
    effective_theory/inclusive_higgs_eft_fast.cpp
    #
    electroweak_corrections/gluon_fusion_ew_coefficients.cpp
    #
    exact_qcd_corrections/nlo_exact_matrix_elements.cpp
    exact_qcd_corrections/inclusive_higgs_exact.cpp
    #
    ggf_manager/ggf_manager.cpp
    #
    mt_expansion/mt_expansion.cpp
    mt_expansion/mt_expansion_gg.cpp
    mt_expansion/mt_expansion_qg.cpp
    mt_expansion/inclusive_mt_expansion.cpp
    ../resummation/threshold_resummation.cpp
    ../resummation/gauss_grid.cpp
    ../resummation/psi.cpp
    ../resummation/wpsipg.cpp
    #
    ../scet_resummation/scetter.cpp
)
target_link_libraries(higgs_lib cppchaplin)
if (DEFINED with_fortran_tests)
add_library(../resummation/fortran_psi wpsipg.f)
endif (DEFINED with_fortran_tests)




