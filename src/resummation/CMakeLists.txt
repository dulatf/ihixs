add_library(resummation
                threshold_resummation.cpp
                gauss_grid.cpp
                psi.cpp
                wpsipg.cpp
                )
if (DEFINED with_fortran_tests)
    add_library(fortran_psi wpsipg.f)
endif (DEFINED with_fortran_tests)
