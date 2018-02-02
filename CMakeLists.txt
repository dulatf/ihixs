


cmake_minimum_required(VERSION 2.6)
project(ihixs)


add_definitions(-Wno-deprecated-declarations) #LHAPDF6 uses features that are deprecated in gcc5


#option(with_google_tests "false")


find_package(OpenMP)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set (CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -DGTEST_HAS_PTHREAD=0 -pthread -Wno-deprecated -Wno-unused-local-typedefs -Wno-sign-compare -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")

find_package(Boost 1.36.0)
if(Boost_FOUND)
message("Boost headers found at : ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})
else()
message("Boost not found, we will trust the user defined: ${BOOST_DIR_USER}")
endif()
if(NOT LHAPDF_DIR)
  message("-- Trying to locate LHAPDF6")
  execute_process(COMMAND lhapdf-config --incdir OUTPUT_VARIABLE LHAPDF_INC_DIR)
  execute_process(COMMAND lhapdf-config --libdir OUTPUT_VARIABLE LHAPDF_LIB_DIR)
  find_library(LHAPDF_LIB_NAMES LHAPDF libLHAPDF HINTS ${LHAPDF_LIB_DIR})
  if(NOT LHAPDF_INC_DIR OR NOT LHAPDF_LIB_DIR OR NOT LHAPDF_LIB_NAMES)
    message("lhapdf-config not found, compiling will probably not work!")
  endif()
message("LHAPDF_INC_DIR : ${LHAPDF_INC_DIR}")
message("LHAPDF_LIB_NAMES: ${LHAPDF_LIB_NAMES}")
else()
  find_path(LHAPDF_INC_DIR LHAPDF/LHAPDF.h HINTS ${LHAPDF_DIR} ${LHAPDF_DIR}/include)
  find_library(LHAPDF_LIB_NAMES LHAPDF libLHAPDF HINTS ${LHAPDF_DIR} ${LHAPDF_DIR}/lib)
  message("LHAPDF_DIR : ${LHAPDF_DIR}")
  message("LHAPDF_INC_DIR : ${LHAPDF_INC_DIR}")
  message("LHAPDF_LIB_NAMES : ${LHAPDF_LIB_NAMES}")
endif()

if(NOT LHAPDF_INC_DIR OR NOT LHAPDF_LIB_NAMES)
  message(FATAL_ERROR "LHAPDF6 not found, aborting...")
endif()

message("-- Trying to locate Cuba")
find_path(CUBA_INC_DIR cuba.h HINTS ${CUBA_DIR_USER} ${CUBA_DIR_USER}/include)
find_library(CUBA_LIB_NAMES libcuba.a libcuba cuba HINTS ${CUBA_DIR_USER} ${CUBA_DIR_USER}/lib)
if(NOT CUBA_INC_DIR OR NOT CUBA_LIB_NAMES)
  message(FATAL_ERROR "Cuba not found, aborting...")
endif()
message("CUBA_INC_DIR : ${CUBA_INC_DIR}")
message("CUBA_LIB_NAMES: ${CUBA_LIB_NAMES}")

include_directories(${PROJECT_SOURCE_DIR}
                    ${PROJECT_SOURCE_DIR}/src/higgs/effective_theory/
                    ${PROJECT_SOURCE_DIR}/src/higgs/electroweak_corrections/
                    ${PROJECT_SOURCE_DIR}/src/higgs/exact_qcd_corrections/
                    ${PROJECT_SOURCE_DIR}/src/higgs/mt_expansion/
                    ${PROJECT_SOURCE_DIR}/src/core/
                    ${PROJECT_SOURCE_DIR}/src/tools/

                    ${PROJECT_SOURCE_DIR}/src/models/
                    ${PROJECT_SOURCE_DIR}/src/chaplin/
                    ${PROJECT_SOURCE_DIR}/src/gtest-1.7.0/include

                    ${PROJECT_SOURCE_DIR}/src/resummation/
                    ${PROJECT_SOURCE_DIR}/src/scet_resummation/

                    ${LHAPDF_INC_DIR}
                    ${BOOST_DIR_USER}
                    ${CUBA_INC_DIR}

		    )

add_subdirectory(src/chaplin)
add_subdirectory(src/higgs)
add_subdirectory(src/tools)
add_subdirectory(src/models)

if (DEFINED with_google_tests)
    add_subdirectory(src/tests)
    add_subdirectory(src/gtest-1.7.0)
    set(GTEST_LIB_NAME gtest)
    message("GTESTS are switched ON")
endif (DEFINED with_google_tests)




add_executable(ihixs
            ${PROJECT_SOURCE_DIR}/ihixs.cpp
            #${PROJECT_SOURCE_DIR}/src/higgs/ggf_manager/ggf_manager.cpp
            )


target_link_libraries(ihixs
cppchaplin
higgs_lib
models
tools
${CUBA_LIB_NAMES}
${LHAPDF_LIB_NAMES}
${GTEST_LIB_NAME}
)



file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/grids)

file(COPY ${PROJECT_SOURCE_DIR}/README.md
    DESTINATION ${CMAKE_BINARY_DIR}/ )

file(COPY ${PROJECT_SOURCE_DIR}/runcard/default.card
    DESTINATION ${CMAKE_BINARY_DIR}/)
file(COPY ${PROJECT_SOURCE_DIR}/runcard/pheno.card
DESTINATION ${CMAKE_BINARY_DIR}/)

set(CMAKE_BUILD_TYPE Release)
