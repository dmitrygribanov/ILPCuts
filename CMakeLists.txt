cmake_minimum_required(VERSION 3.16)

project(ILPCuts VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include(InstallRequiredSystemLibraries)
#include(GNUInstallDirs)
set(ILOCPLEX_PATH /home/istikhobaev/cplex_installation_folder)
set(ILOCPLEX_INCLUDE_PATH ${ILOCPLEX_PATH}/cplex/include ${ILOCPLEX_PATH}/concert/include)
set(ILOCPLEX_INCLUDE_LIB_PATH ${ILOCPLEX_PATH}/cplex/lib/x86-64_linux/static_pic ${ILOCPLEX_PATH}/concert/lib/x86-64_linux/static_pic)

set(BOOST_PATH /home/istikhobaev/boost)
set(BOOST_INCLUDE_PATH ${BOOST_PATH}/include)
set(BOOST_INCLUDE_LIB_PATH ${BOOST_PATH}/lib)

set(FLINT_PATH /home/istikhobaev/old_env/flint)
set(FLINT_INCLUDE_PATH ${FLINT_PATH}/include)
set(FLINT_INCLUDE_LIB_PATH ${FLINT_PATH}/lib)

set(GMP_PATH /home/istikhobaev/old_env/gmp)
set(GMP_INCLUDE_PATH ${GMP_PATH}/include)
set(GMP_INCLUDE_LIB_PATH ${GMP_PATH}/lib)

set(MPFR_PATH /home/istikhobaev/old_env/mpfr_lib)
set(MPFR_INCLUDE_PATH ${MPFR_PATH}/include)
set(MPFR_INCLUDE_LIB_PATH ${MPFR_PATH}/lib)

set(NTL_PATH /home/istikhobaev/ntl)
set(NTL_INCLUDE_PATH ${NTL_PATH}/include)
set(NTL_INCLUDE_LIB_PATH ${NTL_PATH}/lib)

add_executable(CPlexConsole CPlex-console.cpp)
add_executable(CPlexDummy dummy_main.cpp)
add_executable(CPlexNodeExpsCPP CPlex-node-experiments-c++.cpp)
add_executable(CPlexNodeExpsC CPlex-node-experiments-C.cpp)
add_executable(CPlexCutsPower CPlex-cuts-power.cpp)
add_executable(Tests Tests.cpp)
add_executable(CPlexCutsApplier CPlex-cuts-appplier.cpp)

target_include_directories(CPlexConsole PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(CPlexNodeExpsCPP PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(CPlexNodeExpsC PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(CPlexCutsPower PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(Tests PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(CPlexCutsApplier PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")
target_include_directories(CPlexDummy PUBLIC "${ILOCPLEX_INCLUDE_PATH}" "${BOOST_INCLUDE_PATH}" "${FLINT_INCLUDE_PATH}" "${GMP_INCLUDE_PATH}" "${MPFR_INCLUDE_PATH}" "${NTL_INCLUDE_PATH}")

target_link_directories(CPlexConsole PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(CPlexNodeExpsCPP PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(CPlexNodeExpsC PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(CPlexCutsPower PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(Tests PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(CPlexCutsApplier PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")
target_link_directories(CPlexDummy PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}" "${BOOST_INCLUDE_LIB_PATH}" "${FLINT_INCLUDE_LIB_PATH}" "${GMP_INCLUDE_LIB_PATH}" "${MPFR_INCLUDE_LIB_PATH}" "${NTL_INCLUDE_LIB_PATH}")

#removed deps.: -lgivaro -lblas -llinbox
target_link_libraries(CPlexConsole PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexDummy PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexNodeExpsCPP PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexNodeExpsC PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexCutsPower PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(Tests PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexCutsApplier PUBLIC -lflint -lgmp -lgmpxx -lm -lpthread -lntl -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(CPlexConsole PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexDummy PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexNodeExpsCPP PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexNodeExpsC PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexCutsPower PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(Tests PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexCutsApplier PUBLIC OpenMP::OpenMP_CXX)
endif()