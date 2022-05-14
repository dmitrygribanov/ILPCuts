cmake_minimum_required(VERSION 3.16)

project(ILPCuts VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include(InstallRequiredSystemLibraries)
#include(GNUInstallDirs)
set(ILOCPLEX_PATH /opt/ibm/ILOG/CPLEX_Studio201)
set(ILOCPLEX_INCLUDE_PATH ${ILOCPLEX_PATH}/cplex/include ${ILOCPLEX_PATH}/concert/include)
set(ILOCPLEX_INCLUDE_LIB_PATH ${ILOCPLEX_PATH}/cplex/lib/x86-64_linux/static_pic ${ILOCPLEX_PATH}/concert/lib/x86-64_linux/static_pic)

add_executable(CPlexConsole CPlex-console.cpp)
add_executable(CPlexDummy dummy_main.cpp)

target_include_directories(CPlexConsole PUBLIC "${ILOCPLEX_INCLUDE_PATH}")

target_link_directories(CPlexConsole PUBLIC "${ILOCPLEX_INCLUDE_LIB_PATH}")

target_link_libraries(CPlexConsole PUBLIC -lgmp -lgmpxx -lm -lpthread -lilocplex -lconcert -lcplex ${CMAKE_DL_LIBS})
target_link_libraries(CPlexDummy PUBLIC -lflint -lmpfr -lgmp -lgmpxx -lpthread -lntl ${CMAKE_DL_LIBS})

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(CPlexConsole PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(CPlexDummy PUBLIC OpenMP::OpenMP_CXX)
endif()