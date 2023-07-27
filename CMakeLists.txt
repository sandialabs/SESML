#--------------------------------------------------------------------#
#  Copyright (2013) Sandia Corporation. Under the terms of Contract  #
#  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    #
#  retains certain rights in this software.                          #
#                                                                    #
#  For full terms of distribution see the accompanying LICENSE file. #
#  To distribute this file outside of SESML, substitue the full text #
#  of the LICENSE file for this notice.                              #
#--------------------------------------------------------------------#

cmake_minimum_required(VERSION 2.8)

#
# CMakeLists.txt file for the Sandia Equation of State Model Library
#
# C is included for Threads detection for the GTest library
#
project(SESML CXX C)

enable_testing()

#
# Source files
#
set(sesml_sources
DFMin.C FiniteDiffDeriv.C PhaseBoundaryInfo.C Region.C 
EOSData.C EOSFactory.C EOSModel.C EOSParam.C
IAPWS-IF97.C IAPWStransport.C IAPWS95.C
)

#
# Header files
#
set(sesml_headers
DFMin.H FiniteDiffDeriv.H PhaseBoundaryInfo.H Region.H Vars.H
EOSData.H EOSFactory.H EOSModel.H EOSParam.H
IAPWS95.H IAPWS-IF97.H IAPWStransport.H
)

#
# Library
#
add_library(sesml ${sesml_sources})

#
# GTest stuff
#
include_directories(${SESML_GTEST_INCLUDE_DIR})
if(SESML_GTEST_IN_TREE)
  set(SESML_GTEST_LIB gtest)
  set(SESML_GTEST_MAIN gtest_main)
else()
  find_library(SESML_GTEST_LIB gtest ${SESML_GTEST_LIB_DIR})
  find_library(SESML_GTEST_MAIN gtest_main ${SESML_GTEST_LIB_DIR})
endif()
find_package(Threads)

#
# Unit tests
#
macro(makeunittest name)
add_executable(${name}_unittest ${name}_unittest.C)
target_link_libraries(${name}_unittest sesml ${SESML_GTEST_MAIN} ${SESML_GTEST_LIB} ${CMAKE_THREAD_LIBS_INIT})
add_test(NAME ${name}_unittest COMMAND ${name}_unittest)
endmacro(makeunittest)

makeunittest(DFMin)
makeunittest(EOSData)
makeunittest(EOSFactory)
makeunittest(EOSModel)
makeunittest(EOSParam)
makeunittest(FiniteDiffDeriv)
makeunittest(IAPWS-IF97)
makeunittest(IAPWS95)
makeunittest(IAPWStransport)
makeunittest(Region)

#
# install
#
install(FILES ${sesml_headers} DESTINATION include/sesml)
install(TARGETS sesml EXPORT sesml
        LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

#
# Source packaging
#
set(CPACK_SOURCE_GENERATOR "ZIP;TGZ")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "sesml-0.1")

set(CPACK_SOURCE_IGNORE_FILES ".*~" "[.]svn" "build")

include(CPack)
