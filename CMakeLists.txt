PROJECT(LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(LRSpline_VERSION_MAJOR 0)
SET(LRSpline_VERSION_MINOR 4)
SET(LRSpline_VERSION_PATCH 0)
SET(LRSpline_VERSION ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR}.${LRSpline_VERSION_PATCH})

# Generate header with version info
CONFIGURE_FILE(include/LRSpline/LRSpline_version.h.in
               ${PROJECT_SOURCE_DIR}/include/LRSpline/LRSpline_version.h @ONLY)
CONFIGURE_FILE(test/fixupdebs.sh.in
               ${PROJECT_SOURCE_DIR}/test/fixupdebs.sh @ONLY)

ENABLE_LANGUAGE(CXX)

INCLUDE(CheckCXXCompilerFlag)

ADD_DEFINITIONS(-DTIME_LRSPLINE)

IF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
# check if compiler supports c++-0x
  CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
  IF(HAVE_0x)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
  ELSE(HAVE_0x)
    MESSAGE("A compiler with c++-0x support is needed")
    EXIT(1)
  ENDIF(HAVE_0x)
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES GNU)

ENABLE_TESTING()

# Required packages
FIND_PACKAGE(GoTrivariate REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(Boost REQUIRED)

# Required libraries
SET(DEPLIBS ${GoTrivariate_LIBRARIES}
            ${GoTools_LIBRARIES}
            ${BOOST_LIBRARIES})

# Required include directories
SET(INCLUDES ${PROJECT_SOURCE_DIR}/include
             ${GoTools_INCLUDE_DIRS} 
             ${GoTrivariate_INCLUDE_DIRS} 
             ${BOOST_INCLUDES})

INCLUDE_DIRECTORIES(${INCLUDES})

SET(EXECUTABLE_OUTPUT_PATH bin)
SET(LIBRARY_OUTPUT_PATH lib)

IF(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-parentheses")
ENDIF(NOT WIN32)

# Make the LRSpline library
FILE(GLOB_RECURSE LRSPLINE_SRCS ${PROJECT_SOURCE_DIR}/src/*.cpp)
ADD_LIBRARY(LRSpline ${LRSPLINE_SRCS})
TARGET_LINK_LIBRARIES(LRSpline ${DEPLIBS})
SET_TARGET_PROPERTIES(LRSpline PROPERTIES  VERSION ${LRSpline_VERSION}
                      SOVERSION ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR})

# Make some Apps
ADD_EXECUTABLE(TensorComparison ${PROJECT_SOURCE_DIR}/Apps/TensorComparison.cpp)
TARGET_LINK_LIBRARIES(TensorComparison LRSpline ${DEPLIBS})

ADD_EXECUTABLE(PartitionOfUnityTest ${PROJECT_SOURCE_DIR}/Apps/PartitionOfUnityTest.cpp)
TARGET_LINK_LIBRARIES(PartitionOfUnityTest LRSpline ${DEPLIBS})

ADD_EXECUTABLE(RefinementUnchanged ${PROJECT_SOURCE_DIR}/Apps/RefinementUnchanged.cpp)
TARGET_LINK_LIBRARIES(RefinementUnchanged LRSpline ${DEPLIBS})

ADD_EXECUTABLE(StresstestRefinement ${PROJECT_SOURCE_DIR}/Apps/StresstestRefinement.cpp)
TARGET_LINK_LIBRARIES(StresstestRefinement LRSpline ${DEPLIBS})

ADD_EXECUTABLE(StresstestEvaluation ${PROJECT_SOURCE_DIR}/Apps/StresstestEvaluation.cpp)
TARGET_LINK_LIBRARIES(StresstestEvaluation LRSpline ${DEPLIBS})

ADD_EXECUTABLE(TestReadWrite ${PROJECT_SOURCE_DIR}/Apps/TestReadWrite.cpp)
TARGET_LINK_LIBRARIES(TestReadWrite LRSpline ${DEPLIBS})

ADD_EXECUTABLE(Diagonal ${PROJECT_SOURCE_DIR}/Apps/Diagonal.cpp)
TARGET_LINK_LIBRARIES(Diagonal LRSpline ${DEPLIBS})

ADD_EXECUTABLE(MakeWaterfall ${PROJECT_SOURCE_DIR}/Apps/MakeWaterfall.cpp)
TARGET_LINK_LIBRARIES(MakeWaterfall LRSpline ${DEPLIBS})

ADD_EXECUTABLE(LinearIndep ${PROJECT_SOURCE_DIR}/Apps/LinearIndep.cpp)
TARGET_LINK_LIBRARIES(LinearIndep LRSpline ${DEPLIBS})

ADD_EXECUTABLE(drawLRmesh ${PROJECT_SOURCE_DIR}/Apps/drawLRmesh.cpp)
TARGET_LINK_LIBRARIES(drawLRmesh LRSpline ${DEPLIBS})

ADD_EXECUTABLE(TopologyRefinement ${PROJECT_SOURCE_DIR}/Apps/TopologyRefinement.cpp)
TARGET_LINK_LIBRARIES(TopologyRefinement LRSpline ${DEPLIBS})

ADD_EXECUTABLE(BezierExtract ${PROJECT_SOURCE_DIR}/Apps/BezierExtract.cpp)
TARGET_LINK_LIBRARIES(BezierExtract LRSpline ${DEPLIBS})

ADD_EXECUTABLE(refine ${PROJECT_SOURCE_DIR}/Apps/refine.cpp)
TARGET_LINK_LIBRARIES(refine LRSpline ${DEPLIBS})

# # Regression tests
FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/RefinementUnchanged/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/RefinementUnchanged" "${TESTFILE}")
ENDFOREACH()

FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/TensorComparison/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/TensorComparison" "${TESTFILE}")
ENDFOREACH()

FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/Diagonal/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/Diagonal" "${TESTFILE}")
ENDFOREACH()

FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/LinearIndep/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/LinearIndep" "${TESTFILE}")
ENDFOREACH()

FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/TestReadWrite/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/TestReadWrite/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/TestReadWrite" "${TESTFILE}")
ENDFOREACH()

FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/TopologyRefinement/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/TopologyRefinement" "${TESTFILE}")
ENDFOREACH()

# 'install' target
IF(WIN32)
  # TODO
ELSE(WIN32)
  # lib
  INSTALL(TARGETS LRSpline DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
          COMPONENT liblrspline)

  # headers
  FILE(GLOB LRSPLINE_HEADERS include/LRSpline/Basisfunction.h
                             include/LRSpline/Element.h
                             include/LRSpline/Meshline.h
                             include/LRSpline/LRSpline_version.h
                             include/LRSpline/LRSpline.h
                             include/LRSpline/LRSplineSurface.h
                             include/LRSpline/LRSplineVolume.h
                             include/LRSpline/HashSet.h
                             include/LRSpline/MeshRectangle.h)
  INSTALL(FILES ${LRSPLINE_HEADERS}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/include/LRSpline
                COMPONENT liblrspline-dev)
  INSTALL(FILES cmake/Modules/FindLRSpline.cmake
          DESTINATION ${CMAKE_ROOT}/Modules COMPONENT liblrspline-dev)
ENDIF(WIN32)

# Packaging
EXECUTE_PROCESS(COMMAND lsb_release "-sc" OUTPUT_VARIABLE CODENAME)
STRING(REGEX REPLACE "\n" "" SYSTEM_CODENAME ${CODENAME})
EXECUTE_PROCESS(COMMAND lsb_release "-si" OUTPUT_VARIABLE VENDOR)
STRING(REGEX REPLACE "\n" "" SYSTEM_VENDOR ${VENDOR})

IF (SYSTEM_VENDOR STREQUAL "Ubuntu" OR SYSTEM_VENDOR STREQUAL "Debian")
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE ARCHITECTURE)
  STRING(REGEX REPLACE "\n" "" CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${ARCHITECTURE})
  SET(SYSTEM_ARCHITECTURE ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
  SET(CPACK_GENERATOR "DEB")
  SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Arne Morten Kvarving <arne.morten.kvarving@sintef.no>")
  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libgotools-core")
  SET(CPACK_DEB_COMPONENT_INSTALL ON)
ENDIF (SYSTEM_VENDOR STREQUAL "Ubuntu" OR SYSTEM_VENDOR STREQUAL "Debian")

SET(CPACK_PACKAGE_NAME "liblrspline")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LRSplines - the locally-refined spline library")
SET(CPACK_PACKAGE_VERSION "${LRSpline_VERSION}")
SET(CPACK_PACKAGE_VERSION_MAJOR ${LRSpline_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${LRSpline_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${LRSpline_VERSION_PATCH})
SET(CPACK_PACKAGE_FILE_NAME "liblrspline_${CPACK_PACKAGE_VERSION}_${SYSTEM_ARCHITECTURE}-${SYSTEM_CODENAME}")
SET(CPACK_COMPONENTS_ALL liblrspline liblrspline-dev)
#SET(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/COPYING)
SET(CPACK_STRIP_FILES "lib/libLRSpline.so")

INCLUDE(CPack)

ADD_CUSTOM_TARGET(ubuntu 
                  COMMAND make package
                  COMMAND ${PROJECT_SOURCE_DIR}/test/fixupdebs.sh
                  COMMAND echo "All done. Packages are in UbuntuDebs/"
                  COMMENT "Generating and fixing up Ubuntu packages" VERBATIM)
