PROJECT(LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(LRSpline_VERSION_MAJOR 0)
SET(LRSpline_VERSION_MINOR 5)
SET(LRSpline_VERSION_PATCH 0)
SET(LRSpline_VERSION ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR}.${LRSpline_VERSION_PATCH})
SET(TIME_LRSPLINE 1)

# Generate header with version info
CONFIGURE_FILE(include/LRSpline/LRSpline_version.h.in
               ${PROJECT_SOURCE_DIR}/include/LRSpline/LRSpline_version.h @ONLY)

ENABLE_LANGUAGE(CXX)

INCLUDE(CheckCXXCompilerFlag)
INCLUDE(cmake/Modules/UseMultiArch.cmake)

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



# Check that the required C++11 features are available
CHECK_CXX_SOURCE_COMPILES("
                          #include <map>
                          #include <list>
                          int main(void)
                          { 
                            std::map<int, std::list<int> > loop;
                            loop.insert(std::make_pair(1,std::list<int>(1)));
                            std::map<int, std::list<int >>::iterator it = loop.begin();
                            for (int a : it->second)
                              ;
                          }
                          " COMPILER_SATISFIES_CXX11)
IF (NOT COMPILER_SATISFIES_CXX11)
  MESSAGE(FATAL_ERROR "Compiler do not provide C++-11 features required, bailing")
ENDIF(NOT COMPILER_SATISFIES_CXX11)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/boost_1_55_0/boost_1_55_0/stage/lib -I/boost_1_55_0/boost_1_55_0")

ENABLE_TESTING()

# Required packages
#    - none -

# Required include directories
SET(INCLUDES ${PROJECT_SOURCE_DIR}/include)

# Optional packages
FIND_PACKAGE(GoTools)
FIND_PACKAGE(Boost)


IF(GOTOOLS_FOUND)
  # must have all of GoTools or none of it
  FIND_PACKAGE(GoTrivariate REQUIRED)
  OPTION(HAS_GOTOOLS "Link library with GoTools" ON)
  SET(INCLUDES ${INCLUDES}
               ${GoTools_INCLUDE_DIRS} 
               ${GoTrivariate_INCLUDE_DIRS})
  SET(DEPLIBS ${DEPLIBS}
              ${GoTrivariate_LIBRARIES}
              ${GoTools_LIBRARIES})
  ADD_DEFINITIONS(-DHAS_GOTOOLS)
ELSE(GOTOOLS_FOUND)
  OPTION(HAS_GOTOOLS "Link library with GoTools" OFF)
ENDIF(GOTOOLS_FOUND)

IF(Boost_FOUND)
  SET(DEPLIBS  ${DEPLIBS}  ${BOOST_LIBRARIES})
  SET(INCLUDES ${INCLUDES} ${BOOST_INCLUDES} )
  ADD_DEFINITIONS(-DHAS_BOOST)
  OPTION(HAS_BOOST "Link library with Boost" ON)
ELSE(Boost_FOUND)
  OPTION(HAS_BOOST "Link library with Boost" OFF)
ENDIF(Boost_FOUND)


IF(TIME_LRSPLINE)
  ADD_DEFINITIONS(-DTIME_LRSPLINE)
ENDIF(TIME_LRSPLINE)


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
ADD_LIBRARY(LRSpline STATIC ${LRSPLINE_SRCS})
TARGET_LINK_LIBRARIES(LRSpline ${DEPLIBS})
#SET_TARGET_PROPERTIES(LRSpline PROPERTIES  VERSION ${LRSpline_VERSION}
#                      SOVERSION ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR})

# Make some Apps
ADD_EXECUTABLE(reduceContinuity ${PROJECT_SOURCE_DIR}/Apps/reduceContinuity.cpp)
TARGET_LINK_LIBRARIES(reduceContinuity LRSpline ${DEPLIBS})

ADD_EXECUTABLE(PartitionOfUnityTest ${PROJECT_SOURCE_DIR}/Apps/PartitionOfUnityTest.cpp)
TARGET_LINK_LIBRARIES(PartitionOfUnityTest LRSpline ${DEPLIBS})

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

ADD_EXECUTABLE(RefinementUnchanged ${PROJECT_SOURCE_DIR}/Apps/RefinementUnchanged.cpp)
TARGET_LINK_LIBRARIES(RefinementUnchanged LRSpline ${DEPLIBS})

IF(HAS_GOTOOLS)
  ADD_EXECUTABLE(TensorComparison ${PROJECT_SOURCE_DIR}/Apps/TensorComparison.cpp)
  TARGET_LINK_LIBRARIES(TensorComparison LRSpline ${DEPLIBS})
ENDIF(HAS_GOTOOLS)

ADD_EXECUTABLE(TopologyRefinement ${PROJECT_SOURCE_DIR}/Apps/TopologyRefinement.cpp)
TARGET_LINK_LIBRARIES(TopologyRefinement LRSpline ${DEPLIBS})

ADD_EXECUTABLE(LinearIndep ${PROJECT_SOURCE_DIR}/Apps/LinearIndep.cpp)
TARGET_LINK_LIBRARIES(LinearIndep LRSpline ${DEPLIBS})

ADD_EXECUTABLE(drawLRmesh ${PROJECT_SOURCE_DIR}/Apps/drawLRmesh.cpp)
TARGET_LINK_LIBRARIES(drawLRmesh LRSpline ${DEPLIBS})

ADD_EXECUTABLE(BezierExtract ${PROJECT_SOURCE_DIR}/Apps/BezierExtract.cpp)
TARGET_LINK_LIBRARIES(BezierExtract LRSpline ${DEPLIBS})

ADD_EXECUTABLE(refine ${PROJECT_SOURCE_DIR}/Apps/refine.cpp)
TARGET_LINK_LIBRARIES(refine LRSpline ${DEPLIBS})

# # Regression tests
FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/RefinementUnchanged/*.reg")
FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/RefinementUnchanged" "${TESTFILE}")
ENDFOREACH()

IF(HAS_GOTOOLS)
  FILE(GLOB_RECURSE REGRESESSION_TESTFILES "${PROJECT_SOURCE_DIR}/test/TensorComparison/*.reg")
  FOREACH(TESTFILE ${REGRESESSION_TESTFILES})
    ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/test/regtest.sh "${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/TensorComparison" "${TESTFILE}")
  ENDFOREACH()
ENDIF(HAS_GOTOOLS)

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
  install(TARGET LRSplines DESTINATION LRSplines)
  install(FILES ${INCLUDES} DESTINATION LRSplines/include)
  # TODO
ELSE(WIN32)
  # lib
  INSTALL(TARGETS LRSpline DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT liblrspline)

  # For generating the doxy
  ADD_CUSTOM_TARGET(doc doxygen 
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    COMMENT "Generating API documentation" VERBATIM)
  INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)

  INSTALL(DIRECTORY doc/html/ DESTINATION share/doc/liblrsplines1/html)

  # headers
  FILE(GLOB LRSPLINE_HEADERS include/LRSpline/Basisfunction.h
                             include/LRSpline/Element.h
                             include/LRSpline/Meshline.h
                             include/LRSpline/LRSpline_version.h
                             include/LRSpline/LRSpline.h
                             include/LRSpline/LRSplineSurface.h
                             include/LRSpline/LRSplineVolume.h
                             include/LRSpline/Streamable.h
                             include/LRSpline/HashSet.h
                             include/LRSpline/MeshRectangle.h)
  INSTALL(FILES ${LRSPLINE_HEADERS}
                DESTINATION include/LRSpline
                COMPONENT liblrspline-dev)
  INSTALL(FILES cmake/Modules/LRSplineConfig.cmake
    DESTINATION lib/LRSpline COMPONENT liblrspline-dev)
ENDIF(WIN32)

include(CPack)
