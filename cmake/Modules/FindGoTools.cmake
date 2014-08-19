IF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
  SET(GoTools_FIND_QUIETLY TRUE)
ENDIF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)

FIND_PATH(GoTools_INCLUDE_DIRS
  NAMES GoTools/geometry/SplineSurface.h
  PATHS "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  /sima/libs/GoTools/include
)

FIND_LIBRARY(GoTools_LIBRARIES
  NAMES GoToolsCore
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
)

IF(GOTOOLS_HEADER)
  # Check for newer GoTools
  EXECUTE_PROCESS(COMMAND cat "${GoTools_INCLUDE_DIRS}/GoTools/geometry/GoTools.h" OUTPUT_VARIABLE GOTOOLS_HEADER)
  STRING(REGEX MATCH "GO_VERSION_MAJOR ([0-9]+)" GoTools_VERSION_MAJOR ${GOTOOLS_HEADER})
  IF(NOT GoTools_VERSION_MAJOR)
    EXECUTE_PROCESS(COMMAND cat "${GoTools_INCLUDE_DIRS}/GoTools/geometry/GoTools_version.h" OUTPUT_VARIABLE GOTOOLS_HEADER)
    STRING(REGEX MATCH "GO_VERSION_MAJOR ([0-9]+)" GoTools_VERSION_MAJOR ${GOTOOLS_HEADER})
  ENDIF(NOT GoTools_VERSION_MAJOR)
  STRING(REGEX REPLACE "GO_VERSION_MAJOR ([0-9]+)" "\\1" GoTools_VERSION_MAJOR "${GoTools_VERSION_MAJOR}")
  STRING(REGEX MATCH "GO_VERSION_MINOR ([0-9]+)" GoTools_VERSION_MINOR ${GOTOOLS_HEADER})
  STRING(REGEX REPLACE "GO_VERSION_MINOR ([0-9]+)" "\\1" GoTools_VERSION_MINOR "${GoTools_VERSION_MINOR}")
  STRING(REGEX MATCH "GO_VERSION_PATCH ([0-9]+)" GoTools_VERSION_PATCH ${GOTOOLS_HEADER})
  STRING(REGEX REPLACE "GO_VERSION_PATCH ([0-9]+)" "\\1" GoTools_VERSION_PATCH "${GoTools_VERSION_PATCH}")
  
  IF (GoTools_VERSION_MAJOR GREATER 2)
    INCLUDE(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
    IF(HAVE_0x)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
      SET(CMAKE_REQUIRED_FLAGS "-std=gnu++0x")
      CHECK_CXX_SOURCE_COMPILES (
        "#include <memory>
        int main(void) {
          std::shared_ptr<int> bar;
          }" HAVE_SHARED_PTR_0x)
    ENDIF(HAVE_0x)
    IF(NOT HAVE_SHARED_PTR_0x)
      # C++0x shared_ptr is not supported - check for Boost
      FIND_PACKAGE(Boost)
      IF(Boost_FOUND)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
        SET(GoTools_COMMON_INCLUDE_DIRS
            ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
      ELSE(Boost_FOUND)
          MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
      ENDIF(Boost_FOUND)
    ENDIF(NOT HAVE_SHARED_PTR_0x)
  
    # Quirk test: GoTools r10221 changed function signature without bumping
    IF (GoTools_VERSION_MAJOR EQUAL 4 AND GoTools_VERSION_MINOR EQUAL 0
                                      AND GoTools_VERSION_PATCH EQUAL 1)
      INCLUDE(CheckFunctionExists)
      SET(CMAKE_REQUIRED_FLAGS ${GoTools_CXX_FLAGS})
      SET(CMAKE_REQUIRED_LIBRARIES ${GoTools_LIBRARIES})
        SET(CMAKE_REQUIRED_INCLUDES ${GoTools_INCLUDE_DIRS})
        CHECK_CXX_SOURCE_COMPILES("
                                   #include <GoTools/geometry/LoopUtils.h>
                                int main(void)
                                { 
                                  std::vector<shared_ptr<Go::CurveOnSurface> > loop;
                                  Go::LoopUtils::paramIsCCW(loop, 1.0, 1.0)
                                }
                               " GoTools_QUIRK_paramIsCCW)
      ENDIF (GoTools_VERSION_MAJOR EQUAL 4 AND GoTools_VERSION_MINOR EQUAL 0
                                           AND GoTools_VERSION_PATCH EQUAL 1)
  ENDIF (GoTools_VERSION_MAJOR GREATER 2)
  
  INCLUDE(FindPackageHandleStandardArgs)
  IF(GoTools_LIBRARIES)
    find_package_handle_standard_args(GoTools DEFAULT_MSG
                                      GoTools_INCLUDE_DIRS GoTools_LIBRARIES)
  ENDIF(GoTools_LIBRARIES)
ENDIF(GOTOOLS_HEADER)
