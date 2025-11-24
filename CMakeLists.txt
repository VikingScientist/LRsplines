cmake_minimum_required(VERSION 3.22)

set(LRSpline_VERSION_MAJOR 1)
set(LRSpline_VERSION_MINOR 7)
set(LRSpline_VERSION_PATCH 0)
set(LRSpline_VERSION ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR}.${LRSpline_VERSION_PATCH})

project(LRSpline
        VERSION ${LRSpline_VERSION}
        LANGUAGES CXX
)

if(WIN32)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)
endif()

configure_file(include/LRSpline/LRSpline_version.h.in LRSpline_version.h @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in Doxyfile @ONLY)

add_library(LRSpline)
target_compile_features(LRSpline PUBLIC cxx_std_17)

set_target_properties(
  LRSpline
  PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY
      lib
    LIBRARY_OUTPUT_DIRECTORY
      lib
)

find_package(SISL QUIET)
find_package(GoToolsCore QUIET)
find_package(GoTrivariate QUIET)
if(TARGET GoToolsCore AND TARGET GoTrivariate)
  add_library(GoTools::GoToolsCore IMPORTED INTERFACE)
  string(REPLACE "." ";" GoVersion ${GoToolsCore_VERSION})
  list(GET GoVersion 0 GoToolsCore_VERSION_MAJOR)
  target_compile_definitions(GoTools::GoToolsCore
    INTERFACE
      GoTools_VERSION_MAJOR=${GoToolsCore_VERSION_MAJOR}
      GOTOOLS_HAS_BASISDERIVS_SF3=1
  )
  target_link_libraries(GoTools::GoToolsCore INTERFACE GoToolsCore)
  add_library(GoTools::GoTrivariate IMPORTED INTERFACE)
  target_link_libraries(
    GoTools::GoTrivariate
    INTERFACE
      GoTrivariate
      GoTools::GoToolsCore
  )
else()
  find_package(GoTools)
  find_package(GoTrivariate)
  if(GoTools_FOUND AND GoTrivariate_FOUND)
    find_package(GoTools)
    add_library(GoTools::GoToolsCore INTERFACE IMPORTED)
    set(GoTools_LIBRARY_DIRS ${GoTools_LIBRARIES})
    list(FILTER GoTools_LIBRARY_DIRS INCLUDE REGEX -L.*)
    list(TRANSFORM GoTools_LIBRARY_DIRS REPLACE "-L" "")
    list(FILTER GoTools_LIBRARIES EXCLUDE REGEX -L.*)
    target_link_libraries(GoTools::GoToolsCore INTERFACE ${GoTools_LIBRARIES})
    target_include_directories(GoTools::GoToolsCore INTERFACE ${GoTools_INCLUDE_DIRS})
    target_compile_definitions(GoTools::GoToolsCore INTERFACE ${GoTools_CXX_FLAGS} -DGoTools_VERSION_MAJOR=4)
    target_link_directories(GoTools::GoToolsCore INTERFACE ${GoTools_LIBRARY_DIRS})
    include(CheckTypeSize)
    set(CMAKE_EXTRA_INCLUDE_FILES "GoTools/geometry/SplineSurface.h")
    set(CMAKE_REQUIRED_FLAGS -std=c++11)
    check_type_size(Go::BasisDerivsSf3 HAS_BD3 LANGUAGE CXX)
    if(HAS_BD3)
      target_compile_definitions(GoTools::GoToolsCore INTERFACE GOTOOLS_HAS_BASISDERIVS_SF3=1)
    endif()

    add_library(GoTools::GoTrivariate INTERFACE IMPORTED)
    set(GoTrivariate_LIBRARY_DIRS ${GoTrivariate_LIBRARIES})
    list(FILTER GoTrivariate_LIBRARY_DIRS INCLUDE REGEX -L.*)
    list(TRANSFORM GoTrivariate_LIBRARY_DIRS REPLACE "-L" "")
    list(FILTER GoTrivariate_LIBRARIES EXCLUDE REGEX -L.*)
    target_include_directories(GoTools::GoTrivariate INTERFACE ${GoTrivariate_INCLUDE_DIRS})
    target_link_libraries(GoTools::GoTrivariate INTERFACE ${GoTrivariate_LIBRARIES} GoTools::GoToolsCore)
    target_link_directories(GoTools::GoTrivariate INTERFACE ${GoTrivariate_LIBRARY_DIRS})
  endif()
endif()

target_compile_definitions(LRSpline INTERFACE LRSPLINE_HAS_FUNC_TOLERANCE=1)

if(NOT TARGET GoTrivariate)
  message(STATUS "Compiling without GoTools")
endif()

find_package(Boost CONFIG)

if(NOT Boost_FOUND)
  message(STATUS "Compiling without Boost")
endif()

if(NOT WIN32)
  target_compile_definitions(LRSpline PRIVATE TIME_LRSPLINE=1)
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wall HAS_WALL)
if(HAS_WALL)
  target_compile_options(LRSpline PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
endif()

check_cxx_compiler_flag(-Wno-parentheses HAS_PARENTHESES)
if(HAS_PARENTHESES)
  target_compile_options(LRSpline PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-Wno-parentheses>)
endif()

target_sources(
  LRSpline
  PRIVATE
    src/Basisfunction.cpp
    src/Element.cpp
    src/LRSpline.cpp
    src/LRSplineSurface.cpp
    src/LRSplineVolume.cpp
    src/Meshline.cpp
    src/MeshRectangle.cpp
    src/Profiler.cpp
)

# TODO: Remove list when we can drop support for CMake 3.22
set(LRSpline_PUBLIC_HEADERS
  ${PROJECT_BINARY_DIR}/LRSpline_version.h
  include/LRSpline/Basisfunction.h
  include/LRSpline/Element.h
  include/LRSpline/HashSet.h
  include/LRSpline/LRSpline.h
  include/LRSpline/LRSplineSurface.h
  include/LRSpline/LRSplineVolume.h
  include/LRSpline/Meshline.h
  include/LRSpline/MeshRectangle.h
  include/LRSpline/Profiler.h
  include/LRSpline/Streamable.h
)

if(CMAKE_VERSION GREATER_EQUAL 3.23)
  target_sources(
    LRSpline
    PUBLIC
    FILE_SET
      HEADERS
    BASE_DIRS
      ${PROJECT_BINARY_DIR}
      include
    FILES
    ${LRSpline_PUBLIC_HEADERS}
  )
endif()

set_target_properties(
  LRSpline
  PROPERTIES
    VERSION
      ${LRSpline_VERSION}
    SOVERSION
      ${LRSpline_VERSION_MAJOR}.${LRSpline_VERSION_MINOR}
)

if(TARGET GoTools::GoTrivariate)
  target_compile_definitions(LRSpline PUBLIC HAS_GOTOOLS=1)
  target_link_libraries(LRSpline PUBLIC GoTools::GoTrivariate)
endif()

set(APP_LIST
  BezierExtract
  Diagonal
  drawLRmesh
  drawLRmultipatch
  MakeWaterfall
  reduceContinuity
  refine
  RefinementUnchanged
  PartitionOfUnityTest
  StresstestEvaluation
  StresstestRefinement
  testIntegral
  TestMatchingKnots
  TestReadWrite
  TopologyRefinement
)

if(TARGET GoTools::GoTrivariate)
  list(APPEND APP_LIST TensorComparison)
endif()

if(Boost_FOUND)
  list(APPEND APP_LIST LinearIndep)
endif()

foreach(APP ${APP_LIST})
  add_executable(${APP} Apps/${APP}.cpp)
  target_link_libraries(${APP} PRIVATE LRSpline)
  set_target_properties(
    ${APP}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY
        bin
  )
endforeach()
if(TARGET LinearIndep)
  target_link_libraries(LinearIndep PRIVATE Boost::boost)
endif()

# Regression tests
enable_testing()

function(add_tests)
  set(oneValueArgs TARGET WORKING_DIRECTORY)
  set(multiValueArgs TEST_FILES)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  foreach(REG ${PARAM_TEST_FILES})
    get_filename_component(TEST_NAME ${REG} NAME)
    get_filename_component(WD ${REG} DIRECTORY)
    add_test(
      NAME
        ${PARAM_TARGET}+${TEST_NAME}
      COMMAND
        ${PROJECT_SOURCE_DIR}/test/regtest.sh
        $<TARGET_FILE:${PARAM_TARGET}>
        ${PROJECT_SOURCE_DIR}/${REG}
        ${PROJECT_BINARY_DIR}
    )
  endforeach()
endfunction()

add_tests(
  TARGET
    RefinementUnchanged
  TEST_FILES
    test/RefinementUnchanged/article_example.reg
    test/RefinementUnchanged/corner_refine.reg
    test/RefinementUnchanged/cross_1_5_m1.reg
    test/RefinementUnchanged/cross_1_5_m2.reg
    test/RefinementUnchanged/discontinuous.reg
    test/RefinementUnchanged/mergeVol.reg
    test/RefinementUnchanged/one_third_raise_mult.reg
    test/RefinementUnchanged/one_third_start_at.reg
    test/RefinementUnchanged/overlap.reg
    test/RefinementUnchanged/vol_3cube_corner.reg
)
add_tests(
  TARGET
    Diagonal
  TEST_FILES
    test/Diagonal/fullspan_p3_m1.reg
    test/Diagonal/minspan_p3_m3.reg
    test/Diagonal/struct_p3_m1.reg
)
add_tests(
  TARGET
    TestReadWrite
  TEST_FILES
    test/TestReadWrite/p3n50.reg
    test/TestReadWrite/p3q6.reg
    test/TestReadWrite/p3.reg
    test/TestReadWrite/vol_p345.reg
)
add_tests(
  TARGET
    TopologyRefinement
  TEST_FILES
    test/TopologyRefinement/orderDep_2D.reg
    test/TopologyRefinement/orderIndep_2D.reg
    test/TopologyRefinement/vol.reg
)
add_tests(
  TARGET
    testIntegral
  TEST_FILES
    test/Integrals/surf_p33.reg
    test/Integrals/vol_p333.reg
)

if(Boost_FOUND)
  add_tests(
    TARGET
      RefinementUnchanged
    TEST_FILES
      test/RefinementUnchanged/LinDep_elongation.reg
      test/RefinementUnchanged/LinDep_newline.reg
  )
  add_tests(
    TARGET
      LinearIndep
    TEST_FILES
      test/LinearIndep/Pettersen_InitPaper_Overload.reg
      test/LinearIndep/Pettersen_InitPaper.reg
  )
endif()

if(TARGET GoTrivariate)
  add_tests(
    TARGET
      TensorComparison
    TEST_FILES
      test/TensorComparison/p4q5.reg
      test/TensorComparison/vol_p345.reg
  )
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/Templates/LRSplineConfig.cmake.in
  ${PROJECT_BINARY_DIR}/LRSplineConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LRSpline
)

export(TARGETS LRSpline FILE LRSplineTargets.cmake NAMESPACE LRSpline::)

if(WIN32)
  set(INSTALL_LIBDIR LRSpline)
  set(INSTALL_INCLUDE_DIR LRSpline/include)
  set(INSTALL_DOCDIR LRSpline/docs)
else()
  set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
  set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR}/LRSpline)
  set(INSTALL_DOCDIR ${CMAKE_INSTALL_DOCDIR}/LRSpline)
endif()

if(CMAKE_VERSION GREATER_EQUAL 3.23)
  install(
    TARGETS
      LRSpline
    EXPORT
      LRSpline
    DESTINATION
      ${INSTALL_LIBDIR}
    FILE_SET
      HEADERS
    DESTINATION
      ${INSTALL_INCLUDEDIR}
  )
else()
  install(
    TARGETS
      LRSpline
    EXPORT
      LRSpline
    DESTINATION
      ${INSTALL_LIBDIR}
  )
  install(
    FILES
      ${LRSpline_PUBLIC_HEADERS}
    DESTINATION
      ${INSTALL_INCLUDEDIR}
  )
endif()

install(
  EXPORT
    LRSpline
  DESTINATION
    ${INSTALL_LIBDIR}/cmake/LRSpline
  FILE
    LRSplineTargets.cmake
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/LRSplineConfig.cmake
  DESTINATION
    ${INSTALL_LIBDIR}/cmake/LRSpline
)

# For generating the doxy
add_custom_target(doc
  doc
  COMMAND
    doxygen ${PROJECT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY
    ${PROJECT_SOURCE_DIR}
  COMMENT
    "Generating API documentation" VERBATIM
)

install(
  CODE
  "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")"
)

install(
  DIRECTORY
    ${PROJECT_BINARY_DIR}/doc/html
  DESTINATION
    ${INSTALL_DOCDIR}/html
  PATTERN
    *.md5 EXCLUDE
  PATTERN
    *.map EXCLUDE
)
