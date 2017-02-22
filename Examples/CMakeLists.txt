PROJECT(Hello_LRSpline)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Required packages
FIND_PACKAGE(LRSpline REQUIRED)

# Use same definitions as was compiled with LRspline 
LIST(APPEND DEFINITIONS ${LRSpline_DEFINITIONS})

# specify executable
ADD_EXECUTABLE(hello_lrspline Hello_LRSpline.cpp)

# compiling with c++11 allows us to use "for(auto b : some_container)" loop constructs
TARGET_COMPILE_FEATURES(hello_lrspline PRIVATE cxx_range_for)

# link with compiled LRSpline library
TARGET_LINK_LIBRARIES(hello_lrspline LRSpline)

