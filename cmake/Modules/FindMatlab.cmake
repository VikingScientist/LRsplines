find_program(MATLAB_MEX NAMES mex HINTS /usr/lib/matlab2012a/bin)

find_program(MATLAB_CC NAMES gcc-4.4 gcc)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Matlab DEFAULT_MSG MATLAB_MEX MATLAB_CC)
