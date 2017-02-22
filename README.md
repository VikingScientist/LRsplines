# LR B-splines

### Locally Refined B-splines

LR B-splines is a technology which enables users to perform local refinement on B-spline surfaces or volumes, which traditionally have been limited to tensor products. This has a wide variety of applications within computer-aided design (CAD) and computer-aided engineering (CAE). While this library was written with the latter in mind, it is also possible to take use of it in a design environment. For more information consult the website : http://lrbsplines.com

### Getting the code

Simply press the green button on the top right of the github screen that states "Clone or download". On ubuntu (or with git bash on any OS) you can clone the library by typing in 

  `git clone https://github.com/VikingScientist/LRsplines`

### Compiling the code

This is a c++ library and you will have to compile and link this to your applications accordingly. It uses CMake as the compilation tools. There are two optional dependencies which is boost and GoTools. These will expand the capabilities of the library, but are not required for compilation.Boost allows for more linear dependency testing and GoTools is a tensor product B-spline library.

#### Ubuntu 

1. Install all compilation tools that we are going to need

  `sudo apt-get install cmake g++`

2. **[Optional]**: Install boost

  `sudo apt-get install libboost-dev`

3. **[Optional]**: Install GoTools. This is compiled on launchpad. Add https://launchpad.net/~ifem/ (follow the instructions on site) followed by typing

    `sudo apt-get install libgotools-core-dev libgotools-trivariate-dev`

To compile the code, first navigate to the root catalogue of LR-splines, here denoted by `<LRSpline root>` 

1. `cd <LRSpline root>`
2. `mkdir Release`
3. `cd Release`
4. `cmake -DCMAKE_BUILD_TYPE=Release ..`
5. `make`
6. **[Optional]**: `make test`
7. **[Optional]**: `sudo make install`

Point 6 will run a series of test to verify that the library compiled correctly and is running as it should. Point 7 will install this on your system by placing the header-files, library-files and cmake-files in their right place (default: `/usr/local/lib` and `/usr/local/include`): this will in turn make it much easier to compile your own applications which uses this library. Changing any instance of `Release` with `Debug` makes the library compile with debug-flags on.

#### Windows

1. Download and install cmake: https://cmake.org/download/
2. Download and install visual studio: https://www.visualstudio.com/downloads/
3. Compile a visual studio project by using cmake. See https://cmake.org/runningcmake/ for instructions
4. Open the LRspline project in visual studio and build the project.

CMake will generate a big red warning that it cannot find GoTools. This is as expected and it will compile fine without it.

### Using the Code

See the `Examples/` folder for a sample c++ program that uses the LRSpline library. It assumes that CMake is able to find the config files `LRSplineConfig.cmake` as well as all header and library files. If these are not installed on default location (by `sudo make install`), then you will have to specify them manually. `Hello_LRSpline.cpp` is compiled and run by

1. `cd Examples`
2. `cmake .`
3. `make`
4. `./hello_lrspline`
