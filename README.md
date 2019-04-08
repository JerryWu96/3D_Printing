# Computational Fabrication (3D Printing) - Course Project	

I am passionate about 3D Printing, so i decided to give it a try!

This repo contains all the <b>ongoing</b> projects I participated in CS581 Computational Fabrication during 2019 Spring at BU.

By May I am going to finish a complete 3D Printing project, which includes the design of an object using CAS software, and print it out using 3D printers (for example, printers from FormLabs!)

* Currently I finished a warm-up coding assignment. It is a voxelizer which you can use to convert a 3D model(.obj) into its voxelized form. Samples and outputs are also provided!



# Usage

## Install CMake
I used CLion to build the project. Simply import the header files and source files and run it directly.

You can also use terminal/console for the smae purpose with the help of CMake!

To install CMake easier, personally i recommend to use Cmake with GUI if you are new to this field.

If you take it seriously, please use CMake with command lines to get a better understanding of it.

* Hint: Use homebrew to easily install and manage CMake if you prefer MacOS.

* Refer to `https://cmake.org/install/` for detailed installation guide.

## Build with CMake

Attention! `./voxelizer/CMakeLists.txt` is used for building in IDEs (As I mentioned, I used Clion). 

If you want to build using terminal/console, please use `./voxelizer/terminal/CMakeLists.txt`

* MacOS-specific guide:
  * Install CMake by using brew: brew install cmake
  
  * Go to the root directory of the assignment, where you can find a text file called CMakeList.txt.
  
  * Create a folder called "build" with command: `mkdir ./build`
  
  * Then go into "build" folder, run command:  `cmake ../`
  
  * Determine Your Current Directory with command `pwd` if necessary!
  
  * If everything is okay, you can start to make executable from source use command `make` at current directory.

  * Every time after editing any source file with texteditoers(for example, Sublime Text), you need to run `make` command again to build again.. 

* Windows-specific guide:

  * Install MINGW from https://sourceforge.net/projects/mingw/files/latest/download?source=files
  
  * Add  C:\MinGW\bin to path system variable
  
  * Install CMake from https://cmake.org/files/v3.10/cmake-3.10.2-win64-x64.msi
  
  * Put CMakeLists inside source folder
  
  * Change the Line inside CMakeLists to:
  
  CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
  
  project(voxelizer)
  
  file(GLOB_RECURSE HEADER_CODE *.h)
  
  file(GLOB_RECURSE SRC_CODE *.cpp)
  
  ADD_EXECUTABLE(voxelizer ${SRC_CODE} ${HEADER_CODE})
  
  * Create the directory Voxelizer/build
  
  * Run the CMake GUI, give the source path and the build path
  
  * Click: 'configure' and specify as generator the MinGW Makefiles. Press: 'finish' then 'generate'
  
  * Change to the /build directory
  
  * Run: `make OR mingw32-make`
  
  * Troubleshooting: may need to substitute exit(0) with return 0. May need to remove #define NOMINMAX 

## Execute the sample program
To run the sample program on MacOS/linux:

1. Open your terminal/console, go to the directory where the sample executable stays:
  `cd: <WORKINGDIR>`
For instance: `cd /Users/jerrywux/Desktop/`
2. run `./voxelizer <InputMeshFileDir OutputMeshFilenDir>`
For instance: `./voxelizer /Users/jerrywux/Desktop/sphere.obj /Users/jerrywux/Desktop/sphere64.obj`

You are all set!

# Samples (before and after)
For Mac OS, you can use meshlab to view .obj files.

sphere.obj:
<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/sphere%20original.png"></img>

voxels with 32x32x32 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/SPHERE%2032x32x32.png"></img>

voxels with 64x64x64 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/SPHERE%2064x64x64.png"></img>

teapot.obj:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/teapot%20original.png"></img>

voxels with 32x32x32 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/TEAPOT%2032x32x32.png"></img>

voxels with 64x64x64 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/TEAPOT%2064x64x64.png"></img>

bunny.obj:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/bunny%20original.png"></img>

voxels with 32x32x32 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/BUNNY%2032x32x32.png"></img>

voxels with 64x64x64 diameter:

<img src="https://github.com/JerryWu96/3D_Printing/blob/master/Voxelizer/sample%20pics/BUNNY%2064x64x64.png"></img>

## Credits

Boston University Department of Computer Science, CS581, Professor Whiting

This repo was created by:

Xiankang (Jerry) Wu