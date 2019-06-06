This directory must be enriched with two additional files besides this README.

***** On Windows MSVC *****

The two additional files are:
1) mpir.h
2) One of: mpir64.lib (on MS Windows 64bit), mpir32.lib (on MS Windows 32bit)

To create these two files:

Option (1): check if your MSVC version matches one of the subfolders in this diorectory. If so, just copy from there. Otherwise use option 2 below.

Option (2):

1) Clone or download MPIR from https://github.com/BrianGladman/mpir
2) cd mpir-master/msvc/vsXX (XX is one of 13, 15, 17, 19, ...select the one matching your MSVC version)
3) open mpir.sln
4) select Debug/Release and Win32/x64 based on your needs
5) build the project lib_mpir_gc
6) You will find both mpir.h and mpir.lib in msvc/vsXX/lib_mpir_gc/x64/Release/ (if you compiled for x64 in Release mode).
7) Copy these two files in the directory containing this README file.
8) Rename mpir.lib to mpir64.lib (if you compiled for x64, otherwise rename to mpir32.lib).


***** On Mac OSX *****

The two additional files are:
1) mpir.h
2) libmpir.a

To create these two files:
1) brew install mpir 
2) The two files will be placed in /usr/local/lib/libmpir.a and /usr/local/include/mpir.h
3) Copy both of them in the directory containing this README file.


***** On Linux *****

The two additional files are:
1) mpir.h
2) libmpir.a

To create these two files:
1) sudo apt-get install yasm m4 build-essential unzip wget -y
2) wget http://mpir.org/mpir-3.0.0
3) unzip mpir-3.0.0.zip
4) cd mpir-3.0.0
5) ./configure
6) make
7) The two files will be placed in .libs/libmpir.a and ./mpir.h
8) Copy these two files in the directory containing this README file.
