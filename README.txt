----------------------------
TMesh_Kernel
----------------------------

by Marco Attene

Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            

TMesh_Kernel is a C++ library providing a hybrid arithmetic kernel for geometric computation in 3D.
It includes hybrid numbers, 3D vectors and points along with standard operations, geometric predicates, geometric constructors, and intersections.
The hybrid kernel combines the efficiency of standard floating point computation with the robustness of rational arithmetic.
In TMesh_Kernel, geometric predicates are always exact, whereas geometric constructions can be either approximate or exact, depending on the programmer's choice. Programs based on TMesh_Kernel can dynamically switch between rational and floating point arithmetic.
TMesh_Kernel is at the base of a larger project called ImatiSTL.

-------------------
Citation policy
--------------------
You are free to use TMesh_Kernel according to the licensing terms specified at the end of this document.
If you use TMesh_Kernel for research purposes and produces publications, please cite the following paper 
that describes the underlying hybrid paradigm:

> Marco Attene. ImatiSTL - Fast and Reliable Mesh Processing with a Hybrid Kernel.
  Springer LNCS Transactions on Computational Science, Vol. XXIX, pp. 86-96, 2017.

Note that the distinction between 'approximated' and 'filtered' described in the paper has been deprecated
in the code, where the only two possibilities are 'filtered' and 'precise'. This simplification is because
the performance loss of a filtered kernel wrt an approximated kernel is mostly negligible.


-------------------
Documentation
-------------------

The aforementioned paper is the best source of information to understand the hybrid arithmetic paradigm.
However, the paper is not meant to be a comprehensive programmer's guide.
The file src/test_kernel.cpp is a good starting point to start playing with the library.
It has been extensively commented to step-by-step explain how to use the main library features.
Other than that, please look at the comments within the header files or use doxygen to
produce documentation in a more readable format.


-------------------
System Requirements
--------------------

TMesh_Kernel has been tested on 32 and 64 bit PCs running:
 - Microsoft Windows OS with MSVC 12.0 (Visual C++ 2013)
 - Linux with standard gcc/g++ development environment

TMesh_Kernel exploits MPIR to implement exact arithmetic.
MPIR is used to deal with multiprecision and rational coordinates.
Precompiled static .lib files are provided for Windows systems
in the 'mpir' folder. A precompiled static .a file is also provided
for 64-bit Linux systems. If you are using another system you need
to download MPIR at http://www.mpir.org/ and compile the library yourself.
Note that MPIR ***IS NOT*** part of TMesh_Kernel. The pre-compiled libs and
header files are provided along with this distribution in the hope to
ease the use of MPIR within TMesh_Kernel.
MPIR is an independent software licensed under the terms of GNU LGPL v3.
Both the header files and the pre-compiled .lib and .a files included in
'mpir' are produced out of MPIR version 3.0, released on 2017-03-01.
Please look at mpir/README.txt for further information.

-------------------
Building the tree
-------------------

On MS Windows:
TMesh_Kernel can be compiled for both 32 and 64 bit architectures.
It can be compiled in three different modes:
'Fast' = top speed / lowest memory footprint / non-robust (self-contained / no dependencies)
'Hybrid' = slowest / intermediate memory footprint / robust (uses MPIR)
'Lazy' = intermediate speed / high memory footprint / robust (uses MPIR)

In 'Fast' mode, the basic 'coord' type is simply an IEEE standard double.
This mode has been included mainly for benchmarking purposes as it is not hybrid.
In 'Hybrid' mode, the 'coord' type is a polymorphic number that can be either
a standard double or an exact rational number (see the paper for details).
'Lazy' mode is equivalent to 'Hybrid', but the exact rational number is evaluated
only when a predicate requires it.

Double-click on build/msvc/TMesh_Kernel.sln.
Select the proper combination of Configuration/Platform based on your needs.
Configuration can be: Fast, Hybrid, Lazy
Platform can be: Win32, x64
Press f7 to compile the solution.

On Linux/Unix:
> cd build/linux-gcc
> make
> cd ../..
> mv bin64/*.a lib64

This should produce both the library (in lib64/) and a test executable (in bin64/).
Note that the Makefile provided for Linux is only able to compile 64 bit versions.

Compilation options through preprocessor definitions:
1) IS64BITPLATFORM - set to compile a 64bit-compliant library
2) USE_HYBRID_KERNEL - set to enable the hybrid kernel. If set, you need to have Mpir installed.
3) USE_LAZY_KERNEL - set to enable lazy evaluation. If set, USE_HYBRID_KERNEL must be set too.
4) EXTENSIBLE_TMESH - set to enable polymorphism. Usually not needed.


-------------------
Using the library
-------------------

The easiest way to start is to copy the compilation rules from kernel_test.vcxproj (on 
Windows) or from the Makefile provided for Linux-gcc.
In any case, here are general rule to compile your application.

--- include path must include:
$(TMESH_HOME)/include
$(TMESH_HOME)/mpir (only if USE_HYBRID_KERNEL)

--- preprocessor definitions
must be the same used to compile the library (see above)

--- library path must include:
$(TMESH_HOME)/lib (if compiled for 32bit)
$(TMESH_HOME)/lib64 (if compiled for 64bit)
$(TMESH_HOME)/mpir

--- static libraries to be linked
kernel_XXX.lib (where XXX can be either Fast, Hybrid or Lazy, possibly with 64 extension)
mpirXXX.lib (only if USE_HYBRID_KERNEL, XXX is either 32 or 64)
On Linux systems the YYY.lib is replaced with libYYY.a

--- other options to be set
Support for OpenMP might be active

---------------------
Copyright and license
---------------------

TMesh_Kernel
Authors: Marco Attene                                                    

Copyright(C) 2012: IMATI-GE / CNR                                        

IMATI-GE / CNR is Consiglio Nazionale delle Ricerche                     
Istituto di Matematica Applicata e Tecnologie Informatiche               
Genova (Italy)                                                           

TMesh_Kernel is free software; you can redistribute it and/or modify     
it under the terms of the GNU Lesser General Public License as published 
by the Free Software Foundation; either version 3 of the License, or (at 
your option) any later version.                                          

TMesh_Kernel is distributed in the hope that it will be useful, but      
WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser 
General Public License for more details.                                 

You should have received a copy of the GNU Lesser General Public License 
along with the TMesh_Kernel.  If not, see http://www.gnu.org/licenses/.
