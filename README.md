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
If you use TMesh_Kernel for research purposes and produce publications, please cite the following paper 
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
 - Microsoft Windows OS with MSVC
 - Linux with standard gcc/g++ development environment
 - Mac OSX

TMesh_Kernel exploits MPIR to deal with arbitrarily precise rational numbers.
MPIR is an independent software licensed under the terms of GNU LGPL v3.
Please look at mpir/README.txt for further information.

-------------------
Building the tree
-------------------

Before building TMesh_Kernel, you must obtain and build MPIR.
Please look at mpir/README.txt for instructions.

Once done with MPIR, return to the directory containing this README file and type:
```
mkdir build
cd build
cmake ..
```

This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a TMesh_Kernel.sln file
On Linux/OSx, this will produce a Makefile

Thanks to Teseo Schneider for having produced the CMake environment for TMesh_Kernel.

-------------------
Using the library
-------------------

C/C++ code using TMesh_Kernel must be compiled according to the following rules.

**include path must include:**
```
$(TMESH_HOME)/include
$(TMESH_HOME)/mpir
```

**preprocessor definitions**

In all cases: `USE_HYBRID_KERNEL` and `USE_LAZY_KERNEL`
Only for 64bit builds: IS64BITPLATFORM

**library path must include:**
```
$(TMESH_HOME)/lib (if compiled for 32bit)
$(TMESH_HOME)/lib64 (if compiled for 64bit)
$(TMESH_HOME)/mpir
```

**static libraries to be linked**

`kernel_Lazy.lib` (or `kernel_Lazy64.lib` for 64bit builds)
`mpirXXX.lib` (XXX is either 32 or 64)
On Linux/OSX the YYY.lib is replaced with `libYYY.a`

**other options to be set**
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
