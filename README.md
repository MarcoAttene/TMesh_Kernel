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
The file tests/test_kernel.cpp is a good starting point to start playing with the library.
It has been extensively commented to step-by-step explain how to use the main library features.
Other than that, please look at the comments within the header files or use doxygen to
produce documentation in a more readable format.


-------------------
System Requirements
--------------------

TMesh_Kernel has been tested on 64 bit PCs running:
 - Microsoft Windows OS with MSVC
 - Linux with standard gcc/g++ development environment
 - Mac OSX

TMesh_Kernel exploits MPIR to deal with arbitrarily precise rational numbers.
MPIR is an independent software licensed under the terms of GNU LGPL v3.
Please look at mpir/README.txt for further information.

-------------------------------------
Building the tree in LAZY-HYBRID mode
-------------------------------------

Before building TMesh_Kernel, you must obtain and build MPIR.
Please look at mpir/README.txt for instructions.

Once done with MPIR, return to the directory containing this README file,
edit CMakeLists.txt at line 8 to enter the path containing mpir.a (or mpir.lib on Windows).
Then type:
```
mkdir build
cd build
cmake ..
```

This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a TMesh_Kernel.sln file.
On Linux/OSx, this will produce a Makefile. 
Use it as usual to compile TMesh_Kernel.

------------------------------
Building the tree in FAST mode
------------------------------

In this mode, the hybrid kernel is disabled.
MPIR is not needed to compile in this mode.
Just run the following commands:
```
mkdir build
cd build
cmake -DUSE_LAZY_KERNEL=OFF ..
```

and proceed as for the LAZY-HYBRID mode.
**Warning:** When compiled in FAST mode, rational numbers are not supported!
If you need the hybrid arithmetic kernel, compile in LAZY-HYBRID mode.

-------------------------------------
Using the library in LAZY-HYBRID mode
-------------------------------------

C/C++ code using TMesh_Kernel must be compiled according to the following rules.

**include path must include:**
```
$(TMESH_HOME)/include
$(TMESH_HOME)/mpir
```

**preprocessor definitions**
In all cases: `USE_HYBRID_KERNEL`, `USE_LAZY_KERNEL`, `IS64BITPLATFORM`

**library path must include:**
```
$(TMESH_HOME)/lib64
$(TMESH_HOME)/mpir
```

**static libraries to be linked**
`kernel_Lazy64`
`mpir`

**other options to be set**
Support for OpenMP might be active

------------------------------
Using the library in FAST mode
------------------------------
As above, but all the references to mpir are no longer necessary.
Furthermore, the only preprocessor definition to be used is:
`IS64BITPLATFORM`

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
