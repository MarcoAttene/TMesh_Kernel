This directory contains header files with both C and C++ bindings to MPIR-3.0.0.

It also contains pre-compiled static libraries to be linked to your applications.

On MS Windows systems
mpir32.lib is optimized for processors of the Pentium 4 SSE2 family.
mpir64.lib is optimized for processors of the Intel(TM) core-2(TM) family.
mpir32-generic.lib is for generic x86-based architectures.
mpir64-generic.lib is for generic x64-based architectures.

By default, ImatiSTL links to mpir32.lib and mpir64.lib library files.
If you want to create software to be distributed and you want to maximize compatibility, then use the generic versions instead.
To use the generic versions, you need to rename these .lib files accordingly (*-generic.lib -> *.lib).


On Linux systems:
libmpir.a was compiled on a generic 64 bit Ubuntu Linux


All the static libraries included in this distribution were compiled starting from the original MPIR-3.0.0 source code.
This code is freely available for download at:
http://www.mpir.org/

You may replace these pre-compiled versions with more recent and/or customized MPIR versions by simply compiling MPIR
yourself. You just need to replace the .lib (or the .a) files with those produced by your compiler. You also need to replace the
header files with those you used to compile MPIR itself.

Note that the two headers mpir.h and mpirxx.h included in this folder are just redirections to the actual MPIR headers.
