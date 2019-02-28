#############################################################################
# Makefile for building: libTMesh_Kernel.a and libTMesh_Kernel_test
#############################################################################

MAKEFILE      = Makefile

####### Dependencies

MPIR_DIR	= ../../mpir
MPIR_LIBS 	= $(MPIR_DIR)

TMESH_KERNEL_SRC	= ../../src
TMESH_KERNEL_INCLUDE	= ../../include
TMESH_KERNEL_BUILD	= ../../bin64

BUILD_OPTIONS = -DUSE_HYBRID_KERNEL -DUSE_LAZY_KERNEL

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = $(BUILD_OPTIONS)
CFLAGS        = -pipe -g -fPIC -Wall -W -D_REENTRANT $(DEFINES) -fopenmp
CXXFLAGS      = -pipe -g -fPIC -Wall -W -D_REENTRANT $(DEFINES) -fopenmp
INCPATH       = -I$(MPIR_DIR) -I$(TMESH_KERNEL_INCLUDE)
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
DISTNAME      = TMesh_Kernel
AR            = ar rcs


####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES_TEST	= $(TMESH_KERNEL_SRC)/test_kernel.cpp 
OBJECTS_TEST	= $(TMESH_KERNEL_BUILD)/test_kernel.o 

SOURCES       = $(TMESH_KERNEL_SRC)/mixedPredicates.cpp \
		$(TMESH_KERNEL_SRC)/tmesh.cpp \
		$(TMESH_KERNEL_SRC)/point.cpp \
		$(TMESH_KERNEL_SRC)/orientation.cpp \
		$(TMESH_KERNEL_SRC)/expansion.cpp \
		$(TMESH_KERNEL_SRC)/coordinates.cpp 

OBJECTS       = $(TMESH_KERNEL_BUILD)/mixedPredicates.o \
		$(TMESH_KERNEL_BUILD)/tmesh.o \
		$(TMESH_KERNEL_BUILD)/point.o \
		$(TMESH_KERNEL_BUILD)/orientation.o \
		$(TMESH_KERNEL_BUILD)/expansion.o \
		$(TMESH_KERNEL_BUILD)/coordinates.o

DESTDIR       = build
TARGET        = all

first: all
####### Build rules

all: libTMesh_Kernel.a libTMesh_Kernel_test

libTMesh_Kernel_test: $(OBJECTS_TEST) $(TMESH_KERNEL_BUILD)/libTMesh_Kernel.a
	-$(DEL_FILE) $(TARGET)
	$(MKDIR) $(TMESH_KERNEL_BUILD)
	$(CXX) -o $(TMESH_KERNEL_BUILD)/libTMesh_Kernel_test $(OBJECTS_TEST) $(TMESH_KERNEL_BUILD)/libTMesh_Kernel.a $(MPIR_LIBS)/libmpir.a -fopenmp

libTMesh_Kernel.a:  $(OBJECTS) 
	-$(DEL_FILE) $(TARGET)
	$(MKDIR) $(TMESH_KERNEL_BUILD)/
	$(AR) $(TMESH_KERNEL_BUILD)/libTMesh_Kernel.a $(OBJECTS) 
clean: 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) $(OBJECTS_TEST)
	-$(DEL_FILE) *~ core *.core


####### Compile

$(TMESH_KERNEL_BUILD)/%.o : $(TMESH_KERNEL_SRC)/%.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

