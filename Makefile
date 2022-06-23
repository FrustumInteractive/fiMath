MKPATH=$(abspath $(lastword $(MAKEFILE_LIST)))

SOURCEDIR=src
TOPDIR=../..
INCDIR=-I$(TOPDIR)/include
LIBDIR=-L$(TOPDIR)/lib

### Input Source Files
C_FILES=
CXX_FILES=\
	frustum.cpp \
	matrix.cpp \
	plane.cpp

OBJC_FILES=
OBJCXX_FILES=
#C_FILES_LINUX
#CXX_FILES_LINUX
#C_FILES_WIN
#CXX_FILES_WIN
#C_FILES_OSX
#CXX_FILES_OSX
#OBJC_FILES_OSX
#OBJCXX_FILES_OSX

### Targets / Output / Resource folders and files
TARGET=libfiMath.a
CLEAN_TARGET=$(TARGET)
RESOURCES=

### Specify project dependencies here
DEP_WIN32=
DEP_OSX=
DEP_LINUX=
DEPENDENCIES=

### Post Build Steps
POSTBUILD=
POSTBUILDCLEAN=

include $(TOPDIR)/templates/Makefile_include_static_lib
