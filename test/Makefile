MKPATH=$(abspath $(lastword $(MAKEFILE_LIST)))
INCDIR=-I./include -I/usr/include -I../../../include
LIBDIR= -L/usr/lib -L/usr/lib/x86-linux-gnu -L/usr/lib/fglrx -L../../../lib
CXX=g++
CXXFLAGS+=$(INCDIR) -Wall -std=c++0x
OBJCFLAGS=$(INCDIR)
LDFLAGS=$(LIBDIR)
INFLAGS?=

LIBS_LINUX=-lfiMath 
LIBS_W32=-lfiMath
LIBS_OSX=-lfiMath

LIBS=$(LIBS_LINUX)

SOURCEDIR=src
TOPDIR=.

### Input Source Files
C_FILES=
CXX_FILES=\
	main.cpp

### Tool / OS specific substitutions
RM=rm -f
CP=cp
SLASH=/
CMDSEP=;

### Targets / Output
TARGET=test
CLEAN_TARGET=$(TARGET)

### Specify project dependencies here
DEPENDENCIES=..$(SLASH)..$(SLASH)fiMath

### Operating system specific defines, etc.
ifeq ($(OS),Windows_NT)
    CXXFLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
        CCFLAGS += -D AMD64
    endif
    ifeq ($(PROCESSOR_ARCHITECTURE),x86)
        CCFLAGS += -D IA32
    endif

    LIBS=$(LIBS_W32)
    CLEAN_TARGET=$(TARGET).exe
	CXX_FILES+=
    LIBDIR +=
	
ifndef BASH_EMU
    RM=del /F /Q
    CP=copy
    SLASH=\\
    CMDSEP=&
endif

else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CXXFLAGS += -D LINUX
		CXX_FILES +=
    endif
    ifeq ($(UNAME_S),Darwin)
	LIBS=$(LIBS_OSX)
	CXXFLAGS += -D OSX
	CXX_FILES +=
	OBJC_FILES +=
	OBJCXX_FILES +=
	LIBDIR:= \
		-L/usr/lib \
		-L../../../lib
    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        CXXFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        CXXFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        CXXFLAGS += -D ARM
    endif
endif

### Set path for output object files
OBJDIR=obj
OBJECTS=$(patsubst %, $(OBJDIR)/%, $(C_FILES:.c=.o))
OBJECTS+=$(patsubst %, $(OBJDIR)/%, $(CXX_FILES:.cpp=.o))
OBJECTS+=$(patsubst %, $(OBJDIR)/%, $(OBJC_FILES:.m=.o))
OBJECTS+=$(patsubst %, $(OBJDIR)/%, $(OBJCXX_FILES:.mm=.o))
		
all: dependencies init $(TARGET) finish

$(OBJDIR)/%.o: $(SOURCEDIR)/%.c
	@echo ""
	@echo " Compile:"				
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: $(SOURCEDIR)/%.cpp
	@echo ""
	@echo " Compile:"				
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: $(SOURCEDIR)/%.m
	@echo ""
	@echo " Compile:"				
	$(CXX) $(OBJCFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: $(SOURCEDIR)/%.mm
	@echo ""
	@echo " Compile:"				
	$(CXX) $(OBJCFLAGS) -c -o $@ $<

$(TARGET): $(OBJECTS)
	@echo ""
	@echo " Link:"
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

dependencies: $(DEPENDENCIES)
	@echo ""
	@echo " Building prerequesite project: $^"
	$(foreach x, $^, $(MAKE) -C $(x)$(CMDSEP))

init:
	@echo ""
	@echo "*************************************"
	@echo "	Building: $(TARGET)"
	@echo "*************************************"
	
finish: 
	@echo ""
	@echo "*Build completed.*"
	@echo ""

.PHONY: clean

deepclean: $(DEPENDENCIES) clean
	$(foreach x, $(DEPENDENCIES), $(MAKE) -C $(x) clean$(CMDSEP)) 

clean:
	$(RM) $(CLEAN_TARGET)
	$(RM) $(OBJDIR)$(SLASH)*.o
	$(RM) core.*
	$(RM) include$(SLASH)*
	$(RM) src$(SLASH)moc_*
