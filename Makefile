CXX=mpicxx

# You have to specify the directory location where you placed acro-pebbl
PEBBL_ROOT=/home/kagawa/Projects/thesis/installpebbl

# You have to specify the location where you placed openmpi
# Example:
MPI_ROOT=/usr/lib/x86_64-linux-gnu/openmpi

######################################## SYMBOLS ########################################
SYMBOLS=HAVE_CONFIG_H ANSI_HDRS ANSI_NAMESPACES
DEFSYMBOLS=$(patsubst %,-D%,$(SYMBOLS))

####################################### INCLUDES ##########################################
HEADERDIRS=./src $(PEBBL_ROOT)/include $(MPI_ROOT)/include 
INCLUDES=$(patsubst %,-I%,$(HEADERDIRS))

######################################### LIB ############################################
LIBDIRS=$(PEBBL_ROOT)/lib  $(MPI_ROOT)/lib
LIBLOCATIONS=$(patsubst %,-L%,$(LIBDIRS))
LIBS=pebbl mpi mpi_cxx open-rte open-pal # utilib
LIBSPECS=$(patsubst %,-l%,$(LIBS)) 

####################################### FLAGS ##########################################
DEBUGFLAGS=-g -fpermissive -O0
MISCCXXFLAGS= -std=c++11 #-std=c++98

# include
CXXFLAGS=$(DEFSYMBOLS) $(INCLUDES) $(MISCCXXFLAGS) $(DEBUGFLAGS) 

# Library
LDFLAGS=$(DEBUGFLAGS) $(LIBLOCATIONS)   

######################################################################################

SRCDIR=./src
OBJDIR=./obj
_HEADERS=Time.h argRMA.h baseRMA.h dataRMA.h serRMA.h parRMA.h
_SOURCES=main.cpp argRMA.cpp baseRMA.cpp serRMA.cpp parRMA.cpp
_OBJECTS = $(_SOURCES:%.cpp=%.o)
SOURCES = $(patsubst %,$(SRCDIR)/%,$(_SOURCES)) 
HEADERS = $(patsubst %,$(SRCDIR)/%,$(_HEADERS))
OBJECTS = $(patsubst %,$(OBJDIR)/%,$(_OBJECTS))
EXECUTABLE=rma

#####################################################################################

all: $(EXECUTABLE) $(SOURCES)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $^ $(LIBSPECS) -o $@

$(OBJDIR)/%.o:  $(SRCDIR)/%.cpp $(HEDADERS)  
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(OBJDIR)/*.o $(EXECUTABLE) $(SRCDIR)/*~ *~
