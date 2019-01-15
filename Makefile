CXX=mpicxx

# You have to specify the directory location where you placed acro-pebbl
ACRO_ROOT=<DIRECTORY>/acro-pebbl 

# You have to specify the location where you placed openmpi
# Example:
MPI_ROOT=/opt/openmpi/1.10.1

######################################## SYMBOLS ########################################
SYMBOLS=HAVE_CONFIG_H ANSI_HDRS ANSI_NAMESPACES
DEFSYMBOLS=$(patsubst %,-D%,$(SYMBOLS))

####################################### INCLUDES ##########################################
HEADERDIRS=. $(ACRO_ROOT)/include $(ACRO_ROOT)/include/pebbl $(ACRO_ROOT)/include/utilib $(MPI_ROOT)/include 
INCLUDES=$(patsubst %,-I%,$(HEADERDIRS))

######################################### LIB ############################################
LIBDIRS=$(ACRO_ROOT)/lib  $(MPI_ROOT)/lib
LIBLOCATIONS=$(patsubst %,-L%,$(LIBDIRS))
LIBS=pebbl utilib mpi mpi_cxx open-rte open-pal
LIBSPECS=$(patsubst %,-l%,$(LIBS)) 

####################################### FLAGS ##########################################
DEBUGFLAGS=-g -fpermissive -O0
MISCCXXFLAGS= -std=c++98

# include
CXXFLAGS=$(DEFSYMBOLS) $(INCLUDES) $(MISCCXXFLAGS) $(DEBUGFLAGS) 

# Library
LDFLAGS=$(DEBUGFLAGS) $(LIBLOCATIONS)   

######################################################################################

SRCDIR=./src
OBJDIR=./
_HEADERS=serRMA.h parRMA.h
_SOURCES=main.cpp serRMA.cpp parRMA.cpp
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
