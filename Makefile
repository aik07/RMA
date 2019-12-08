# MAKEFILE

CXX=mpicxx

RMA_DIR=/home/kagawa/Projects/thesis/code/RMA
PEBBL_DIR=/home/kagawa/Projects/thesis/code/pebbl/installpebbl
MPI_ROOT=/usr/lib/x86_64-linux-gnu/openmpi #/opt/openmpi/1.10.1

######################################### SYMBOLS ########################################
SYMBOLS=HAVE_CONFIG_H ANSI_HDRS ANSI_NAMESPACES
DEFSYMBOLS=$(patsubst %, -D%, $(SYMBOLS))

######################################## INCLUDES ##########################################
HEADERDIRS=$(RMA_DIR)/src $(PEBBL_DIR)/include $(MPI_ROOT)/include
INCLUDES=$(patsubst %,-I%,$(HEADERDIRS))

######################################### LIB ############################################
LIBDIRS=$(PEBBL_DIR)/lib $(MPI_ROOT)/lib
LIBLOCATIONS=$(patsubst %,-L%,$(LIBDIRS))
LIBS=pebbl mpi mpi_cxx open-rte open-pal
LIBSPECS=$(patsubst %,-l%,$(LIBS))

########################################## FLAGS ##########################################
DEBUGFLAGS=-g -fpermissive -O0
MISCCXXFLAGS= -std=c++11  #98

# include
CXXFLAGS=$(DEFSYMBOLS) $(INCLUDES) $(MISCCXXFLAGS) $(DEBUGFLAGS)
# Library
LDFLAGS=$(DEBUGFLAGS) $(LIBLOCATIONS)

#####################################################################################

HDRDIR=./src
SRCDIR=./src
OBJDIR=./obj
_HEADERS=Time.h driverRMA.h argRMA.h baseRMA.h dataRMA.h \
         serRMA.h parRMA.h # greedyRMA.h
_SOURCES=driverRMA.cpp argRMA.cpp baseRMA.cpp dataRMA.cpp \
         serRMA.cpp parRMA.cpp  driver.cpp #greedyRMA.cpp
_OBJECTS=$(_SOURCES:.cpp=.o)

SOURCES = $(patsubst %, $(SRCDIR)/%, $(_SOURCES))
HEADERS = $(patsubst %, $(HDRDIR)/%, $(_HEADERS))
OBJECTS = $(patsubst %, $(OBJDIR)/%, $(_OBJECTS))
EXECUTABLE=rma

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LIBSPECS) -o $@

$(OBJDIR)/%.o:  $(SRCDIR)/%.cpp $(HEDADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(EXECUTABLE)
