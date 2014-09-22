
IDIR =./include
CC=g++
CFLAGS=-I$(IDIR) -I/usr/include/gdal -lcurl -lboost_program_options -lboost_iostreams

ADIR=./apps
SDIR=./src

ODIR=./obj
LDIR =./lib

LIBS=-lm -lgdal

_APP_SRC = pp2g.cpp
APP_SRC = $(patsubst %,$(ADIR)/%,$(_APP_SRC))
APP_OBJ=  $(patsubst %.cpp,$(ODIR)/%.o,$(_APP_SRC))


_LIB_SRC = GridFile.cpp GridMap.cpp InCoreInterp.cpp Interpolation.cpp OutCoreInterp.cpp
LIB_SRC = $(patsubst %,$(SDIR)/%,$(_LIB_SRC))
LIB_OBJ=  $(patsubst %.cpp,$(ODIR)/%.o,$(_LIB_SRC))

SRC = $(APP_SRC) $(LIB_SRC)

#depend: .depend
#.depend: $(SRC)
#	rm -f ./.depend
#	$(CC) $(CFLAGS) -MM $^ -MF  ./.depend;
#include .depend


pp2g: $(APP_OBJ) $(LIB_OBJ)
	g++ -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: $(ADIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o ./.depend pp2p