XX=g++
CXXFLAGS=-std=c++14 -W -Wall -Werror -pedantic -Wextra

ODIR=obj
LDIR=lib
SODIR=../lib
IDIR=include
TDIR=../test
LIBNAME=pairwise.so

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

LIBS=-lm

_OBJS = main.o match_scorers.o global_alignment.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS = match_scorers.h global_alignment.h model.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

$(ODIR)/%.o: $(LDIR)/%.cpp $(TESTSDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fPIC

all: makelib copyfiles

makelib: $(OBJS) 
	$(CXX) $(LDIR)/integration.cpp -lboost_python -lpython2.7 -I /usr/include/python2.7/ -o $(SODIR)/$(LIBNAME) $(OBJS) -shared $(CXXFLAGS) $(LIBS) -fPIC -static-libstdc++

# copy compiled library to test folder for testing
copyfiles:
	cp $(SODIR)/$(LIBNAME) $(TDIR)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.so
