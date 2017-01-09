XX=g++
CXXFLAGS=-std=c++14 -W -Wall -Werror -pedantic -Wextra

ODIR=obj
LDIR=lib
SODIR=../lib
IDIR=include
TDIR=../test
LIBNAME=pairwise

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

LIBS=-lm

_OBJS = integration.o match_scorers.o global_alignment.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS = match_scorers.h global_alignment.h model.h debug.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

all: makelib copyfiles

makelib: $(OBJS) 
	@mkdir -p $(SODIR)
	$(CXX) -lboost_python -lpython2.7 -I /usr/include/python2.7/ -o $(SODIR)/$(LIBNAME).so $(OBJS) -shared $(CXXFLAGS) $(LIBS) -fPIC -static-libstdc++

$(ODIR)/%.o: $(LDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fPIC

# copy compiled library to test folder for testing
copyfiles:
	cp $(SODIR)/$(LIBNAME).so $(TDIR)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(SODIR)*.so 
