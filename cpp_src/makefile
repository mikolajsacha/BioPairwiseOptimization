CXX=g++
CXXFLAGS=-std=c++14 -W -Wall -Werror -pedantic -Wextra

ODIR=obj
LDIR=lib
IDIR=include
TESTSDIR=tests
INTDIR=integration

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

LIBS=-lm

_OBJS = main.o match_scorers.o global_alignment.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS = match_scorers.h global_alignment.h model.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

$(ODIR)/%.o: $(LDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

makeall: maketests makedll
	
maketests: $(OBJS)
	$(CXX) -o runtests $(OBJS) $(CXXFLAGS) $(LIBS) -v

makedll: $(OBJS) 
	$(CXX) $(INTDIR)/integration.cpp -lboost_python -lpython2.7 -I /usr/include/python2.7/ -o pairwise.so -shared $(CXXFLAGS) $(LIBS) -fPIC -static-libstdc++ -v
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ *.so
