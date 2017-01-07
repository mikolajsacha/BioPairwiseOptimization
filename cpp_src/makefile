CXX=g++
CXXFLAGS=-std=c++14 -W -Wall -Werror -pedantic -Wextra

ODIR=obj
LDIR=lib
IDIR=include
TESTSDIR=tests

LIBS=-lm

_OBJS = main.o match_scorers.o global_alignment.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS = match_scorers.h global_alignment.h model.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

$(ODIR)/%.o: $(LDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

makeall: $(OBJS)
	$(CXX) -o runtests $(OBJS) $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
