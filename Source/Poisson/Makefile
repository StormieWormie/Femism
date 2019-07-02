PROGRAM=PoissonParameters Poisson Poisson2

all : $(PROGRAM)

% : %.cpp
	c++ -Wall -std=c++14 -o $@ $<

clean :
	rm -f *~ *.o $(PROGRAM)
