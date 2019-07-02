PROGRAM= FEM_U2 Membrane MembraneParameters Membrane4 Membrane5 MembraneParameters5 U2_Parameters

all : $(PROGRAM)

% : %.cpp
	c++ -Wall -std=c++14 -o $@ $<

clean :
	rm -f *~ *.o $(PROGRAM)
