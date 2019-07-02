PROGRAM=demo

all : $(PROGRAM)

$(PROGRAM) : $(PROGRAM).cpp
	c++ -Wall -std=c++14 -o $(PROGRAM) $(PROGRAM).cpp

clean :
	rm -f *~ *.o $(PROGRAM)
