#include "../../Library/FEM.hpp"

#include "../../Library/MembraneFunctions.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile;
    gridfile.open("grid.txt");
    elementfile.open("elements.txt");
    print("test1");
    DynTriangulatedDoughnut Model(4,4,0,1);
    print("test2");
    gridfile<<Model.evalnodes();
    elementfile<<Model.evalelements();
    gridfile.close();
    elementfile.close();
    return 0;
}