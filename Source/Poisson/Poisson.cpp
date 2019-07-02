#include "../../Library/FEM.hpp"

#include "../../Library/Poisson.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(5),dim2(5);
    int btype[1] = {1};
    std::string savelocation("./");
    if(argc>1){
        savelocation = argv[1];
        if(argc>3){
            dim1 = std::atoi(argv[2]);
            dim2 = std::atoi(argv[3]);
            if(argc>4){
                CONST_F = std::atof(argv[4]);
            }
        }
    }
    gridfile.open(savelocation+"grid.txt");
    elementfile.open(savelocation+"elements.txt");
    datafile.open(savelocation+"data.txt");
    parameterfile.open(savelocation+"parameters.txt");
    resultfile.open(savelocation+"result.txt");

    DynamicVector solution(dim1*dim2);
    solution.setConstant(0);
    DynPlaneTriangleFem<1> Model(dim1,dim2,0,1,0,1);
    Model.Build_Vector(Vector_Map);
    Model.Build_Matrix(Matrix_Map);
    Model.Set_Boundary(solution,btype);
    solution += Model.SparseLU();

    for(int i=0;i<dim1*dim2;i++){
        datafile<<Model.grid.nodes[i].eval()+" "+std::to_string(solution(i))+"\n";
    }
    gridfile<<Model.grid.evalnodes();
    elementfile<<Model.grid.evalelements();
    parameterfile<<dim1<<" "<<dim2<<" "<<CONST_F<<"\n";
    resultfile<<solution.maxCoeff()<<"\n";
    gridfile.close();
    elementfile.close();
    datafile.close();
    parameterfile.close();
    resultfile.close();
    return 0;
}