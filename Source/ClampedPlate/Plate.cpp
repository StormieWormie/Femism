#include "../../Library/FEM.hpp"

#include "../../Library/ClampedPlate.hpp"



int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(5),dim2(5);
    int btype[3] = {0,0,0};
    std::string savelocation("./");

    if(argc>1){
        savelocation = argv[1];
        if(argc>3){
            dim1 = std::atoi(argv[2]);
            dim2 = std::atoi(argv[3]);
            if(argc>5){
                CONST_P = std::atof(argv[4]);
                PENALTY = std::atof(argv[5]);
            }
        }
    }
    int Nodecount(dim1*dim2);
    gridfile.open(savelocation+"grid.txt");
    elementfile.open(savelocation+"elements.txt");
    datafile.open(savelocation+"data.txt");
    parameterfile.open(savelocation+"parameters.txt");
    resultfile.open(savelocation+"result.txt");

    DynamicVector solution(3*dim1*dim2);
    solution.setConstant(0);
    DynPlaneTriangleFem<3> Model(dim1,dim2,0,1,0,1);
    Model.Build_Vector(Vector_Map);
    Model.Build_Matrix(Matrix_Map);
    Model.Set_Boundary(solution,btype);
    solution += Model.SparseLU();

    for(int i=0;i<dim1*dim2;i++){
        datafile<<Model.grid.nodes[i].eval()+
        " "+std::to_string(solution(i))+
        " "+std::to_string(solution(i+Nodecount))+
        " "+std::to_string(solution(i+2*Nodecount))+
        "\n";
    }
    gridfile<<Model.grid.evalnodes();
    elementfile<<Model.grid.evalelements();
    parameterfile<<dim1<<" "<<dim2<<" "<<CONST_P<<" "<<PENALTY<<"\n";
    resultfile<<solution.head(Nodecount).maxCoeff()<<"\n";
    gridfile.close();
    elementfile.close();
    datafile.close();
    parameterfile.close();
    resultfile.close();
    return 0;
}