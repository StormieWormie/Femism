#include "../../Library/FEM.hpp"

#include "../../Library/Poisson.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(10),dim2(10);
    int btype[1] = {0};
    std::string savelocation("./");
    CONST_F = 0.1;
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
    gridfile.open(savelocation+"/grid.txt");
    elementfile.open(savelocation+"/elements.txt");
    datafile.open(savelocation+"/data.txt");
    parameterfile.open(savelocation+"/parameters.txt");
    resultfile.open(savelocation+"/result.txt");

    DynamicVector solution(dim1*dim2);
    solution.setConstant(0);
    DynPlaneTriangleFem<1> Model(dim1,dim2,0,1,0,1);
    Model.Build_Vector(Vector_Map);
    Model.Build_Matrix(Matrix_Map);
    
    //print("first",Model.grid.internal);
    //print("out",Model.grid.external);
    DynamicVectori tempinternal((dim1-2)*(dim2-1));
    DynamicVectori tempexternal(2*dim2+dim1-2);
    tempinternal.head((dim1-2)*(dim2-2)) = Model.grid.internal;
    tempinternal.tail(dim1-2) = Model.grid.external.tail(dim1-1).head(dim1-2);
    tempexternal.head((dim1-1)+2*(dim2-1)) = Model.grid.external.head((dim1-1)+2*(dim2-1));
    tempexternal.tail(1) = Model.grid.external.tail(1);
    //Model.grid.internal.resize(Model.grid.internal.size()+dim1-2);
    //print("tempin",tempinternal);
    //print("tempout",tempexternal);
    Model.grid.internal = tempinternal;
    Model.grid.external = tempexternal;
    //print("inafter",Model.grid.internal);
    //print("outafter",Model.grid.external);
    //print("done");
    Model.Set_Boundary(solution,btype);
    //print(Model.S.rows(), Model.fS.rows());
    //print(Model.perm.rows());

    solution += Model.SparseLU();
    //print("here?");
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