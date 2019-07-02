#include "../../Library/FEM.hpp"

#include "../../Library/MembraneFunctions.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(10),dim2(10);
    const int funccount = 3;
    int btype[funccount] = {0,0,0};
    std::string savelocation("./");
    PENALTY = 1000;
    double sigma = 0.5;
    if(argc>3){
        savelocation = argv[1];
        dim1 = std::atoi(argv[2]);
        dim2 = std::atoi(argv[3]);
        PENALTY = std::atof(argv[4]);
        sigma = std::atof(argv[5]);
        print("input accepted");
    }
    GLOBAL_d0 = 0.5*(2*sigma-1)+1;
    GLOBAL_d1 = 0.5*(2*sigma-1);
    GLOBAL_e = 4;
    GLOBAL_b = 1;
    GLOBAL_c0 = -1;
    GLOBAL_l0 = 2;
    GLOBAL_a = 0.5*(2*sigma-1);
    GLOBAL_c1 = -1;
    gridfile.open(savelocation+"/grid.txt");
    elementfile.open(savelocation+"/elements.txt");
    datafile.open(savelocation+"/data.txt");
    parameterfile.open(savelocation+"/parameters.txt");
    resultfile.open(savelocation+"/result.txt");

    DynamicVector solution(funccount*dim1*dim2);
    
    
    solution.setConstant(0);
    print("Setting up the system");
    DynCylinderTriangleFem<funccount> Model(dim1,dim2,0,100);
    print("Building problem matrix");
    Model.Build_Matrix(U4_MatrixMap);
    print("Building penalty matrix");
    Model.Build_Matrix(Matrix_Map_PGeorge);
    print("Building problem vector");
    Model.Build_Vector(U4_VectorMap);
    /*
    double offset1(0.25), offset2(0.2*M_PI), psi1, psi2, dtheta(2*M_PI/double(dim1));
    print("Setting boundaries");
    //Set U to elliptical
    
    for(int i=0;i<dim1;i++){
        psi1 = atan(ratio*tan(i*dtheta+offset1));
        psi2 = atan(ratio*tan(i*dtheta+offset2));
        solution(Model.grid.external(i)) = sqrt(ratio*pow(cos(psi1),2)+ 1/ratio * pow(sin(psi1),2)) - 1;
        solution(Model.grid.external(dim1+i)) = sqrt(ratio*pow(cos(psi2),2)+ 1/ratio * pow(sin(psi2),2)) - 1;
    }
    solution(dim1*dim2+Model.grid.external(0)) = (solution(Model.grid.external(1))-solution(Model.grid.external(dim1-1)))/(2*dtheta);
    solution(dim1*dim2+Model.grid.external(dim1-1)) = (solution(Model.grid.external(0))-solution(Model.grid.external(dim1-2)))/(2*dtheta);
    solution(dim1*dim2+Model.grid.external(dim1)) = (solution(Model.grid.external(dim1+1))-solution(Model.grid.external(2*dim1-1)))/(2*dtheta);
    solution(dim1*dim2+Model.grid.external(2*dim1-1)) = (solution(Model.grid.external(dim1))-solution(Model.grid.external(2*dim1-2)))/(2*dtheta);
    
    for(int i=1;i<dim1-1;i++){
        solution(dim1*dim2+Model.grid.external(i)) = (solution(Model.grid.external(i+1))-solution(Model.grid.external(i-1)))/(2*dtheta);
        solution(dim1*dim2+Model.grid.external(dim1+i)) = (solution(Model.grid.external(dim1+i+1))-solution(Model.grid.external(dim1+i-1)))/(2*dtheta);
    }
    

    for(int i=0;i<Model.grid.external.size();i++){

        // SET VALUE FOR PROTEINS ON BOUNDARY
        // solution(3*dim1*dim2+Model.grid.external(i)) = 1;
    }
    */
    
    //print("buzz0");
    print("Psych! actually building boundary");
    Model.Set_Boundary(solution, btype);
    //print("buzz1");
    //print("solsize",solution.size());
    print("Solving...");
    solution += Model.SparseLU();
    //print("solvedsize",solutionadd.size());
    //print("buzz2");
    for(int i=0;i<dim1*dim2;i++){
        datafile<<Model.grid.nodes[i].eval()+" "
        <<solution(i)<<" "
        <<solution(i+dim1*dim2)<<" "
        <<solution(i+2*dim1*dim2)<<"\n";
    }
    
    
    gridfile<<Model.grid.evalnodes();
    elementfile<<Model.grid.evalelements();
    parameterfile<<
    "location: "<<savelocation<<
    " dim1: "<<dim1<<
    " dim2: "<<dim2<<
    " sigma "<<sigma<<
    " Penalty: "<<PENALTY<<
    "\n";
    resultfile;
    gridfile.close();
    elementfile.close();
    datafile.close();
    parameterfile.close();
    resultfile.close();
    return 0;
}