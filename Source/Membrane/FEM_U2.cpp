#include "../../Library/FEM.hpp"

#include "../../Library/MembraneFunctions.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(10),dim2(10);
    const int funccount = 1;
    int btype[funccount] = {0};
    std::string savelocation("./");
    if(argc>3){
        savelocation = argv[1];
        dim1 = std::atoi(argv[2]);
        dim2 = std::atoi(argv[3]);
        GLOBAL_a = std::atof(argv[4]);
        GLOBAL_b = std::atof(argv[5]);
        GLOBAL_d0 = std::atof(argv[6]);
        GLOBAL_d1 = std::atof(argv[7]);
        print("input accepted");
    }
    else{
        exit(0);
    }
    
    gridfile.open(savelocation+"/grid.txt");
    elementfile.open(savelocation+"/elements.txt");
    datafile.open(savelocation+"/data.txt");
    parameterfile.open(savelocation+"/parameters.txt");
    resultfile.open(savelocation+"/result.txt");

    DynamicVector solution(funccount*dim1*dim2);
    
    
    solution.setConstant(0);
    print("Setting up the system");
    DynPlaneTriangleFem<funccount> Model(dim1,dim2,0,1,0,1);
    print("Building problem matrix");
    Model.Build_Matrix(U2_MatrixMap);
    print("Building problem vector");
    Model.Build_Vector(U2_VectorMap);
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
        <<solution(i)<<"\n";
    }
    
    
    gridfile<<Model.grid.evalnodes();
    elementfile<<Model.grid.evalelements();
    parameterfile<<
    "location: "<<savelocation<<
    " dim1: "<<dim1<<
    " dim2: "<<dim2<<
    " Penalty: "<<PENALTY<<
    " a: "<<GLOBAL_a<<
    " b: "<<GLOBAL_b<<
    " d0: "<<GLOBAL_d0<<
    " d1: "<<GLOBAL_d1<<
    "\n";
    resultfile;
    gridfile.close();
    elementfile.close();
    datafile.close();
    parameterfile.close();
    resultfile.close();
    return 0;
}