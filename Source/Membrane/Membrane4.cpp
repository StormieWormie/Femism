#include "../../Library/FEM.hpp"

#include "../../Library/MembraneFunctions.hpp"

int main(int argc, char** argv){
    std::ofstream gridfile, elementfile, datafile, parameterfile, resultfile;
    int dim1(5),dim2(5);
    int btype[4] = {0,0,0,1};
    std::string savelocation("./");
    int A1Choice = 0;
    int A2Choice = 0;
    double ratio(1.00005);
    if(argc>1){
        savelocation = argv[1];
        if(argc>3){
            dim1 = std::atoi(argv[2]);
            dim2 = std::atoi(argv[3]);
            if(argc>4){
                ALPHA = std::atof(argv[4]);
                GAMMA = std::atof(argv[5]);
                LAMBDA = std::atof(argv[6]);
                SIGMA = std::atof(argv[7]);
                PENALTY = std::atof(argv[8]);
                A1Choice = std::atoi(argv[9]);
                A2Choice = std::atoi(argv[10]);
                if(argc>11){
                    ratio = 1.0+std::atof(argv[11]);
                    print(ratio);
                }
            }
        }
    }
    print(A1Choice, A2Choice);
    gridfile.open(savelocation+"/grid.txt");
    elementfile.open(savelocation+"/elements.txt");
    datafile.open(savelocation+"/data.txt");
    parameterfile.open(savelocation+"/parameters.txt");
    resultfile.open(savelocation+"/result.txt");

    DynamicVector solution(4*dim1*dim2);
    
    
    solution.setConstant(0);
    print("Setting up the system");
    DynCylinderTriangleFem<4> Model(dim1,dim2,0,100);
    print("Building a matrix");
    Model.Build_Matrix(Matrix_Map_NA);
    Model.Build_Matrix(Matrix_Map_P);
    switch (A1Choice)
        {
        case 0:
            print("Using A1a\n");
            Model.Build_Matrix(Matrix_Map_A1a);
            break;
        case 1:
            print("Using A1ba\n");
            Model.Build_Matrix(Matrix_Map_A1ba);
            break;
        case 2:
            print("Using A1bb\n");
            Model.Build_Matrix(Matrix_Map_A1bb);
            break;
        
        default:
            break;
        }
    switch (A2Choice)
        {
        case 0:
            print("Using A2a\n");
            Model.Build_Matrix(Matrix_Map_A2a);
            break;
        case 1:
            print("Using A2ba\n");
            Model.Build_Matrix(Matrix_Map_A2ba);
            break;
        case 2:
            print("Using A2bb\n");
            Model.Build_Matrix(Matrix_Map_A2bb);
            break;
        
        default:
            break;
        }

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
        <<solution(i+2*dim1*dim2)<<" "
        <<solution(i+3*dim1*dim2)<<"\n";
    }
    
    
    gridfile<<Model.grid.evalnodes();
    elementfile<<Model.grid.evalelements();
    parameterfile<<
    "location: "<<savelocation<<
    " dim1: "<<dim1<<
    " dim2: "<<dim2<<
    " A1choice: "<<A1Choice<<
    " A2choise: "<<A2Choice<<
    " Alpha: "<<ALPHA<<
    " Gamma: "<<GAMMA<<
    " Lambda: "<<LAMBDA<<
    " Sigma: "<<SIGMA<<
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