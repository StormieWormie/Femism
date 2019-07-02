#ifndef FEM_HPP
#define FEM_HPP


#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>


#include "Tools.hpp"
#include "Discretization.hpp"




//Actual code
//Parent class
template<typename Grid_type, int Variable_count,int Node_count, int Element_count, int etype>
class FiniteElementMethod{
public:
//Attributes
    Eigen::Matrix<double,Variable_count*Node_count,1> f;
    DynamicVector ff;
    SparseMatrix S,fS,perm,permT;
    int Vcount, Ncount, Ecount;
    Grid_type grid;
//Constructors

    FiniteElementMethod():f(),ff(),S(Variable_count*Node_count,Variable_count*Node_count),fS(),perm(Variable_count*Node_count,Variable_count*Node_count),permT(Variable_count*Node_count,Variable_count*Node_count),Vcount(Variable_count),Ncount(Node_count),Ecount(Element_count){
        f.setZero();
    }
    /*
    FiniteElementMethod(double LowerZ, double HigherZ):f(),ff(),S(Variable_count*Node_count,Variable_count*Node_count),fS(),perm(Variable_count*Node_count,Variable_count*Node_count),permT(Variable_count*Node_count,Variable_count*Node_count),Vcount(Variable_count),Ncount(Node_count),Ecount(Element_count),grid(LowerZ,HigherZ){
        f.setZero();
    }
    */
    FiniteElementMethod(Grid_type given_grid):f(),ff(),S(Variable_count*Node_count,Variable_count*Node_count),fS(),perm(Variable_count*Node_count,Variable_count*Node_count),permT(Variable_count*Node_count,Variable_count*Node_count),Vcount(Variable_count),Ncount(Node_count),Ecount(Element_count),grid(given_grid){
        f.setZero();
    }
//Methods
    Eigen::Matrix<double,etype*Variable_count,1> get_vector(ElementToVector function[Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,1> result;
        
        for(int i=0;i<Variable_count;i++){
            result.segment(i*etype,etype) = function[i](grid.elements[k]);
        }
        return result;
    }
    Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> get_matrix(ElementToMatrix function[Variable_count][Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> result;
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> get_penalty(ElementToMatrix function[Variable_count][Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> result;
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> get_boundary_matrix(ElementToMatrix function[Variable_count][Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> result;
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    Eigen::Matrix<double,Variable_count*Node_count,1>& Build_Vector(ElementToVector function[Variable_count]){
        Eigen::Matrix<double,etype*Variable_count,1> element_vector;
        for(int e=0;e<Element_count;e++){
            element_vector = get_vector(function, e);
            for(int I=0;I<Variable_count;I++){
                for(int i=0;i<etype;i++){
                    f(grid.elements[e].associate[i]+I*(Node_count))+=element_vector(i+etype*I);
                }
            }
        }
        return f;
    }
    SparseMatrix& Build_Matrix(ElementToMatrix function[Variable_count][Variable_count]){
        S.setZero();
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> element_matrix;
        for(int e=0;e<Element_count;e++){
            element_matrix = get_matrix(function,e);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[e].associate[i]+I*(Node_count),grid.elements[e].associate[j]+J*(Node_count))+=element_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
        }
        return S;
    }
    SparseMatrix& Build_Penalty(ElementToMatrix function[Variable_count][Variable_count]){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> penalty_matrix;
        for(int e=0;e<Element_count;e++){
            penalty_matrix = get_penalty(function,e);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[e].associate[i]+I*(Node_count),grid.elements[e].associate[j]+J*(Node_count))+=penalty_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
        }
        return S;
    }
    SparseMatrix& Build_NeumannB(ElementToMatrix Mfunction[Variable_count][Variable_count], ElementToVector Ffunction[Variable_count]){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> penalty_matrix;
        Eigen::Matrix<double,etype*Variable_count,1> penalty_vector;
        for(int e=0;e<grid.externalE.size();e++){
            penalty_matrix = get_boundary_matrix(Mfunction,grid.externalE[e]);
            penalty_vector = get_vector(Ffunction,grid.externalE[e]);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[grid.externalE[e]].associate[i]+I*Node_count,grid.elements[grid.externalE[e]].associate[j]+J*Node_count) += penalty_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
            for(int I=0;I<Variable_count;I++){
                for(int i=0;i<etype;i++){
                    f(grid.elements[grid.externalE[e]].associate[i]+I*Node_count)+=penalty_vector(i+etype*I);
                }
            }
        }
        return S;
    }
    void AllDirichlet(DynamicVector bc){
        int isize(grid.internal.size());
        int esize(grid.external.size());
        for(int v=0;v<Variable_count;v++){
            for(int i=0;i<isize;i++){
                perm.insert(i+v*isize,grid.internal(i)+v*Node_count)=1;
                permT.insert(grid.internal(i)+v*Node_count,i+v*isize)=1;
            }
        }
        for(int v=0;v<Variable_count;v++){
            for(int i=0;i<esize;i++){
                perm.insert(i+Variable_count*isize+v*esize,grid.external(i)+v*Node_count)=1;
                permT.insert(grid.external(i)+v*Node_count,i+Variable_count*isize+v*esize)=1;
    
            }
        }
        SparseMatrix tempS = perm*S*permT;
        DynamicVector tempf(perm*f), tempbc(perm*bc);
        /*
        print("S:");
        print(S);
        print("perm");
        print(perm);
        print("tempS");
        print(tempS);
        */
        for(int i=0;i<Variable_count*esize;i++){
            for(int j=0;j<Variable_count*isize;j++){
                tempf(j) -= tempbc(i+Variable_count*isize)*tempS.coeffRef(j,i+Variable_count*isize);
            }
        }
        fS = tempS.block(0,0,Variable_count*isize,Variable_count*isize);
        ff = tempf.head(Variable_count*isize);
    }
    SparseMatrix& Build_Perm(DynamicVector bc, int btype[Variable_count]){
        /*
        Boundary types are as follows:
        Dirichlet: 0
        Neumann: 1
        Robin: 2 !!!not programmed in!!!
        */
        int isize(grid.internal.size());
        int esize(grid.external.size());
        int offset(0);
        for(int v=0;v<Variable_count;v++){
            if(btype[v]==0){
                for(int i=0;i<isize;i++){
                    perm.insert(i+offset,grid.internal(i)+v*Node_count)=1;
                    permT.insert(grid.internal(i)+v*Node_count,i+offset)=1;
                }
                offset+=isize;
            }
            else if(btype[v]==1){
                for(int i=0;i<Node_count;i++){
                    perm.insert(i+offset,i+v*Node_count)=1;
                    permT.insert(i+v*Node_count,i+offset)=1;
                }
                offset+=Node_count;
            }
        }
        int Isize = offset;
        for(int v=0;v<Variable_count;v++){
            if(btype[v]==0){
                for(int i=0;i<esize;i++){
                    perm.insert(i+offset,grid.external(i)+v*Node_count)=1;
                    permT.insert(grid.external(i)+v*Node_count,i+offset)=1;
                }
                offset+=esize;
            }
            else if(btype[v]==1){
            }
        }
        int Esize = offset-Isize;

        SparseMatrix tempS = perm*S*permT;
        DynamicVector tempf(perm*f), tempbc(perm*bc);
        /*
        print("bc");
        print(bc);
        print("tempbc");
        print(tempbc);
        */
       /*
        print("Isize: ");
        print(Isize);
        print("Esize: ");
        print(Esize);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        */
        //Setting Dirichlet conditions
        for(int i=0;i<Esize;i++){
            for(int j=0;j<Isize;j++){
                tempf(j) -= tempbc(i+Isize)*tempS.coeffRef(j,i+Isize);
            }
        }

        fS = tempS.block(0,0,Variable_count*isize,Variable_count*isize);
        ff = tempf.head(Variable_count*isize);
        return perm;
    }



    template<typename Solver>
    DynamicVector Solve(){
        Solver solver;
        DynamicVector Solution(Variable_count*grid.internal.size()),result(Variable_count*Node_count);
        result.tail(Variable_count*grid.external.size()).setZero();
        solver.compute(fS);
        if(solver.info()!=Eigen::Success){
            print("Decomposition failed!");
            return result;
        }
        Solution = solver.solve(ff);
        if(solver.info()!=Eigen::Success){
            print("Solving failed!");
            return result;
        }
        print("Solve complete!");
        result.head(Variable_count*grid.internal.size()) = Solution;
        return permT*result;
    }

    DynamicVector SparseLU(){
        return Solve<Eigen::SparseLU<SparseMatrix>>();
    }

    DynamicVector BiCGSTAB(){
        return Solve<Eigen::BiCGSTAB<SparseMatrix>>();
    }

};

template<int Variable_count, int dim1, int dim2>
class PlaneTriangleFem : public FiniteElementMethod<TriangulatedGrid<dim1,dim2>,Variable_count,dim1*dim2, 2*(dim1-1)*(dim2-1), 3>{
public:
    PlaneTriangleFem(double LowerX, double HigherX, double LowerY, double HigherY):FiniteElementMethod<TriangulatedGrid<dim1,dim2>,Variable_count,dim1*dim2, 2*(dim1-1)*(dim2-1), 3>(){
        this->grid = TriangulatedGrid<dim1,dim2>(LowerX,HigherX,LowerY,HigherY);
    }
};

template<int Variable_count, int dim1, int dim2>
class CylinderTriangleFem : public FiniteElementMethod<TriangulatedCylinder<dim1,dim2>,Variable_count, dim1*dim2, 2*dim1*(dim2-1), 3>{
public:
//Attributes
//Constructors

    CylinderTriangleFem(double LowerZ, double HigherZ):FiniteElementMethod<TriangulatedCylinder<dim1,dim2>,Variable_count, dim1*dim2, 2*dim1*(dim2-1), 3>(){
        this->grid = TriangulatedCylinder<dim1,dim2>(LowerZ,HigherZ);
    }

};


template<typename Grid_type, int Variable_count, int etype>
class DynFiniteElementMethod{
public:
//Attributes
    int Node_count, Element_count, Isize;
    DynamicVector f, ff;
    SparseMatrix S,fS,perm,permT;
    Grid_type grid;
//Constructors
    DynFiniteElementMethod(int n, int e):
    Node_count(n),Element_count(e),Isize(0),
    f(Variable_count*Node_count),ff(),
    S(Variable_count*Node_count,Variable_count*Node_count),fS(),perm(Variable_count*Node_count,Variable_count*Node_count),permT(Variable_count*Node_count,Variable_count*Node_count)
    {f.setZero();S.setZero();}
//Methods
    DynamicVector get_vector(ElementToVector function[Variable_count], int k){
        DynamicVector result(Variable_count*etype);
        for(int i=0;i<Variable_count;i++){
            result.segment(i*etype,etype) = function[i](grid.elements[k]);
        }
        return result;
    }
    DynamicMatrix get_matrix(ElementToMatrix function[Variable_count][Variable_count], int k){
        DynamicMatrix result(Variable_count*etype,Variable_count*etype);
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    /*
    Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> get_penalty(ElementToMatrix function[Variable_count][Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> result;
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    */
    /*
    Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> get_boundary_matrix(ElementToMatrix function[Variable_count][Variable_count], int k){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> result;
        for(int i=0;i<Variable_count;i++){
            for(int j=0;j<Variable_count;j++){
                result.block(i*etype,j*etype,etype,etype) = function[i][j](grid.elements[k]);
            }
        }
        return result;
    }
    */
    DynamicVector& Build_Vector(ElementToVector function[Variable_count]){
        DynamicVector element_vector(Variable_count*etype);
        for(int e=0;e<Element_count;e++){
            element_vector = get_vector(function, e);
            for(int I=0;I<Variable_count;I++){
                for(int i=0;i<etype;i++){
                    f(grid.elements[e].associate[i]+I*(Node_count))+=element_vector(i+etype*I);
                }
            }
        }
        return f;
    }
    SparseMatrix& Build_Matrix(ElementToMatrix function[Variable_count][Variable_count]){
        DynamicMatrix element_matrix(Variable_count*etype,Variable_count*etype);
        for(int e=0;e<Element_count;e++){
            element_matrix = get_matrix(function,e);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[e].associate[i]+I*(Node_count),grid.elements[e].associate[j]+J*(Node_count))+=element_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
        }
        return S;
    }
    /*
    SparseMatrix& Build_Penalty(ElementToMatrix function[Variable_count][Variable_count]){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> penalty_matrix;
        for(int e=0;e<Element_count;e++){
            penalty_matrix = get_penalty(function,e);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[e].associate[i]+I*(Node_count),grid.elements[e].associate[j]+J*(Node_count))+=penalty_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
        }
        return S;
    }
    */
    /*
    SparseMatrix& Build_NeumannB(ElementToMatrix Mfunction[Variable_count][Variable_count], ElementToVector Ffunction[Variable_count]){
        Eigen::Matrix<double,etype*Variable_count,etype*Variable_count> penalty_matrix;
        Eigen::Matrix<double,etype*Variable_count,1> penalty_vector;
        for(int e=0;e<grid.externalE.size();e++){
            penalty_matrix = get_boundary_matrix(Mfunction,grid.externalE[e]);
            penalty_vector = get_vector(Ffunction,grid.externalE[e]);
            for(int I=0;I<Variable_count;I++){
                for(int J=0;J<Variable_count;J++){
                    for(int i=0;i<etype;i++){
                        for(int j=0;j<etype;j++){
                            S.coeffRef(grid.elements[grid.externalE[e]].associate[i]+I*Node_count,grid.elements[grid.externalE[e]].associate[j]+J*Node_count) += penalty_matrix(i+etype*I,j+etype*J);
                        }
                    }
                }
            }
            for(int I=0;I<Variable_count;I++){
                for(int i=0;i<etype;i++){
                    f(grid.elements[grid.externalE[e]].associate[i]+I*Node_count)+=penalty_vector(i+etype*I);
                }
            }
        }
        return S;
    }
    */
    /*
    void AllDirichlet(DynamicVector bc){
        int isize(grid.internal.size());
        int esize(grid.external.size());
        for(int v=0;v<Variable_count;v++){
            for(int i=0;i<isize;i++){
                perm.insert(i+v*isize,grid.internal(i)+v*Node_count)=1;
                permT.insert(grid.internal(i)+v*Node_count,i+v*isize)=1;
            }
        }
        for(int v=0;v<Variable_count;v++){
            for(int i=0;i<esize;i++){
                perm.insert(i+Variable_count*isize+v*esize,grid.external(i)+v*Node_count)=1;
                permT.insert(grid.external(i)+v*Node_count,i+Variable_count*isize+v*esize)=1;
    
            }
        }
        SparseMatrix tempS = perm*S*permT;
        DynamicVector tempf(perm*f), tempbc(perm*bc);
        
        for(int i=0;i<Variable_count*esize;i++){
            for(int j=0;j<Variable_count*isize;j++){
                tempf(j) -= tempbc(i+Variable_count*isize)*tempS.coeffRef(j,i+Variable_count*isize);
            }
        }
        fS = tempS.block(0,0,Variable_count*isize,Variable_count*isize);
        ff = tempf.head(Variable_count*isize);
    }
    */
    SparseMatrix& Set_Boundary(DynamicVector bc, int btype[Variable_count]){
        /*
        Boundary types are as follows:
        Dirichlet: 0
        Neumann: 1
        Robin: 2 !!!not programmed in!!!
        */
        int isize(grid.internal.size());
        int esize(grid.external.size());
        int offset(0);
        //print("BUZZ0");
        for(int v=0;v<Variable_count;v++){
            if(btype[v]==0){
                for(int i=0;i<isize;i++){
                    perm.insert(i+offset,grid.internal(i)+v*Node_count)=1;
                    permT.insert(grid.internal(i)+v*Node_count,i+offset)=1;
                }
                offset+=isize;
            }
            else if(btype[v]==1){
                for(int i=0;i<Node_count;i++){
                    perm.insert(i+offset,i+v*Node_count)=1;
                    permT.insert(i+v*Node_count,i+offset)=1;
                }
                offset+=Node_count;
            }
        }
        //print("BUZZ1");
        Isize = offset;
        for(int v=0;v<Variable_count;v++){
            if(btype[v]==0){
                for(int i=0;i<esize;i++){
                    perm.insert(i+offset,grid.external(i)+v*Node_count)=1;
                    permT.insert(grid.external(i)+v*Node_count,i+offset)=1;
                }
                offset+=esize;
            }
            else if(btype[v]==1){
            }
        }
        int Esize = offset-Isize;
        //print("BUZZ2");
        SparseMatrix tempS = perm*S*permT;
        //print("hi!");
        //print("Permsize", perm.rows(), perm.cols());
        DynamicVector tempf(perm*f), tempbc(perm*bc);
        //print("afterperms");
        /*
        print("bc");
        print(bc);
        print("tempbc");
        print(tempbc);
        */
       /*
        print("Isize: ");
        print(Isize);
        print("Esize: ");
        print(Esize);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        */
        //Setting Dirichlet conditions
        for(int i=0;i<Esize;i++){
            for(int j=0;j<Isize;j++){
                tempf(j) -= tempbc(i+Isize)*tempS.coeffRef(j,i+Isize);
            }
        }
        //print("Block!");
        //print("Isize",Isize);
        //print("tempS size",tempS.rows(),tempS.cols());
        fS = tempS.block(0,0,Isize,Isize);
        //print("head!");
        ff = tempf.head(Isize);
        //print("BUZZ3");
        return perm;
    }



    template<typename Solver>
    DynamicVector Solve(){
        Solver solver;
        DynamicVector Solution(Variable_count*grid.internal.size()),result(Variable_count*Node_count);
        result.tail(Variable_count*grid.external.size()).setZero();
        solver.compute(fS);
        if(solver.info()!=Eigen::Success){
            print("Decomposition failed!");
            return result;
        }
        Solution = solver.solve(ff);
        if(solver.info()!=Eigen::Success){
            print("Solving failed!");
            return result;
        }
        print("Solve complete!");
        //TODO DEFINITELY NOT 25!!!!!
        //print("solsize",Solution.size());
        //print("Isize",Isize);
        result.head(Isize) = Solution;
        print(Isize);
        return permT*result;
    }

    DynamicVector SparseLU(){
        return Solve<Eigen::SparseLU<SparseMatrix>>();
    }

    DynamicVector BiCGSTAB(){
        return Solve<Eigen::BiCGSTAB<SparseMatrix>>();
    }

    DynamicVector ConGrad(){
        return Solve<Eigen::ConjugateGradient<SparseMatrix, Eigen::Upper>>();
    }

    DynamicVector LDLT(){
        return Solve<Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper>>();
    }

    
};

template<int Variable_count>
class DynPlaneTriangleFem : public DynFiniteElementMethod<DynTriangulatedGrid,Variable_count,3>{
public:
    DynPlaneTriangleFem(int dim1, int dim2, double LowerX, double HigherX, double LowerY, double HigherY):
    DynFiniteElementMethod<DynTriangulatedGrid,Variable_count, 3>(dim1*dim2,2*(dim1-1)*(dim2-1)){
        this->grid = DynTriangulatedGrid(dim1,dim2,LowerX,HigherX,LowerY,HigherY);
    }
};
//todo INPUT OF PARENT
template<int Variable_count>
class DynCylinderTriangleFem : public DynFiniteElementMethod<DynTriangulatedCylinder,Variable_count, 3>{
public:
//Attributes
//Constructors
    DynCylinderTriangleFem(int dim1, int dim2, double LowerZ, double HigherZ):DynFiniteElementMethod<DynTriangulatedCylinder,Variable_count, 3>(dim1*dim2, 2*dim1*(dim2-1)){
        this->grid = DynTriangulatedCylinder(dim1,dim2,LowerZ,HigherZ);
    }
};
template<int Variable_count>
class DynDoughnutTriangleFem : public DynFiniteElementMethod<DynTriangulatedDoughnut, Variable_count, 3>{
public:
    DynDoughnutTriangleFem(int dim1, int dim2, double LowerZ, double HigherZ):
    DynFiniteElementMethod<DynTriangulatedDoughnut, Variable_count, 3>(dim1*dim2, 2*dim1*dim2){
        this->grid = DynTriangulatedDoughnut(dim1,dim2,LowerZ,HigherZ);
    }
};

#endif //FEM_HPP