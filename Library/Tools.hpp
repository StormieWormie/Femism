#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <iostream>
#include <fstream>
#include <Eigen/Sparse>

#include <vector>

#include "Discretization.hpp"

#define BUG print("BUZZ!");
//Useful functions
const bool inteven(int i){
    return (i==2*(i/2));
}
double psipsi(int i, int j){
    if(i==j){
        return 1.0/6.0;
    }
    else{
        return 1.0/12.0;
    }
}
template<typename T>
void print(T in){
    std::cout<<in<<"\n";
}
template<typename T, typename ... Ts>
void print(T in, Ts ... ins){
    print(in);print(ins...);
}

template<typename T>
T dynexp(T base, int exp){
    if(exp==0){
        return 1;
    }
    else{
        return base*dynexp(base,exp-1);
    }
}

std::vector<double> linspaceparameter(double min, double max, int size){
    std::vector<double> result(size);
    double step = (max-min)/(size-1);
    for(int i=0;i<size;i++){
        result[i] = min+i*step;
    }
    return result;
}

std::vector<double> expspaceparameter(double base, int min, int max){
    std::vector<double> result(max-min+1);
    double val(1);
    if(min<0){
        //This can be done smarter but I dont have time right now
        for(int i=0;i>min;i--){
            val /= base;
        }
        for(int i=0;i<max-min+1;i++){
            result[i] = val;
            val*=base;
        }
    }
    else{
        for(int i=0;i<min;i++){
            val *= base;
        }
        for(int i=0;i<max-min+1;i++){
            result[i] = val;
            val*=base;
        }
    }
    return result;
}

//Type definitions
typedef Eigen::VectorXd                 DynamicVector;
typedef Eigen::VectorXi                 DynamicVectori;
typedef Eigen::MatrixXd                 DynamicMatrix;
typedef Eigen::SparseMatrix<double>     SparseMatrix;
typedef Eigen::SparseMatrix<int>        SparseMatrixi;

//functionals
typedef DynamicMatrix (*ElementToMatrix) (Element<3>);
typedef DynamicVector (*ElementToVector) (Element<3>);
/*Put in Element to data functions!*/


//Empty matrix creation
DynamicMatrix Nought(Element<3> TR){
    DynamicMatrix result(3,3);
    result.setZero();
    return result;
}
//Empty vector creation
DynamicVector Vought(Element<3> TR){
    DynamicVector result(3);
    result.setZero();
    return result;
}

#endif //TOOLS_HPP