#ifndef POISSON_HPP
#define POISSON_HPP

#include "Tools.hpp"

//Constants
double CONST_F;

//Global checker
void global_checker(){
    std::cout<<"F: "<<CONST_F<<"\n";
}


//Element vectors
DynamicVector Element_Vector(Element<3> TR){
    DynamicVector result(3);
    for(int i=0;i<3;i++){
        result(i) = CONST_F*TR.Area/3.0;
    }
    return result;
}
//Element vectors
DynamicMatrix Element_Matrix(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}

//Mapping collective
ElementToMatrix Matrix_Map[1][1]{
    {Element_Matrix}
};
ElementToVector Vector_Map[1]{
    Element_Vector
};
#endif