#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "Tools.hpp"

//Constants
double CONST_P,
PENALTY;

//Global checker
void global_checker(){
    std::cout<<"F: "<<CONST_P<<"\n";
}
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

//Element vectors
DynamicVector Fu(Element<3> TR){
    DynamicVector result(3);
    for(int i=0;i<3;i++){
        result(i) = CONST_P*TR.Area/3.0;
    }
    return result;
}
//Element vectors
DynamicMatrix Muu(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=PENALTY*(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}
DynamicMatrix Muv(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=-PENALTY*TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}
DynamicMatrix Muw(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=-PENALTY*TR.c[i]*TR.Area/3.0;
        }
    }
    return result;
}
DynamicMatrix Mvu(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=-PENALTY*TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}
DynamicMatrix Mvv(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=(TR.b[i]*TR.b[j]+PENALTY*psipsi(i,j))*TR.Area;
        }
    }
    return result;
}
DynamicMatrix Mvw(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=TR.c[j]*TR.b[i]*TR.Area;
        }
    }
    return result;
}
DynamicMatrix Mwu(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=-PENALTY*TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}
DynamicMatrix Mwv(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=TR.b[j]*TR.c[i]*TR.Area;
        }
    }
    return result;
}
DynamicMatrix Mww(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)=(TR.b[i]*TR.b[j]+PENALTY*psipsi(i,j))*TR.Area;
        }
    }
    return result;
}


//Mapping collective
ElementToMatrix Matrix_Map[3][3]{
    {Muu,Muv,Muw},
    {Mvu,Mvv,Mvw},
    {Mwu,Mwv,Mww}
};
ElementToVector Vector_Map[3]{
    Fu,
    Vought,
    Vought
};
#endif //FUNCTIONS_HPP