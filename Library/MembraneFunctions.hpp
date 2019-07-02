#ifndef MEMBRANE_HPP
#define MEMBRANE_HPP

#include "Tools.hpp"

//Constants(ish)
double
    ALPHA,
    GAMMA,
    LAMBDA,
    SIGMA,
    PENALTY;

//Globals for the PHD
const int DIMENSIONS = 2;
double
    GLOBAL_a(0),
    GLOBAL_b(0),
    GLOBAL_c0(0),
    GLOBAL_c1(0),
    GLOBAL_d0(0),
    GLOBAL_d1(0),
    GLOBAL_e(0),
    GLOBAL_e0(0),
    GLOBAL_e1(0),
    GLOBAL_f(0),
    GLOBAL_g0(0),
    GLOBAL_g1(0),
    GLOBAL_h0(0),
    GLOBAL_h1(0),
    GLOBAL_i(0),
    GLOBAL_j(0),
    GLOBAL_k0(0),
    GLOBAL_k1(0),
    GLOBAL_l0(0),
    GLOBAL_l1(0),
    GLOBAL_m(0);
    

//general terms
DynamicMatrix square(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

DynamicMatrix dx0square(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[i]*TR.b[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix dx1square(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.c[i]*TR.c[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix yidx0j(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix yjdx0i(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}
DynamicMatrix yidx1j(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix yjdx1i(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.c[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix dx0jdx1i(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[j]*TR.c[i]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix dx0idx1j(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[i]*TR.c[j]*TR.Area;
        }
    }
    return result;
}

DynamicVector single(Element<3> TR){
    DynamicVector result(3);
    for(int i=0;i<3;i++){
        result(i) = TR.Area/3.0;
    }
    return result;
}

DynamicVector dx0(Element<3> TR){
    DynamicVector result(3);
    for(int i=0;i<3;i++){
        result(i) = TR.b[i]*TR.Area;
    }
    return result;
}

DynamicVector dx1(Element<3> TR){
    DynamicVector result(3);
    for(int i=0;i<3;i++){
        result(i) = TR.c[i]*TR.Area;
    }
    return result;
}
//George problem

DynamicMatrix GUU(Element<3> TR){    
    return GLOBAL_b*square(TR)+GLOBAL_d0*dx0square(TR)+GLOBAL_d1*dx1square(TR);
}

DynamicMatrix GUV(Element<3> TR){
    return GLOBAL_l0*yidx0j(TR);
}

DynamicMatrix GUW(Element<3> TR){
    return GLOBAL_l1*yidx1j(TR);
}

DynamicMatrix GUF(Element<3> TR){
    return GLOBAL_i*square(TR);
}

DynamicMatrix GVU(Element<3> TR){
    return GLOBAL_l0*yjdx0i(TR);
}

DynamicMatrix GVV(Element<3> TR){
    return (GLOBAL_e+GLOBAL_e0)*dx0square(TR);
}

DynamicMatrix GVW(Element<3> TR){
    return GLOBAL_e*dx0idx1j(TR);
}

DynamicMatrix GVF(Element<3> TR){
    return GLOBAL_k0*yjdx0i(TR);
}

DynamicMatrix GWU(Element<3> TR){
    return GLOBAL_l1*yjdx1i(TR);
}

DynamicMatrix GWV(Element<3> TR){
    return GLOBAL_e*dx0jdx1i(TR);
}

DynamicMatrix GWW(Element<3> TR){
    return (GLOBAL_e+GLOBAL_e1)*dx1square(TR);
}

DynamicMatrix GWF(Element<3> TR){
    return GLOBAL_k1*yjdx1i(TR);
}

DynamicMatrix GFU(Element<3> TR){
    return GLOBAL_i*square(TR);
}

DynamicMatrix GFV(Element<3> TR){
    return GLOBAL_k0*yidx0j(TR);
}

DynamicMatrix GFW(Element<3> TR){
    return GLOBAL_k1*yidx1j(TR);
}

DynamicMatrix GFF(Element<3> TR){    
    return GLOBAL_f*square(TR)+GLOBAL_h0*dx0square(TR)+GLOBAL_h1*dx1square(TR);
}

DynamicVector GU(Element<3> TR){
    return GLOBAL_a*single(TR);
}

DynamicVector GV(Element<3> TR){
    return GLOBAL_c0*dx0(TR);
}

DynamicVector GW(Element<3> TR){
    return GLOBAL_c1*dx1(TR);
}

DynamicVector GF(Element<3> TR){
    return GLOBAL_m*single(TR);
}

ElementToMatrix U4F2_MatrixMap[4][4]{
    {GUU, GUV, GUW, GUF},
    {GVU, GVV, GVW, GVF},
    {GWU, GWV, GWW, GWF},
    {GFU, GFV, GFW, GFF}
};

ElementToVector U4F2_VectorMap[4]{GU, GV, GW, GF};

ElementToMatrix U4_MatrixMap[3][3]{
    {GUU, GUV, GUW},
    {GVU, GVV, GVW},
    {GWU, GWV, GWW}
};

ElementToVector U4_VectorMap[3]{GU, GV, GW};

ElementToMatrix U2F2_MatrixMap[2][2]{
    {GUU, GUF},
    {GFU, GFF}
};

ElementToVector U2F2_VectorMap[2]{GU,GF};

ElementToMatrix U2_MatrixMap[1][1]{{GUU}};

ElementToVector U2_VectorMap[1]{GU};

//ElementToMatrix** Matrixbuilder(double a, doub)
//unambiguous Matrix terms:

DynamicMatrix NAUU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = (1-SIGMA)*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}





DynamicMatrix NAUF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAVV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[i]*TR.b[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAVW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.b[i]*TR.c[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAWV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.c[i]*TR.b[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAWW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = TR.c[i]*TR.c[j]*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAFU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

DynamicMatrix NAFF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            //TODO 
            result(i,j) = (GAMMA*(TR.b[i]*TR.b[j]+TR.c[i]+TR.c[j])+(ALPHA+LAMBDA*LAMBDA)*psipsi(i,j))*TR.Area;
        }
    }
    return result;
}
//PENALTY matrix! to  link u to v and w
DynamicMatrix PUU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = PENALTY*(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}

DynamicMatrix PUV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -PENALTY*TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix PUW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -PENALTY*TR.c[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix PVU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -PENALTY*TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix PVV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = PENALTY*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

DynamicMatrix PWU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -PENALTY*TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix PWW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = PENALTY*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

//Ambiguous term 1: -u L(u)(sigma-2)
DynamicMatrix A1aUV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -0.5*(SIGMA-2)*TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A1aUW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -0.5*(SIGMA-2)*TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A1aVU(Element<3> TR){
    DynamicMatrix result(3,3);;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -0.5*(SIGMA-2)*TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A1aWU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j)= -0.5*(SIGMA-2)*TR.c[i]*TR.Area/3.0;;
        }
    }
    return result;
}

DynamicMatrix A1baUU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = (SIGMA-2)*(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}

DynamicMatrix A1bbVV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = (SIGMA-2)*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

DynamicMatrix A1bbWW(Element<3> TR){
    DynamicMatrix result(3,3);;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = (SIGMA-2)*psipsi(i,j)*TR.Area;
        }
    }
    return result;
}

//Ambiguous term 2: LAMBDA phi L(u)
DynamicMatrix A2aVF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2aWF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*TR.c[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2aFV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2aFW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = LAMBDA*TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2baUF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA *(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}

DynamicMatrix A2baFU(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA *(TR.b[i]*TR.b[j]+TR.c[i]*TR.c[j])*TR.Area;
        }
    }
    return result;
}

DynamicMatrix A2bbVF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA*TR.b[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2bbWF(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA*TR.c[j]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2bbFV(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA*TR.b[i]*TR.Area/3.0;
        }
    }
    return result;
}

DynamicMatrix A2bbFW(Element<3> TR){
    DynamicMatrix result(3,3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            result(i,j) = -LAMBDA*TR.c[i]*TR.Area/3.0;
        }
    }
    return result;
}

//Matrix Accumulation:
ElementToMatrix Matrix_Map_NA[4][4] = {
    {NAUU, Nought, Nought, NAUF},
    {Nought, NAVV, NAVW, Nought},
    {Nought, NAWV, NAWW, Nought},
    {NAFU, Nought, Nought, NAFF}
};

ElementToMatrix Matrix_Map_P[4][4] = {
    {PUU, PUV, PUW, Nought},
    {PVU, PVV, Nought, Nought},
    {PWU, Nought, PWW, Nought},
    {Nought, Nought, Nought, Nought}
};

ElementToMatrix Matrix_Map_PGeorge[3][3] = {
    {PUU, PUV, PUW},
    {PVU, PVV, Nought},
    {PWU, Nought, PWW}
};
//Ambiguity 1
ElementToMatrix Matrix_Map_A1a[4][4]{
    {Nought, A1aUV, A1aUW, Nought},
    {A1aVU, Nought, Nought, Nought},
    {A1aWU, Nought, Nought, Nought},
    {Nought, Nought, Nought, Nought}
};
ElementToMatrix Matrix_Map_A1ba[4][4]{
    {A1baUU, Nought, Nought, Nought},
    {Nought, Nought, Nought, Nought},
    {Nought, Nought, Nought, Nought},
    {Nought, Nought, Nought, Nought}
};
ElementToMatrix Matrix_Map_A1bb[4][4]{
    {Nought, Nought, Nought, Nought},
    {Nought, A1bbVV, Nought, Nought},
    {Nought, Nought, A1bbWW, Nought},
    {Nought, Nought, Nought, Nought}
};

//Ambiguity 2
ElementToMatrix Matrix_Map_A2a[4][4]{
    {Nought, Nought, Nought, Nought},
    {Nought, Nought, Nought, A2aVF},
    {Nought, Nought, Nought, A2aWF},
    {Nought, A2aFV, A2aFW, Nought}
};
ElementToMatrix Matrix_Map_A2ba[4][4]{
    {Nought, Nought, Nought, A2baUF},
    {Nought, Nought, Nought, Nought},
    {Nought, Nought, Nought, Nought},
    {A2baFU, Nought, Nought, Nought}
};
ElementToMatrix Matrix_Map_A2bb[4][4]{
    {Nought, Nought, Nought, Nought},
    {Nought, Nought, Nought, A2bbVF},
    {Nought, Nought, Nought, A2bbWF},
    {Nought, A2bbFV, A2bbFW, Nought}
};

//Hidden part
double PHI;


#endif //MEMBRANE_HPP