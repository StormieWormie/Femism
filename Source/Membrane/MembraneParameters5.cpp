#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../Library/Tools.hpp"

int main(){
    

    
    std::vector<double>
    sigma(linspaceparameter(0.45,0.55,11)),
    penalty({1000, 1000000});

    
    int iter(0);
    std::ofstream parameterfile;
    parameterfile.open("parameters.txt");
    parameterfile << std::setfill('0');
    for(double s : sigma){
        for(double p : penalty){

            parameterfile<<
            s<<" "<<p
            <<"\n";

            iter++;
        }
    }

    parameterfile.close();
    
    return 0;
}