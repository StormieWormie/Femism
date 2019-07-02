#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../Library/Tools.hpp"

int main(){
    

    
    std::vector<double>
    dim1({5, 10}),
    sigma(linspaceparameter(0.45,0.55,11)),
    penalty({1000, 1000000});

    
    int iter(0);
    std::ofstream parameterfile;
    parameterfile.open("parameters.txt");
    parameterfile << std::setfill('0');
    for(int d1 : dim1){
    for(double s : sigma){
    for(double p : penalty){

        parameterfile<<
        "dataset"<<iter<<" "<<
        d1<<" "<<(int)13.78*d1+1<<" "<<
        p<<" "<<s
        <<"\n";

        iter++;
    }}}

    parameterfile.close();
    
    return 0;
}