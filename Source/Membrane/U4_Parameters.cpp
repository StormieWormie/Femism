#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../Library/Tools.hpp"

int main(){
    

    
    std::vector<double>
    dim({5, 10, 20}),
    A({0,10,20}),
    B({0}),
    DX({1});

    
    int iter(0);
    std::ofstream parameterfile;
    parameterfile.open("parameters.txt");
    parameterfile << std::setfill('0');
    for(int d : dim){
    for(double a : A){
    for(double b : B){
    for(double dx : DX){

        parameterfile<<
        "dataset"<<iter<<" "<<
        d<<" "<<d<<" "<<
        a<<" "<<
        b<<" "<<
        dx<<" "<<
        dx
        <<"\n";

        iter++;
    }}}}

    parameterfile.close();
    
    return 0;
}