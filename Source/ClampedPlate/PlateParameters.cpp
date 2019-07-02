#include <iostream>
#include <fstream>

#include "../../Library/Tools.hpp"

int main(int argc, char** argv){
    std::ofstream parameterfile;
    parameterfile.open("parameters.txt");
    int count = 0;
    int n = 10;
    int a = 2;
    int penalty = 10;
    int pinc = 10;
    
    for(int i=1;i<2;i++){
        for(int j=0;j<30;j++){
            parameterfile<<"Dataset"<<count
            <<" "<<n*dynexp(a,i)
            <<" "<<n*dynexp(a,i)
            <<" "<<10
            <<" "<<dynexp(2,j)
            <<"\n";
            count++;
        }
    }
    
    parameterfile.close();
    return 0;
}