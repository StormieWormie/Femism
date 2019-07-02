#include <iostream>
#include <fstream>

#include "../../Library/Tools.hpp"

int main(int argc, char** argv){
    std::ofstream parameterfile;
    parameterfile.open("parameters.txt");
    int count = 0;
    int n = 10;
    int a = 2;
    for(int i=1;i<5;i++){
            parameterfile<<"Dataset"<<count<<" "<<n<<" "<<n<<" "<<5<<"\n";
            count++;
            n = a*n;
    }
    parameterfile.close();
    return 0;
}