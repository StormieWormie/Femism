#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <cmath>
#include "Tools.hpp"

//Individual components for discretization:
class Node{
public:
//Attributes
    int index;
    double First, Second;
    std::string name1, name2;
//Constructors
    Node():index(0),First(0),Second(0),name1("x"),name2("y"){}
    Node(double in1,double in2):index(0),First(in1),Second(in2),name1("x"),name2("y"){}
    Node(std::string n1, double d1, std::string n2, double d2):index(0),First(d1),Second(d2),name1(n1),name2(n2){}
    Node(std::string n1, std::string n2):index(0),First(0),Second(0),name1(n1),name2(n2){}
    Node(int i, double in1,double in2):index(i),First(in1),Second(in2),name1("x"),name2("y"){}
    Node(int i, std::string n1, double d1, std::string n2, double d2):index(i),First(d1),Second(d2),name1(n1),name2(n2){}
    Node(int i, std::string n1, std::string n2):index(0),First(0),Second(i),name1(n1),name2(n2){}
    

//Destructor
//Operators
    double operator()(int index){
        if(index==0){return First;}
        else if(index==1){return Second;}
        else{std::cout<<"Node index out of range\n";return 0.0/0.0;}
    }
    double& operator[](int index){
        //double temp(0.0/0.0);
        if(index==0){return First;}
        else if(index==1){return Second;}
        else{std::cout<<"Node index out of range\n";return First;}
    }
//Methods
    std::string info(){
        return name1+": "+std::to_string(First)+", "+name2+": "+std::to_string(Second)+"\n";
    }
    std::string eval(){
        return std::to_string(First)+" "+std::to_string(Second);
    }
    std::pair<double,double> pos(){
        return {First,Second};
    }
};

template<int size>
class Element{
public:
//Atttributes
    int associate[size];
    double a[size], b[size], c[size], Area;
    Node midpoint;
//Constructor
    Element():associate(),a(),b(),c(),Area(0),midpoint(){}
//Operators
//Methods
    template<int i>
    double phi(std::pair<double,double> position){
        return a[i] + b[i] * position.first + c[i] * position.second;
    }

    std::string info(){
        std::string output;
        output= "associated nodes: {"+std::to_string(associate[0]);
        for(int i=1;i<size;i++){
            output+= ","+std::to_string(associate[i]);
        }
        output+="}\nTest functions:\n";
        for(int i=0;i<size;i++){
            output+="TF"+std::to_string(i)+": "+std::to_string(a[i])+" + "+std::to_string(b[i])+"*"+midpoint.name1+" + "+std::to_string(c[i])+"*"+midpoint.name2+"\n";
        }
        output+="Area: "+std::to_string(Area)+"\n";
        output+="midpoint:\n"+midpoint.info();
        return output;
    }

    std::string eval(){
        std::string output;
        output += midpoint.eval();
        return output;
    }

    std::string shorthand(){
        return std::to_string(Area)+"\n";
    }
};

class Triangle : public Element<3>{
public:
//Attributes
//Constructors
    Triangle():Element<3>(){}
    Triangle(int n0, int n1, int n2, Node N0, Node N1, Node N2){
        associate[0] = n0;associate[1] = n1;associate[2] = n2;
        double D = (N1(0)-N0(0))*(N2(1)-N0(1)) - (N1(1)-N0(1))*(N2(0)-N0(0));
        Area = fabs(D/2.0);
        double temp = 1/D;
        b[0] = temp*(N1(1)-N2(1)); b[1] = temp*(N2(1)-N0(1)); b[2] = temp*(N0(1) - N1(1));
        c[0] = temp*(N2(0)-N1(0)); c[1] = temp*(N0(0)-N2(0)); c[2] = temp*(N1(0) - N0(0));
        a[0] = 1 - b[0]*N0(0) - c[0]*N0(1);
        a[1] = 1 - b[1]*N1(0) - c[1]*N1(1);
        a[2] = 1 - b[2]*N2(0) - c[2]*N2(1);
        midpoint = Node((N0(0)+N1(0)+N2(0))/3.0,(N0(1)+N1(1)+N2(1))/3);
    }
    Triangle(Node N0, Node N1, Node N2):Triangle(N0.index,N1.index,N2.index,N0,N1,N2){}
    Triangle(int n[3], Node N[3]):Triangle(n[0],n[1],n[2],N[0],N[1],N[2]){}
    Triangle(Node N[3]):Triangle(N[0].index, N[1].index, N[2].index, N[0],N[1],N[2]){}
//Operators
//Methods
};

//Static memory space discretization
template<int n,int e, int etype>
class Discretization{
public:
//Attributes
    Node nodes[n];
    Element<etype> elements[e];
    Eigen::VectorXi external, internal;
//Constructors
    Discretization():nodes(),elements(),external(),internal(){}
    Discretization(int ex, int in):nodes(),elements(),external(ex),internal(in){}
    Discretization(Node innodes[n], Element<etype> inelements[e]):Discretization(){
        for(int i=0;i<n;i++){
            nodes[i] = innodes[i];
        }
        for(int i=0;i<e;i++){
            elements[i] = inelements[i];
        }
    }
//Operators
//Methods
    std::string info(){
        std::string output(std::to_string(n)+" Nodes\n");
        for(int i=0;i<n;i++){
            output+=nodes[i].info();
        }
        return output;
    }
    std::string evalnodes(){
        std::string output;
        for(int i=0;i<n;i++){
            output+=nodes[i].eval()+"\n";
        }
        return output;
    }

    std::string evalelements(){
        std::string output;
        for(int i=0;i<e;i++){
            output+=elements[i].eval()+"\n";
        }
        return output;
    }

    std::string shorthand(){
        std::string output;
        for(int i=0;i<e;i++){
            output+=elements[i].shorthand();
        }
        return output;
    }
};

//Flat surface with basic boundaries
template<int dim1, int dim2>
class TriangulatedGrid : public Discretization<dim1*dim2,2*(dim1-1)*(dim2-1),3>{
public:
//Attributes
//Constructors
    TriangulatedGrid(){}
    TriangulatedGrid(double Lowerdim1, double Upperdim1, double Lowerdim2, double Upperdim2):
    Discretization<dim1*dim2,2*(dim1-1)*(dim2-1),3>(2*(dim1+dim2-2),(dim1-2)*(dim2-2))
    {
        double ddim1((Upperdim1-Lowerdim1)/((double)(dim1)-1.5)), ddim2((Upperdim2-Lowerdim2)/((double)(dim2-1)));
        bool switchstate(true);
        int offset,offset1;
        for(int j=0;j<dim2;j++){
            offset = j*dim1;
            if(switchstate){
                this->nodes[offset] = Node(offset,Lowerdim1,j*ddim2+Lowerdim2);
                this->nodes[offset+1] = Node(offset+1,Lowerdim1+(ddim1/2.0),j*ddim2+Lowerdim2);
                for(int i=2;i<dim1;i++){
                    this->nodes[offset+i] = Node(offset+i,((double)i-0.5)*ddim1+Lowerdim1,j*ddim2+Lowerdim2);
                }
                switchstate=false;
            }
            else{
                for(int i=0;i<dim1-1;i++){
                    this->nodes[offset+i] = Node(offset+i,i*ddim1+Lowerdim1,j*ddim2+Lowerdim2);
                }
                this->nodes[(j+1)*dim1-1] = Node((j+1)*dim1-1,Upperdim1, j*ddim2+Lowerdim2);              
                switchstate=true;
            }
        }
        switchstate=true;
        for(int j=0;j<dim2-1;j++){
            offset=j*(dim1-1);
            //Could me replaced by offset + j but makes filling in unreadable.
            //This is why nodes store their own integer, passing the index is also allowed but makes the code unreadable again
            if(switchstate){
                offset1=j*dim1;
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+i] = Triangle(this->nodes[offset1+i],this->nodes[offset1+i+1],this->nodes[offset1+i+dim1]);
                    //print(2*offset+i);
                }
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+dim1-1+i] = Triangle(this->nodes[offset1+dim1+i],this->nodes[offset1+dim1+i+1],this->nodes[offset1+i+1]);
                    //print(2*offset+dim1-1+i);
                }
                switchstate=false;
            }
            else{
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+i] = Triangle(this->nodes[offset1+dim1+i],this->nodes[offset1+dim1+i+1],this->nodes[offset1+2*dim1+i+1]);
                    //print(2*offset+i);
                }
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+dim1-1+i] = Triangle(this->nodes[offset1+2*dim1+i],this->nodes[offset1+2*dim1+i+1],this->nodes[offset1+dim1+i]);
                    //print(2*offset+dim1-1+i);
                }
                switchstate=true;
            }
        }
        for(int i=0;i<dim1;i++){
            this->external(i)=i;
        }
        for(int i=0;i<dim2-2;i++){
            this->external(dim1+2*i)=dim1+i*dim1;
            this->external(dim1+2*i+1)=2*dim1+i*dim1-1;
        }
        for(int i=0;i<dim1;i++){
            this->external(dim1+2*(dim2-2)+i)=dim1*(dim2-1)+i;
        }
        int ex(0);
        int count(0);
        bool safe=true;
        for(int i=0;i<dim1*dim2;i++){
            if(safe){
                if(i==this->external(ex)){
                    ex++;
                    if(ex==this->external.size()){
                        safe=false;
                    }
                }
                else{
                    this->internal(count)=i;
                    count++;
                }
            }
            else{
                this->internal(count)=i;
                count++;
            }
        }
    }
//Operators
//Methods
};
template<int dim1, int dim2>
class TriangulatedCylinder : public Discretization<dim1*dim2, 2*dim1*(dim2-1), 3>{
public:
    //Attributes
    //Constructors
    TriangulatedCylinder(){}
    TriangulatedCylinder(double LowerZ, double UpperZ)
    :Discretization<dim1*dim2, 2*dim1*(dim2-1), 3>(2*dim1,(dim2-2)*dim1)
    {
        double ddim1(2*M_PI/((double)dim1)), ddim2((UpperZ-LowerZ)/((double)(dim2-1))),
        tempdim1, tempdim2(LowerZ), smoll1(ddim1/2.0);
        bool switchstate(true);
        int offset;
        for(int j=0;j<dim2;j++){
            offset=j*dim1;
            if(switchstate){
                tempdim1=0;
                switchstate=false;
            }
            else{
                tempdim1 = ddim1/2.0;
                switchstate=true;
            }
            for(int i=0;i<dim1;i++){
                this->nodes[offset+i] = Node(offset+i,"theta",tempdim1,"z",tempdim2);
                tempdim1+=ddim1;
            }
            tempdim2+=ddim2;
        }
        switchstate=true;
        for(int j=0;j<dim2-1;j++){
            offset = j*dim1;
            if(switchstate){
                //print("case 0");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+2*dim1-1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1],Node(offset,"theta", 2*M_PI,"z",this->nodes[offset].Second),this->nodes[offset+2*dim1-1]);
                //print("case 1");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+1+i<<","<<offset+1+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+1+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1],Node(offset+dim1,"theta", 2*M_PI + smoll1, "z", this->nodes[offset+dim1].Second),Node(offset, "theta", 2*M_PI, "z", this->nodes[offset].Second));
                switchstate=false;
            }
            else{
                //print("case 2");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+dim1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1], Node(offset, "theta", 2*M_PI+smoll1, "z", this->nodes[offset].Second), Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second));
                //print("case 3");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+i+1<<","<<offset+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset+dim1-1<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1], Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second), this->nodes[offset+dim1-1]);
                switchstate=true;
            }
        }
        for(int i=0;i<dim1;i++){
            this->external(i)=i;
            this->external(2*dim1-1-i)=dim1*dim2-1-i;
        }
        int ex(0);
        int count(0);
        bool safe=true;
        for(int i=0;i<dim1*dim2;i++){
            if(safe){
                if(i==this->external(ex)){
                    ex++;
                    if(ex==this->external.size()){
                        safe=false;
                    }
                }
                else{
                    this->internal(count)=i;
                    count++;
                }
            }
            else{
                this->internal(count)=i;
                count++;
            }
        }
        //std::cout<<(this->externalE.size())<<"\n";
        /*
        for(int i=0;i<dim1;i++){
            this->externalE(i)=i;
            this->externalE(2*dim1-1-i)=2*dim1*(dim2-1)-1-i;
        }
        //std::cout<<this->externalE<<"\n";
        //std::cout<<"\n";
        //std::cout<<(this->internalE.size())<<"\n";
        ex=0;count=0;safe=true;
        for(int i=0;i<2*dim1*(dim2-1);i++){
            if(safe){
                if(i==this->externalE(ex)){
                    ex++;
                    if(ex==this->externalE.size()){
                        safe=false;
                    }
                }
                else{
                    this->internalE(count)=i;
                    count++;
                }
            }
            else{
                this->internalE(count)=i;
                count++;
            }
        }
        */
        //std::cout<<this->internalE<<"\n";
    }

};



//Dynamic memory size discretization
template<int etype>
class DynDiscretization{
public:
//Attributes
    int n,e;
    Node* nodes;
    Element<etype>* elements;
    Eigen::VectorXi external, internal;
//Constructors
    DynDiscretization():n(),e(),nodes(),elements(),external(),internal(){}
    DynDiscretization(int N,int E):n(N),e(E),nodes(new Node[n]),elements(new Element<etype>[e]),external(),internal(){}
    DynDiscretization(int N,int E,int ex, int in):n(N),e(E),nodes(new Node[n]),elements(new Element<etype>[e]),external(ex),internal(in){}
//Operators
//Methods
    std::string info(){
        std::string output(std::to_string(n)+" Nodes\n");
        for(int i=0;i<n;i++){
            output+=nodes[i].info();
        }
        return output;
    }
    std::string evalnodes(){
        std::string output;
        for(int i=0;i<n;i++){
            output+=nodes[i].eval()+"\n";
        }
        return output;
    }

    std::string evalelements(){
        std::string output;
        for(int i=0;i<e;i++){
            output+=elements[i].eval()+"\n";
        }
        return output;
    }

    std::string shorthand(){
        std::string output;
        for(int i=0;i<e;i++){
            output+=elements[i].shorthand();
        }
        return output;
    }
};

class DynTriangulatedGrid : public DynDiscretization<3>{
public:
//Attributes
    int dim1,dim2;
//Constructors
    DynTriangulatedGrid():dim1(),dim2(){}
    DynTriangulatedGrid(int d1,int d2):DynDiscretization<3>(d1*d2,2*(d1-1)*(d2-1),2*(d1+d2-2),(d1-2)*(d2-2)),
    dim1(d1),dim2(d2){}
    
    DynTriangulatedGrid(int d1,int d2,double Lowerdim1, double Upperdim1, double Lowerdim2, double Upperdim2):
    DynTriangulatedGrid(d1,d2)
    {
        double ddim1((Upperdim1-Lowerdim1)/((double)(dim1)-1.5)), ddim2((Upperdim2-Lowerdim2)/((double)(dim2-1)));
        bool switchstate(true);
        int offset,offset1;
        for(int j=0;j<dim2;j++){
            offset = j*dim1;
            if(switchstate){
                this->nodes[offset] = Node(offset,Lowerdim1,j*ddim2+Lowerdim2);
                this->nodes[offset+1] = Node(offset+1,Lowerdim1+(ddim1/2.0),j*ddim2+Lowerdim2);
                for(int i=2;i<dim1;i++){
                    this->nodes[offset+i] = Node(offset+i,((double)i-0.5)*ddim1+Lowerdim1,j*ddim2+Lowerdim2);
                }
                switchstate=false;
            }
            else{
                for(int i=0;i<dim1-1;i++){
                    this->nodes[offset+i] = Node(offset+i,i*ddim1+Lowerdim1,j*ddim2+Lowerdim2);
                }
                this->nodes[(j+1)*dim1-1] = Node((j+1)*dim1-1,Upperdim1, j*ddim2+Lowerdim2);              
                switchstate=true;
            }
        }
        switchstate=true;
        for(int j=0;j<dim2-1;j++){
            offset=j*(dim1-1);
            //Could me replaced by offset + j but makes filling in unreadable.
            //This is why nodes store their own integer, passing the index is also allowed but makes the code unreadable again
            if(switchstate){
                offset1=j*dim1;
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+i] = Triangle(this->nodes[offset1+i],this->nodes[offset1+i+1],this->nodes[offset1+i+dim1]);
                    //print(2*offset+i);
                }
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+dim1-1+i] = Triangle(this->nodes[offset1+dim1+i],this->nodes[offset1+dim1+i+1],this->nodes[offset1+i+1]);
                    //print(2*offset+dim1-1+i);
                }
                switchstate=false;
            }
            else{
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+i] = Triangle(this->nodes[offset1+dim1+i],this->nodes[offset1+dim1+i+1],this->nodes[offset1+2*dim1+i+1]);
                    //print(2*offset+i);
                }
                for(int i=0;i<dim1-1;i++){
                    this->elements[2*offset+dim1-1+i] = Triangle(this->nodes[offset1+2*dim1+i],this->nodes[offset1+2*dim1+i+1],this->nodes[offset1+dim1+i]);
                    //print(2*offset+dim1-1+i);
                }
                switchstate=true;
            }
        }
        for(int i=0;i<dim1;i++){
            this->external(i)=i;
        }
        for(int i=0;i<dim2-2;i++){
            this->external(dim1+2*i)=dim1+i*dim1;
            this->external(dim1+2*i+1)=2*dim1+i*dim1-1;
        }
        for(int i=0;i<dim1;i++){
            this->external(dim1+2*(dim2-2)+i)=dim1*(dim2-1)+i;
        }
        int ex(0);
        int count(0);
        bool safe=true;
        for(int i=0;i<dim1*dim2;i++){
            if(safe){
                if(i==this->external(ex)){
                    ex++;
                    if(ex==this->external.size()){
                        safe=false;
                    }
                }
                else{
                    this->internal(count)=i;
                    count++;
                }
            }
            else{
                this->internal(count)=i;
                count++;
            }
        }
    }
    
//Operators
//Methods
};

class DynTriangulatedCylinder : public DynDiscretization<3>{
public:
    //Attributes
    //Constructors
    DynTriangulatedCylinder(){}
    DynTriangulatedCylinder(int dim1, int dim2):
    DynDiscretization<3>(dim1*dim2, 2*dim1*(dim2-1),2*dim1,(dim2-2)*dim1){}
    DynTriangulatedCylinder(int dim1, int dim2, double LowerZ, double UpperZ):DynTriangulatedCylinder(dim1,dim2)
    {
        double ddim1(2*M_PI/((double)dim1)), ddim2((UpperZ-LowerZ)/((double)(dim2-1))),
        tempdim1, tempdim2(LowerZ), smoll1(ddim1/2.0);
        bool switchstate(true);
        int offset;
        for(int j=0;j<dim2;j++){
            offset=j*dim1;
            if(switchstate){
                tempdim1=0;
                switchstate=false;
            }
            else{
                tempdim1 = ddim1/2.0;
                switchstate=true;
            }
            for(int i=0;i<dim1;i++){
                this->nodes[offset+i] = Node(offset+i,"theta",tempdim1,"z",tempdim2);
                tempdim1+=ddim1;
            }
            tempdim2+=ddim2;
        }
        switchstate=true;
        for(int j=0;j<dim2-1;j++){
            offset = j*dim1;
            if(switchstate){
                //print("case 0");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+2*dim1-1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1],Node(offset,"theta", 2*M_PI,"z",this->nodes[offset].Second),this->nodes[offset+2*dim1-1]);
                //print("case 1");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+1+i<<","<<offset+1+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+1+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1],Node(offset+dim1,"theta", 2*M_PI + smoll1, "z", this->nodes[offset+dim1].Second),Node(offset, "theta", 2*M_PI, "z", this->nodes[offset].Second));
                switchstate=false;
            }
            else{
                //print("case 2");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+dim1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1], Node(offset, "theta", 2*M_PI+smoll1, "z", this->nodes[offset].Second), Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second));
                //print("case 3");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+i+1<<","<<offset+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset+dim1-1<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1], Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second), this->nodes[offset+dim1-1]);
                switchstate=true;
            }
        }
        for(int i=0;i<dim1;i++){
            this->external(i)=i;
            this->external(2*dim1-1-i)=dim1*dim2-1-i;
        }
        int ex(0);
        int count(0);
        bool safe=true;
        for(int i=0;i<dim1*dim2;i++){
            if(safe){
                if(i==this->external(ex)){
                    ex++;
                    if(ex==this->external.size()){
                        safe=false;
                    }
                }
                else{
                    this->internal(count)=i;
                    count++;
                }
            }
            else{
                this->internal(count)=i;
                count++;
            }
        }
        
    }

};


class DynTriangulatedDoughnut : public DynDiscretization<3>{
public:
    DynTriangulatedDoughnut(){}
    DynTriangulatedDoughnut(int dim1, int dim2):
    DynDiscretization<3>(dim1*dim2, 2*dim1*dim2, 0, dim1*dim2){}
    DynTriangulatedDoughnut(int dim1, int dim2, double LowerZ, double UpperZ):
    DynTriangulatedDoughnut(dim1,dim2)
    {
        double ddim1(2*M_PI/((double)dim1)), ddim2((UpperZ-LowerZ)/((double)(dim2))),
        tempdim1, tempdim2(LowerZ), smoll1(ddim1/2.0);
        bool switchstate(true);
        int offset;
        //std::cout<<"test3\n";
        for(int j=0;j<dim2;j++){
            offset=j*dim1;
            if(switchstate){
                tempdim1=0;
                switchstate=false;
            }
            else{
                tempdim1 = ddim1/2.0;
                switchstate=true;
            }
            for(int i=0;i<dim1;i++){
                //std::cout<<offset+i<<"\n";
                this->nodes[offset+i] = Node(offset+i,"theta",tempdim1,"z",tempdim2);
                tempdim1+=ddim1;
            }
            tempdim2+=ddim2;
        }
        switchstate=true;
        for(int j=0;j<dim2-1;j++){
            offset = j*dim1;
            if(switchstate){
                //print("case 0");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+2*dim1-1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1],Node(offset,"theta", 2*M_PI,"z",this->nodes[offset].Second),this->nodes[offset+2*dim1-1]);
                //print("case 1");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+1+i<<","<<offset+1+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+1+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1],Node(offset+dim1,"theta", 2*M_PI + smoll1, "z", this->nodes[offset+dim1].Second),Node(offset, "theta", 2*M_PI, "z", this->nodes[offset].Second));
                switchstate=false;
            }
            else{
                //print("case 2");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+i);
                    //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+1+i<<std::endl;
                    this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+1+i]);
                }
                //print(2*offset+dim1-1);
                //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+dim1<<std::endl;
                this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1], Node(offset, "theta", 2*M_PI+smoll1, "z", this->nodes[offset].Second), Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second));
                //print("case 3");
                for(int i=0;i<dim1-1;i++){
                    //print(2*offset+dim1+i);
                    //std::cout<<offset+dim1+i<<","<<offset+dim1+i+1<<","<<offset+i<<std::endl;
                    this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+i]);
                }
                //print(2*(offset+dim1)-1);
                //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset+dim1-1<<std::endl;
                this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1], Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second), this->nodes[offset+dim1-1]);
                switchstate=true;
            }
        }
        offset = (dim2-1)*(dim1);
        //print("case 2");
        for(int i=0;i<dim1-1;i++){
            //print(2*offset+i);
            //std::cout<<offset+i<<","<<offset+1+i<<","<<offset+dim1+1+i<<std::endl;
            //this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],this->nodes[offset+dim1+1+i]);
            this->elements[2*offset+i] = Triangle(this->nodes[offset+i],this->nodes[offset+1+i],Node(i+1, "theta", this->nodes[i+1].First , "z", UpperZ));
        }
        //print(2*offset+dim1-1);
        //std::cout<<offset+dim1-1<<","<<offset<<","<<offset+dim1<<std::endl;
        //this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1], Node(offset, "theta", 2*M_PI+smoll1, "z", this->nodes[offset].Second), Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second));
        this->elements[2*offset+dim1-1] = Triangle(this->nodes[offset+dim1-1], Node(offset, "theta", 2*M_PI+smoll1, "z", this->nodes[offset].Second), Node(0, "theta", 2*M_PI, "z", UpperZ));

        //print("case 3");
        for(int i=0;i<dim1-1;i++){
            //print(2*offset+dim1+i);
            //std::cout<<offset+dim1+i<<","<<offset+dim1+i+1<<","<<offset+i<<std::endl;
            //this->elements[2*offset+dim1+i] = Triangle(this->nodes[offset+dim1+i],this->nodes[offset+dim1+1+i],this->nodes[offset+i]);
            this->elements[2*offset+dim1+i] = Triangle(Node(i, "theta", this->nodes[i].First, "z", UpperZ),Node(i+1, "theta", this->nodes[i+1].First , "z", UpperZ),this->nodes[offset+i]);
        }
        //print(2*(offset+dim1)-1);
        //std::cout<<offset+2*dim1-1<<","<<offset+dim1<<","<<offset+dim1-1<<std::endl;
        //this->elements[2*(offset+dim1)-1] = Triangle(this->nodes[offset+2*dim1-1], Node(offset+dim1, "theta", 2*M_PI, "z", this->nodes[offset+dim1].Second), this->nodes[offset+dim1-1]);
        this->elements[2*(offset+dim1)-1] = Triangle(Node(dim1-1, "theta", this->nodes[dim1-1].First, "z", UpperZ), Node(0, "theta", 2*M_PI, "z", UpperZ), this->nodes[offset+dim1-1]);
        for(int i=0;i<dim1*dim2;i++){
            internal(i) = i;
        }
    }

};





#endif //DISCRETIZATION_HPP