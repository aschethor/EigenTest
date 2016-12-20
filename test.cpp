#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include "Argon.h"
#include "VelocityVerletIntegrator.h"
#include "MDtracker.h"
#include <ctime>

using namespace std;
using namespace Eigen;

void test1(){
    Vector3f u,v;
    IOFormat CleanFmt(5, 0, ", ", "\n", "[", "]");
    u<<1,0,0;
    v<<0,1,0;
    cout<<u.transpose()<<endl;
    MatrixXf M(4,3);
    M<<1,2,3,
            4,5,6,
            7,8,9,
            1,2,3;
    cout<<M<<endl<<M*u<<endl;
    VectorXf w(4);
    w<<0.5,0.25,0.5,0.312354634;
    cout<<(w.asDiagonal()*M).format(CleanFmt)<<endl;
    cout<<w.asDiagonal()*M<<endl;
    cout<<"Hello World!"<<endl;
}

void test2(){
    MatrixXd M(4,3);
    M<<1,2,3,
            4,5,6,
            7,8,9,
            1,2,3;
    VectorXd w(4);
    w<<0.5,0.25,0.5,0.312354634;
    cout<<w.asDiagonal()*M<<endl;
    Vector3d u;
    u<<1,2,3;
    cout<<u;
}

void test3(int a=5){
    /*
    cout<<a<<endl;
    MatrixXd m = MatrixXd::Zero(2,3);
    m(0,1)=1;
    m(1,2)=2;
    cout<<m<<endl;
    for(int i=0;i<m.cols();i++){
        for(int j=0;j<m.rows();j++){
            cout<<m.data()[i*m.rows()+j]<<" ";
        }
        cout<<endl;
    }*/
    cout<<" 0%";
    cout<<"\b\b\b";
    cout<<"99%";
}

void test4(){
    Argon argon(2,2,2,1);
    argon.print();
    argon.evaluate();
    cout<<argon.Ekin()<<endl;
    cout<<argon.T()<<endl;
    VelocityVerletIntegrator VVInt(&argon);
    VVInt.equilibrate(1000,0.5);
    argon.evaluate();
    cout<<argon.T()<<endl;
    cout<<"-----------------------------------"<<endl;
    for(int i=0;i<10000;i++){
        VVInt.integrate();
        cout<<argon.T()<<endl;
    }
    /*
    cout<<MatrixXd::Random(4,4)<<endl;
    VectorXd v(4);
    v.setConstant(0.5);
    v = v*3;
    cout<<v<<endl;
    cout<<MatrixXd::Ones(4,3)*Vector3d(1,2,3).asDiagonal()<<endl;*/
}

void test5(){
    Argon argon(100,Vector3d(5,5,5),2.49);
    argon.evaluate();
    argon.print();
    cout<<"Ekin: "<<argon.Ekin()<<endl;
    cout<<"Epot: "<<argon.Epot()<<endl;
    cout<<"T   : "<<argon.T()<<endl;
    VelocityVerletIntegrator VVInt(&argon);
    cout<<"Integrator generated"<<endl;
    VVInt.equilibrate(1000,0.5);
    argon.evaluate();
    argon.print();
    cout<<"Ekin: "<<argon.Ekin()<<endl;
    cout<<"Epot: "<<argon.Epot()<<endl;
    cout<<"T   : "<<argon.T()<<endl;
    cout<<"-----------------------------------"<<endl;
    for(int i=0;i<10000;i++){
        VVInt.integrate();
        cout<<i<<": "<<argon.T()<<endl;
    };
    argon.evaluate();
    argon.print();
    cout<<"Ekin: "<<argon.Ekin()<<endl;
    cout<<"Epot: "<<argon.Epot()<<endl;
    cout<<"T   : "<<argon.T()<<endl;
}

void test6(){
    Argon argon(100,Vector3d(5,5,5),2.49);
    VelocityVerletIntegrator VVInt(&argon);
    VVInt.equilibrate(1000,0.5);
    argon.evaluate();
    argon.print();
    cout<<"Ekin: "<<argon.Ekin()<<endl;
    cout<<"Epot: "<<argon.Epot()<<endl;
    cout<<"T   : "<<argon.T()<<endl;
    VVInt.integrate();
    argon.evaluate();
    argon.print();
    cout<<"Ekin: "<<argon.Ekin()<<endl;
    cout<<"Epot: "<<argon.Epot()<<endl;
    cout<<"T   : "<<argon.T()<<endl;
}

void test7(){
    Argon argon(8,8,8,1);
    VelocityVerletIntegrator VVI(&argon);
    cout<<"equilibrating...";
    VVI.equilibrate(500,0.5);
    cout<<"done\nintegrating.";
    //FILE* file = fopen("C:\\Users\\Ni\\Desktop\\test2.pdb","w");
    FILE* file = fopen("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\output\\test2.pdb","w");
    for(int i=0;i<1000;i++){
        fprintf(file,"%s\nENDMDL\n",argon.toPDB().c_str());
        for(int j=0;j<10;j++)VVI.integrate();
        cout<<i*10<<endl;
    }
    cout<<"\ndone\n";
    fclose(file);
}

void test8(){
    long start = clock();
    for(int i=0;i<1000000000;i++){
        int a = 0;
        int b = 1;
        a = b;
        //char* blubb = "uhrihfiergi";
        if(i%100000==0)
        cout<<i<<endl;
    }
    cout<<"duration: "<<((double)(clock()-start))/CLOCKS_PER_SEC<<endl;
}

void test9(){
    Argon argon(864,Eigen::Vector3d(10.2,10.2,10.2));
    //Argon argon(108,Eigen::Vector3d(5.1,5.1,5.1));
    VelocityVerletIntegrator VVI(&argon);
    VVI.equilibrate(1000,0.3);
    MDtracker mDtracker("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\outputtest\\",
        &VVI,1000,10,"N864_T0,3_S10,2");//Number Temperature BoxSize
    mDtracker.track();
    cout<<"Temperature: "<<argon.T()<<endl;
}

void test10(){
    vector<double> xSpace;
    vector<complex<double>> kSpace;
    xSpace.resize(100);
    for(int i=0;i<100;i++){
        xSpace[i] = pow(-1,i)*sin(2*M_PI*i/100);
    }
    kSpace.resize(100);
    Eigen::FFT<double> fft;
    fft.fwd(kSpace,xSpace);
    for(int i=0;i<100;i++){
        cout<<i<<": "<<xSpace[i]<<"  "<<abs(kSpace[i])<<endl;
    }
}

void test11(){
    //Argon argon(8,8,8);
    Argon argon(108,Eigen::Vector3d(5.1,5.1,5.1));
    VelocityVerletIntegrator VVI(&argon);
    VVI.equilibrate(500,0.5);
    MDtracker mDtracker("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\outputtest",
                        &VVI,2000*8,10);
    mDtracker.track();
    cout<<"Temperature: "<<argon.T()<<endl;
}