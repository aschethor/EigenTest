//
// Created by Ni on 16.10.2016.
//

#include "Argon.h"
#include <iostream>
#include <ctime>

#define SQR(x) (x*x)
#define CUBE(x) (x*x*x)

using namespace Eigen;
using namespace std;

Argon::Argon(int n,Eigen::Vector3d Box,double Rcutoff,BoundaryCondition bc){
    this->n = n;
    this->Box = Box;
    this->Rcutoff = Rcutoff;
    Ecutoff = 4.0*(1.0/CUBE(SQR(SQR(Rcutoff)))-1.0/CUBE(SQR(Rcutoff)));
    this->bc = bc;
    Force = MatrixXd::Zero(3,n);
    m = VectorXd(n);
    m.setConstant(mass);
    srand((unsigned int) time(0));
    x = Box.asDiagonal()*(MatrixXd::Random(3,n)+MatrixXd::Ones(3,n))/2;
    v = MatrixXd::Random(3,n);
    zeroMomentum();
}

/**
 * Initializes fcc - lattice of Argon atoms.
 * @param nx        lattice size in x direction
 * @param ny        lattice size in y direction
 * @param nz        lattice size in z direction
 * @param a         lattice spacing
 * @param vmax      maximum velocity of randomly initialized atom-velocities
 * @param Rcutoff   cutoff-radius of lennard jones potential
 * @param bc        boundary condition
 * @return
 */
Argon::Argon(int nx,int ny,int nz,double a,double vmax,double Rcutoff,BoundaryCondition bc){
    n = 4*nx*ny*nz;
    this->Box = Vector3d(nx*a,ny*a,nz*a);
    this->Rcutoff = Rcutoff;
    Ecutoff = 4.0*(1.0/CUBE(SQR(SQR(Rcutoff)))-1.0/CUBE(SQR(Rcutoff)));
    this->bc = bc;
    Force = MatrixXd::Zero(3,n);
    m = VectorXd(n);
    m.setConstant(mass);
    x = MatrixXd::Zero(3,n);
    int i=0;
    for(int x = 0;x<nx;x++)
        for(int y = 0;y<ny;y++)
            for(int z=0;z<nz;z++){
                this->x(0, i) = x * a;
                this->x(1, i) = y * a;
                this->x(2, i) = z * a;
                i++;
                this->x(0, i) = (x+0.5) * a;
                this->x(1, i) = (y+0.5) * a;
                this->x(2, i) = z * a;
                i++;
                this->x(0, i) = (x+0.5) * a;
                this->x(1, i) = y * a;
                this->x(2, i) = (z+0.5) * a;
                i++;
                this->x(0, i) = x * a;
                this->x(1, i) = (y+0.5) * a;
                this->x(2, i) = (z+0.5) * a;
                i++;
            }
    x = x+MatrixXd::Ones(3,n)*a/2;
    v = MatrixXd::Random(3,n)*vmax;
    zeroMomentum();
}

void Argon::evaluate() {
    Energy_pot = 0;
    //Force = MatrixXd::Zero(3,n);
    for(int i=0;i<n;i++)frc(i,0)=frc(i,1)=frc(i,2)=0;
    Pressure = 0;
    double drx,dry,drz;
    double BoxX = Box(0),BoxY=Box(1),BoxZ=Box(2);
    double Ff,r2i,r6i;

    for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
            drx = pos(i,0)-pos(j,0);
            dry = pos(i,1)-pos(j,1);
            drz = pos(i,2)-pos(j,2);


            drx-=BoxX*rint(drx/BoxX);
            dry-=BoxY*rint(dry/BoxY);
            drz-=BoxZ*rint(drz/BoxZ);

            r2i = SQR(drx)+SQR(dry)+SQR(drz);
            if(r2i<SQR(Rcutoff)){
                r2i = 1.0/r2i;
                r6i = CUBE(r2i);

                Energy_pot+=4.0*r6i*(r6i-1)-Ecutoff;
                Ff = 48.0*r6i*(r6i-0.5);
                Pressure+=Ff;
                Ff*=r2i;

                frc(i,0)+=Ff*drx;
                frc(i,1)+=Ff*dry;
                frc(i,2)+=Ff*drz;

                frc(j,0)-=Ff*drx;
                frc(j,1)-=Ff*dry;
                frc(j,2)-=Ff*drz;
            }
        }
    }
    Energy_kin = 0.5*(v.array()*v.array()).colwise().sum().matrix()*m.array().inverse().matrix();
    Pressure/=3.0*BoxX*BoxY*BoxZ;
}

double Argon::Epot() {
    return Energy_pot;
}

double Argon::Ekin() {
    return Energy_kin;
}

Eigen::MatrixXd Argon::F() {
    return Force;
}

double Argon::P(){
    return Pressure;
}

double Argon::T(){
    return 2.0*Energy_kin/(3*n-3);
}

string Argon::toPDB(){
    char str[100];
    string ret="";
    for(int i=0;i<n;i++){
        sprintf(str,"%s%7d%s%12d%s%8.3lf%8.3lf%8.3lf\n",
                "ATOM",1," Ar",1," xx ",pos(i,0),pos(i,1),pos(i,2));
        ret.append(str);
    }
    return ret;
}

void Argon::print(){
    MatrixXd output(n,7);
    output<<m,x.transpose(),v.transpose();
    cout<<output.format(IOFormat(5, 0, "\t", "\n"))<<endl;
}

void Argon::zeroMomentum() {
    VectorXd v0 = v.rowwise().sum()/n;
    v = v-v0.asDiagonal()*MatrixXd::Ones(3,n);
}

double& Argon::pos(int n,int i){
    return x.data()[n*3+i];
}

double& Argon::vel(int n,int i){
    return v.data()[n*3+i];
}

double& Argon::frc(int n,int i){
    return Force.data()[n*3+i];
}