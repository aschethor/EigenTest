//
// Created by Ni on 16.10.2016.
//

#ifndef EIGENTEST_ARGON_H
#define EIGENTEST_ARGON_H

#include "MDsystem.h"

class Argon : public MDsystem{
    const double sigma = 1;         ///< sigma in LJ-units
    const double mass = 1;          ///< mass in LJ-units
    double Rcutoff;                 ///< cutoff radius
    double Ecutoff;                 ///< Energy of LJ-Potential at cutoff radius
    double Energy_pot;              ///< potential energy
    double Energy_kin;              ///< kinetic energy
    Eigen::MatrixXd Force;          ///< forces on particles
    double Pressure;                ///< pressure

public:
    Argon(int n,Eigen::Vector3d Box,double Rcutoff=3,BoundaryCondition bc=periodic);
    Argon(int nx,int ny,int nz,double a,double vmax=1,double Rcutoff=2.5,BoundaryCondition bc=periodic);

    void evaluate();
    double Epot();
    double Ekin();
    Eigen::MatrixXd F();
    double P();
    double T();
    std::string toPDB();

    void print();

private:
    void zeroMomentum();

public:
    double& pos(int n,int i);
    double& vel(int n,int i);
    double& frc(int n,int i);
};

#endif //EIGENTEST_Argon_H
