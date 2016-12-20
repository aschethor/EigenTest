//
// Created by Ni on 16.10.2016.
//

#ifndef EIGENTEST_MDSYSTEM_H
#define EIGENTEST_MDSYSTEM_H

#include <Eigen/Dense>
#include <stdio.h>

//TODO: add Time

class MDsystem {
public:
    enum BoundaryCondition{
        periodic
        //, mirror          //evtl. later...
    };
public:
    int n;                              ///< number of particles
    Eigen::VectorXd m;                  ///< masses of particles       (size: n)
    Eigen::MatrixXd x;                  ///< positions of particles    (size: 3 x n)
    Eigen::MatrixXd v;                  ///< velocities of particles   (size: 3 x n)
    Eigen::Vector3d Box;                ///< size of Box in x/y/z direction
    BoundaryCondition bc;               ///< boundary condition
    virtual void evaluate() = 0;        ///< evaluate E / F / P
    virtual double Epot() = 0;          ///< potential energy
    virtual double Ekin() = 0;          ///< kinetic energy
    virtual Eigen::MatrixXd F() = 0;    ///< forces
    virtual double P() = 0;             ///< pressure
    virtual double T() = 0;             ///< temperature
    virtual std::string toPDB() = 0;    ///< PDB-representation of the system

    /**
     * put particles back into box
     */
    void backToBox(){
        if(bc==periodic) {
            double BoxX = Box(0), BoxY = Box(1), BoxZ = Box(2);
            for (int i = 0; i < n; i++) {
                x.data()[i * 3 + 0] = fmod(x.data()[i * 3 + 0],BoxX);
                x.data()[i * 3 + 1] = fmod(x.data()[i * 3 + 1],BoxY);
                x.data()[i * 3 + 2] = fmod(x.data()[i * 3 + 2],BoxZ);
                if(x.data()[i*3+0]<0)x.data()[i*3+0]+=BoxX;
                if(x.data()[i*3+1]<0)x.data()[i*3+1]+=BoxY;
                if(x.data()[i*3+2]<0)x.data()[i*3+2]+=BoxZ;
            }
        }
    }

    double Etot(){
        return Epot()+Ekin();
    }

};

#endif //EIGENTEST_MDSYSTEM_H
