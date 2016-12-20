//
// Created by Ni on 16.10.2016.
//

#ifndef EIGENTEST_MDINTEGRATOR_H
#define EIGENTEST_MDINTEGRATOR_H

#include "MDsystem.h"
#include <iostream>
#include <cmath>

class MDintegrator {
public:
    MDsystem* MD;   ///< Molecular dynamic system that shall be integrated
    double deltaT;  ///< timestep for integration

    virtual void integrate() = 0;

    /**
     * equilibrate system
     * @param nSteps number of equilibration steps
     * @param T target temperature of system
     */
    void equilibrate(int nSteps,double T){
        std::cout<<"equilibrating...\n";
        for(int i=0;i<nSteps;i++){
            integrate();
            MD->v = MD->v*sqrt(T/MD->T());
        }
        std::cout<<"done.\n";
    }
};


#endif //EIGENTEST_MDINTEGRATOR_H
