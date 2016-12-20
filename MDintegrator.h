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
    MDsystem* MD;
    double deltaT;

    virtual void integrate() = 0;

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
