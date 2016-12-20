//
// Created by Ni on 16.10.2016.
//

#ifndef EIGENTEST_VELOCITYVERLETINTEGRATOR_H
#define EIGENTEST_VELOCITYVERLETINTEGRATOR_H


#include "MDintegrator.h"

class VelocityVerletIntegrator : public MDintegrator{
public:
    VelocityVerletIntegrator(MDsystem* MD,double deltaT = 0.001);
    void integrate();
};

#endif //EIGENTEST_VELOCITYVERLETINTEGRATOR_H
