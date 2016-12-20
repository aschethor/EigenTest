//
// Created by Ni on 16.10.2016.
//

#include "VelocityVerletIntegrator.h"

using namespace Eigen;

VelocityVerletIntegrator::VelocityVerletIntegrator(MDsystem* MD,double deltaT) {
    this->MD = MD;
    MD->evaluate();
    this->deltaT = deltaT;
}

void VelocityVerletIntegrator::integrate() {
    MD->x = MD->x+MD->v*deltaT+0.5*MD->F()*MD->m.asDiagonal().inverse()*deltaT*deltaT;
    MD->v = MD->v+0.5*MD->F()*MD->m.asDiagonal().inverse()*deltaT;
    MD->evaluate();
    MD->v = MD->v+0.5*MD->F()*MD->m.asDiagonal().inverse()*deltaT;
    MD->backToBox();
}