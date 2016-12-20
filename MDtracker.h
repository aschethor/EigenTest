//
// Created by Ni on 16.10.2016.
//

#ifndef EIGENTEST_MDTRACKER_H
#define EIGENTEST_MDTRACKER_H

#include "MDintegrator.h"
#include "Histogram.h"
#include <vector>

class MDtracker {
private:
    MDintegrator* MD;
    Eigen::MatrixXd* x_track;       //positions of particles    (size: 3 x n x nSamples)
    Eigen::MatrixXd* v_track;       //velocities of particles   (size: 3 x n x nSamples)
    int nSamples;                   //number of Samples
    int SampleRate;                 //number of integration steps between Samples
    FILE * pdbFile;
    FILE * systemFile;
    FILE * gofrFile;
    FILE * sofkFile;
    FILE * msDisplFile;
    FILE * vvAutoFile;
    std::vector<double> gofr;       //gofr: r goes from 0 to min(BoxSize of MDsystem)/2
    int gofrSize;                   //number of gofr-Bins
    const double maxS = 30;
    std::vector<double> sofk;       //sofk: k goes from 0 to maxS
    int sofkSize;                   //number of sofk-Bins
    std::vector<double> msDispl;    //mean square Displacement (size: nSamples - 1)
    std::vector<double> vvAuto;     //velocity-velocity auto correlation (size: nSamples - 1)

    void calculateGofr();
    void calculateSofk();
    void calculateMsDispl();
    void calculateVvAuto();

    std::vector<double> E_pot;
    std::vector<double> E_kin;
    std::vector<double> E_tot;
    std::vector<double> P;
    std::vector<double> T;

public:
    MDtracker(std::string resultFolder,MDintegrator* MD,int nSamples,int SampleRate=10,std::string suffix = "");
    ~MDtracker();
    void track();
};


#endif //EIGENTEST_MDTRACKER_H
