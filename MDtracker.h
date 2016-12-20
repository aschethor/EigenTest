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
    MDintegrator* MD;               ///< Moleculardynamic integrator that shall be tracked
    Eigen::MatrixXd* x_track;       ///< positions of particles    (size: 3 x n x nSamples)
    Eigen::MatrixXd* v_track;       ///< velocities of particles   (size: 3 x n x nSamples)
    int nSamples;                   ///< number of Samples
    int SampleRate;                 ///< number of integration steps between Samples
    FILE * pdbFile;                 ///< filepointer for PDB-File
    FILE * systemFile;              ///< filepointer for System-File
    FILE * gofrFile;                ///< filepointer for g(r)-File
    FILE * sofkFile;                ///< filepointer for S(k)-File
    FILE * msDisplFile;             ///< filepointer for mean square displacement File
    FILE * vvAutoFile;              ///< filepointer for velocity-velocity autocorrelation File
    int gofrSize;                   ///< number of gofr-Bins
    std::vector<double> gofr;       ///< g(r): r goes from 0 to min(BoxSize of MDsystem)/2 (vector size: gofrsize)
    int sofkSize;                   ///< number of sofk-Bins
    const double maxK = 30;         ///< maximum value for k in S(k)
    std::vector<double> sofk;       ///< S(k): k goes from 0 to maxK (vector size: sofksize)
    std::vector<double> msDispl;    ///< mean square Displacement (vector size: nSamples - 1)
    std::vector<double> vvAuto;     ///< velocity-velocity auto correlation (vector size: nSamples - 1)

    void calculateGofr();
    void calculateSofk();
    void calculateMsDispl();
    void calculateVvAuto();

public:
    MDtracker(std::string resultFolder,MDintegrator* MD,int nSamples,int SampleRate=10,std::string suffix = "");
    ~MDtracker();
    void track();
};


#endif //EIGENTEST_MDTRACKER_H
