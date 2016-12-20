//
// Created by Ni on 18.10.2016.
//

#ifndef EIGENTEST_HISTOGRAM_H
#define EIGENTEST_HISTOGRAM_H
#include <vector>

class Histogram {
private:
    int nBins;
    double binSize;
    double binStart;
    int nSamples = 0;
    int toLowSamples = 0;
    int toHighSamples = 0;
    std::vector<int> bins;

public:
    Histogram(int n,double binSize,double binStart=0);
    int getSize(){return nBins;};
    void addSample(int pos);
    void addSample(double pos);
    int getBin(int n);
    double getNormalizedBin(int n);
    double getPos(int n);
    double getBinSize(){ return binSize;}
    int getToHigh();
    int getToLow();
    std::vector<double> getNormalizedBins();
    Histogram& operator+=(Histogram& h);
};


#endif //EIGENTEST_HISTOGRAM_H
