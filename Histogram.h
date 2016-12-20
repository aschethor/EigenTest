//
// Created by Ni on 18.10.2016.
//

#ifndef EIGENTEST_HISTOGRAM_H
#define EIGENTEST_HISTOGRAM_H
#include <vector>

class Histogram {
private:
    int nBins;                  ///< number of Bins
    double binSize;             ///< size of Bins
    double binStart;            ///< position of first bin
    int nSamples = 0;           ///< number of registered samples
    int toLowSamples = 0;       ///< number of registered samples that were to low for the histogram
    int toHighSamples = 0;      ///< number of registered samples that were to high for the histogram
    std::vector<int> bins;      ///< bins (vector size: nBins)

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