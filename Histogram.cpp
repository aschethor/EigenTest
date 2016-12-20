//
// Created by Ni on 18.10.2016.
//

#include "Histogram.h"

using namespace std;

Histogram::Histogram(int n,double binSize,double binStart){
    nBins = n;
    this->binSize = binSize;
    this->binStart = binStart;
    bins.resize(n);
}

void Histogram::addSample(int pos){
    bins[pos]++;
}

void Histogram::addSample(double pos){
    if(((int)((pos-binStart)/binSize))<0){
        toLowSamples++;
        return;
    }
    if(((int)((pos-binStart)/binSize))>=nBins){
        toHighSamples++;
        return;
    }
    bins[((int)((pos-binStart)/binSize))]++;
    nSamples++;
}

int Histogram::getBin(int n){
    return bins[n];
}

double Histogram::getNormalizedBin(int n){
    return ((double)bins[n])/nSamples;
}

double Histogram::getPos(int n){
    return n*binSize+binStart;
}

int Histogram::getToHigh(){
    return toHighSamples;
}

int Histogram::getToLow(){
    return toLowSamples;
}

Histogram& Histogram::operator+=(Histogram& h){
    if(nBins!=h.nBins)return *this;
    if(binSize!=h.binSize)return *this;
    if(binStart!=h.binStart)return *this;
    nSamples+=h.nSamples;
    toLowSamples+=h.toLowSamples;
    toHighSamples+=h.toHighSamples;
    for(int i=0;i<nBins;i++){
        bins[i]+=h.bins[i];
    }
    return *this;
}

vector<double> Histogram::getNormalizedBins(){
    vector<double> normalizedBins;
    normalizedBins.resize(nBins);
    for(int i=0;i<nBins;i++)normalizedBins[i]=getNormalizedBin(i);
}