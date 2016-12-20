//
// Created by Ni on 16.10.2016.
//

#include "MDtracker.h"
#include <unsupported/Eigen/FFT>
#include <math.h>

using namespace std;

MDtracker::MDtracker(std::string resultFolder,MDintegrator* MD, int nSamples, int SampleRate,std::string suffix){
    pdbFile = fopen((resultFolder+"track"+suffix+".pdb").c_str(),"w");
    gofrFile = fopen((resultFolder+"gofr"+suffix+".dat").c_str(),"w");
    sofkFile = fopen((resultFolder+"sofk"+suffix+".dat").c_str(),"w");
    msDisplFile = fopen((resultFolder+"msDispl"+suffix+".dat").c_str(),"w");
    vvAutoFile = fopen((resultFolder+"vvAuto"+suffix+".dat").c_str(),"w");
    systemFile = fopen((resultFolder+"system"+suffix+".dat").c_str(),"w");
    this->nSamples = nSamples;
    this->SampleRate = SampleRate;
    this->MD = MD;
    gofrSize = 1000;
    gofr.resize(gofrSize);
    sofkSize = 1000;
    sofk.resize(sofkSize);
    msDispl.resize(nSamples);
    vvAuto.resize(nSamples);
    x_track = new Eigen::MatrixXd[nSamples];
    v_track = new Eigen::MatrixXd[nSamples];
}

MDtracker::~MDtracker() {
    fclose(pdbFile);
    fclose(gofrFile);
    fclose(sofkFile);
    fclose(msDisplFile);
    fclose(vvAutoFile);
    fclose(systemFile);
    delete[] x_track;
    delete[] v_track;
}

void showProgress(int A,int ofB=100){
    if(A==0&&ofB==100)cout<<"   ";
    if(A==100&&ofB==100){
        cout<<"\b\b\bdone."<<endl;
        return;
    }
    int percent = (100*A)/ofB;
    cout<<"\b\b\b";
    if(percent<10)cout<<" ";
    cout<<percent<<"%";
}

void MDtracker::calculateGofr() {
    double drx,dry,drz;
    const double BoxX = MD->MD->Box(0),BoxY = MD->MD->Box(1),BoxZ = MD->MD->Box(2);
    double dr = min(MD->MD->Box(0),min(MD->MD->Box(1),MD->MD->Box(2)))/2/gofrSize;
    Histogram histogram(gofrSize,dr);
    for(int i=0;i<nSamples;i++){
        for(int j=0;j<MD->MD->n;j++){
            for(int k=j+1;k<MD->MD->n;k++){
                //calculate distances
                drx = x_track[i](0,j)-x_track[i](0,k);
                dry = x_track[i](1,j)-x_track[i](1,k);
                drz = x_track[i](2,j)-x_track[i](2,k);
                drx-=BoxX*rint(drx/BoxX);
                dry-=BoxY*rint(dry/BoxY);
                drz-=BoxZ*rint(drz/BoxZ);
                histogram.addSample(sqrt(drx*drx+dry*dry+drz*drz));
            }
        }
    }

    double rho0 = (MD->MD->n-1.0)/(MD->MD->Box(0)*MD->MD->Box(1)*MD->MD->Box(2));
    rho0*=(nSamples*MD->MD->n*4*M_PI*dr*dr*dr/2.0);

    gofr[0]=0;
    for(int i=1;i<gofrSize;i++)gofr[i]=histogram.getBin(i)/(rho0*i*i);

    for(int i=0;i<gofrSize;i++)
        fprintf(gofrFile,"%f %f\n",i*dr,gofr[i]);
}

void MDtracker::calculateSofk(){
    double rho0 = (MD->MD->n-1.0)/(MD->MD->Box(0)*MD->MD->Box(1)*MD->MD->Box(2));
    double dk = maxS/sofkSize;
    double dr = min(MD->MD->Box(0),min(MD->MD->Box(1),MD->MD->Box(2)))/2/gofrSize;
    std::vector<double> fofr;
    fofr.resize(gofrSize);
    for(int i=0;i<gofrSize;i++)fofr[i]=gofr[i]-1;
    rho0*=4*M_PI*dr;
    for(int i=0;i<sofkSize;i++){
        double integral = 0;
        for(int j=0;j<gofrSize;j++){
            if(i==0)integral+=fofr[j]*j*dr*j;
            else integral+=fofr[j]*j*sin(dk*i*dr*j)/(dk*i);
        }
        sofk[i]=1+rho0*integral*dr;
    }

    for(int i=0;i<sofkSize;i++)
        fprintf(sofkFile,"%f %f\n",i*dk,sofk[i]);
}

void MDtracker::calculateMsDispl() {
    double drx,dry,drz;
    const double BoxX = MD->MD->Box(0),BoxY = MD->MD->Box(1),BoxZ = MD->MD->Box(2);
    double msDisplSum;
    for(int i=0;i<nSamples;i++){                //loop over dt
        msDisplSum = 0;
        for(int j=0;j<nSamples-i;j++){          //loop over samples
            for(int k=0;k<MD->MD->n;k++) {      //loop over particles
                //calculate distances
                drx = x_track[j](0, k) - x_track[j+i](0, k);
                dry = x_track[j](1, k) - x_track[j+i](1, k);
                drz = x_track[j](2, k) - x_track[j+i](2, k);
                drx -= BoxX * rint(drx / BoxX);
                dry -= BoxY * rint(dry / BoxY);
                drz -= BoxZ * rint(drz / BoxZ);
                msDisplSum += drx*drx+dry*dry+drz*drz;
            }
        }
        msDispl[i]=msDisplSum/MD->MD->n/(nSamples-i);
    }

    double dt = MD->deltaT*SampleRate;
    for(int i=0;i<nSamples;i++)
        fprintf(msDisplFile,"%f %f\n",i*dt,msDispl[i]);
}

void MDtracker::calculateVvAuto() {
    double vvAutoSum;
    for(int i=0;i<nSamples;i++){
        vvAutoSum = 0;
        for(int j=0;j<nSamples-i;j++){
            for(int k=0;k<MD->MD->n;k++){
                for(int l=0;l<3;l++){
                    vvAutoSum+=v_track[j](l,k)*v_track[j+i](l,k);
                }
            }
        }
        vvAuto[i]=vvAutoSum/MD->MD->n/(nSamples-i);
    }

    double dt = MD->deltaT*SampleRate;
    for(int i=0;i<nSamples;i++)
        fprintf(vvAutoFile,"%f %f\n",i*dt,vvAuto[i]);
}

void MDtracker::track() {
    cout<<"integrating MD-system..."<<endl;

    showProgress(0);
    for(int i=0;i<nSamples;i++){
        showProgress(i,nSamples);
        fprintf(pdbFile,"%s\nENDMDL\n",MD->MD->toPDB().c_str());
        fprintf(systemFile,"%f %f %f %f %f %f\n",i*SampleRate*MD->deltaT,
                MD->MD->Ekin(),MD->MD->Epot(),MD->MD->Etot(),MD->MD->T(),MD->MD->P());
        x_track[i] = MD->MD->x;
        v_track[i] = MD->MD->v;
        for(int j=0;j<SampleRate;j++)
            MD->integrate();
    }
    showProgress(100);

    cout<<"calculating g(r)..."<<endl;
    calculateGofr();
    cout<<"done.\ncalculating s(k)..."<<endl;
    calculateSofk();
    cout<<"done.\ncalculating mean square displacement..."<<endl;
    calculateMsDispl();
    cout<<"done.\ncalculating velocity-velocity auto correlation..."<<endl;
    calculateVvAuto();
    cout<<"done."<<endl;
}