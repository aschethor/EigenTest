#include "test.cpp"
#include <conio.h>

std::string to_string(double a){
    std::ostringstream sstream;
    sstream << a;
    return sstream.str();
}

void original(){
    printf("=========== original ===========\n");
    Argon argon(6,6,6,1.7048);
    VelocityVerletIntegrator VVI(&argon,0.003);
    VVI.equilibrate(1000,0.7867);
    MDtracker mDtracker("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\output\\",
                        &VVI,1000,10,"N6x6x6_T0,7867_A1,7048");//Number Temperature BoxSize
    mDtracker.track();
}

void diffTemp(){
    for(int i=0;i<10;i++){
        printf("=========== different Temperatures: %d of 10 ===========\n",i+1);
        Argon argon(6,6,6,1.7048);
        VelocityVerletIntegrator VVI(&argon,0.003);
        VVI.equilibrate(1000,0.5+0.05*i);
        MDtracker mDtracker("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\output\\",
                            &VVI,1000,10,"N6x6x6_T"+to_string(0.5+i*0.05)+"_A1,7048");//Number Temperature BoxSize
        mDtracker.track();
    }
}

void diffDens(){
    for(int i=0;i<10;i++){
        printf("=========== different Densities: %d of 10 ===========\n",i+1);
        Argon argon(6,6,6,1.2+0.1*i);
        VelocityVerletIntegrator VVI(&argon,0.003);
        VVI.equilibrate(1000,0.7867);
        MDtracker mDtracker("C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\output\\",
                            &VVI,1000,10,"N6x6x6_T0,7867_A"+to_string(1.2+0.1*i));//Number Temperature BoxSize
        mDtracker.track();
    }
}

void diffBoxSize(){
    for(int i=0;i<5;i++) {
        printf("=========== different Boxsizes: %d of 5 ===========\n",i+1);
        Argon argon(4+i, 4+i, 4+i, 1.7048);
        VelocityVerletIntegrator VVI(&argon, 0.003);
        VVI.equilibrate(1000, 0.7867);
        MDtracker mDtracker(
                "C:\\Users\\Ni\\Desktop\\Uni\\9.Semester\\Computer Simulation of Physical Systems I\\output\\",
                &VVI, 1000, 10, "N"+to_string(4+i)+"x"+to_string(4+i)+"x"+to_string(4+i)+"_T0,7867_A1,7048");//Number Temperature BoxSize
        mDtracker.track();
    }
}

int main() {
    //original();
    //diffTemp();
    diffDens();
    diffBoxSize();
    getch();
    return 0;
}