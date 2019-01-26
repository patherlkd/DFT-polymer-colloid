#include "useful.h"
#include "colDFT.h"
#include "pstream.h"
#include <ctime>
#include "CLI11.hpp"
#include <iostream>

int main(int argc, char *argv[]) {

    CLI::App dftapp{"Density functional theory (DFT) for polymers + colloids in a film"};

    DFT sim();
    
    db epc1=0.0;
    
    dftapp.add_option("--epc1",epc1," Polymer & colloid type 1 cohesion strength ");
    
    CLI11_PARSE(dftapp,argc,argv);
    
  /*  sim.set_eps(EPP);
    sim.set_h(HEIGHT);
    sim.set_b(B);
    sim.set_conv_fact(0.0357611);
    sim.set_dia(DIA);
    sim.set_Np(NP);
    sim.set_Nm(NM);
    //sim.set_D(0.16666*sqr_d(1.0)); 

    cout << "Nm = " << NM << '\n';

    sim.set_ds(12000);
    //  sim.set_dz(NZ,NM*DIA*0.5); 
    sim.set_dz(NZ, HEIGHT);
    sim.set_A(LENGTH * LENGTH);
    sim.set_gamma(GAMMA);
    sim.set_dt(DT);
*/
    time_t START = time(NULL);
  //  sim.evolve();
    time_t END = time(NULL);

    double duration = difftime(END, START);
    std::cout << "FINISHED in " << duration << " seconds " << '\n';

    return 0;
}
