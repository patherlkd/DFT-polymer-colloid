#include "useful.cpp"
#include "colDFT.h"
#include "pstream.h"
#include <ctime>

int main(int argc, char *argv[])
{

  /****** GRAPH showing/ debug options
0 = density + meanfield
1 = c(i) and greens functions
2 = weight functions
3 = weighted densities
4 = functional derivative of free energy density 

  *******/
  DFT sim(0,"./siminfo/SIMINFO","./meanfield/MEANFIELD","./profiles/PROFILES");
  //piccard sim(0);
  sim.set_eps(EPP);
  sim.set_h(HEIGHT);
  sim.set_b(B);
  sim.set_conv_fact(0.0357611);
  sim.set_dia(DIA);
  sim.set_mddia(DIA); 
  sim.set_Np(NP); 
  sim.set_Nm(NM); 
  //sim.set_D(0.16666*sqr_d(1.0)); 

  cout << "Nm = "<< NM <<'\n';
  
  sim.set_ds(12000); 
  //  sim.set_dz(NZ,NM*DIA*0.5); 
  sim.set_dz(NZ,HEIGHT); 
  sim.set_A(LENGTH*LENGTH);
  sim.set_gamma(GAMMA);
  sim.set_dt(DT); 
  
  time_t START = time(NULL);   
  sim.evolve();
  time_t END = time(NULL);   
  
  double duration = difftime(END,START);
  cout<< "FINISHED in "<<duration<<" seconds "<<'\n';
  
  return 0;
}
