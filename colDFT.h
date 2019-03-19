/*

POLYMER COLLOID DFT

Luke Kristopher Davis:  ucapkda@ucl.ac.uk

 */

#ifndef DFT_H
#define DFT_H

#include <fstream>
#include <string>
#include "useful.h"

class DFT {
public:

    DFT();
    ~DFT();
    
    void set_poly_dens_filename(std::string s);
    void set_col1_dens_filename(std::string s);
    void set_meanfield_filename(std::string s);
    void set_external_pot_filename(std::string s);
    void set_system_out_filename(std::string s);

    
    void set_potential_mode(unsigned int pm);
    void set_polymers_off();
    void set_wall_strength(db ws);
    void set_epp(db epp);
    void set_lambdapp(db lpp);
    void set_epc1(db epc1);
    void set_lambdapc1(db lpc);
    void set_ec1c1(db);
    void set_lambdac1c1(db);
    
    void set_dia(db dia);
    void set_Np(unsigned int N);
    void set_Nm(unsigned int N);
    void set_D(db coef);
    void set_b(db B);
    void set_ds(unsigned int N);
    void set_dz(unsigned int N, db Z);
    void set_A(db a);
    void set_gamma(db g);
    void set_dt(db t);
    void set_DT(db DT);
    void set_tether(unsigned int t);
    void set_H_solver();
    void set_chem1(db ch);
    void set_colbulk1();
    void set_ncolloids1(unsigned int nc);
    void set_rc1(db rc);
    void set_col1_init_cut(db);
    
    void test_dinos_potential(unsigned int Nz,db dz, db eppij, db dij, db lambdaij);
    void test_lukes_potential(unsigned int Nz,db dz, db eppij, db dij, db lambdaij);
    
    
    db get_dz();

    
    virtual void setup();
    // Churns whole scheme
    virtual void evolve();

private:
    unsigned int Np; // number of polymers
    unsigned int Nm; // number of monomers per polymer
    unsigned int Ns = 0; // number of s points
    unsigned int Nz = 0; // number of spatial points
    unsigned int iter; // iteration tracker
    unsigned int potential_mode; // to use Dinos or my shorter ranged potential
    bool polymers_off = false;
    
         // pair potentials
    db dinos_potential(db Z, db eppij, db dij, db lambdaij);
    db lukes_potential(db Z, db eppij, db dij, db lambdaij);
    db LEA(db);
    
    void init_field(db);
    void init_coldensity1(unsigned int);

    virtual void solveGs();
    virtual void comp_dens();
    void comp_POT();
    void update_mf();
    void update_col1();
    void CS(int);
    void norm();
    db comp_free_energy();
    db simp(int);
    db attractive(int);
    db comp_att_term(int zplace, vec& densityj, db eppij, db ri, db rj, db lambdaij); // New general attractive term in free energy

    // Fundamental Measure Theory functions

    void comp_FMT_pol();
    void comp_FMT_col1();
    void RF();

    void comp_n_pol();
    void comp_n_col1();
    void comp_dphi_pol();
    void comp_dphi_col1();
    db WBF(unsigned int);
    db CHF(unsigned int);
    db correct(db, db);

    // geometric weights 

    db w0(db, db);
    db w1(db, db);
    db w2(db, db);
    db w3(db, db);
    db wv1(db, db);
    db wv2(db, db);

    // for plotting purposes

    void export_data();

    std::string poly_dens_filename;
    std::string col1_dens_filename;
    std::string meanfield_filename;
    std::string external_pot_filename;
    std::string system_out_filename;

    std::ofstream system_out_file;
    std::ofstream poly_dens_file;
    std::ofstream col1_dens_file;
    std::ofstream meanfield_file;
    std::ofstream external_pot_file;

    // other sim stuff

    db D; //diffusion coefficient
    db H; //used in solving for G's
    db h; //height
    db A; // Area of Lz plane
    db graft; //grafting density
    db ds; //time slice size
    db dz; // space slice size

    db b; // kuhn length
    db dia; // diameter of monomer
    db r; // radius of monomer
    db rc1; // radius of colloid 1


    db gamma; // convergence criterion
    db hs; // one body direct correlation function
    db dt; // timestep for steepest descent (pol)
    db DT; // timestep for colloids
    db wall_strength; // energetic strength of hard wall
    db att; // attactive term contribution

    db epp; // strength for cohesion (polymer-polymer)
    db epc1; // same but for polymer-colloid
    db ec1c1; // colloid 1 colloid 1 interaction
    

    db lambdapp; //range for gaussian attractive (polymer-polymer)
    db lambdapc1; // range for attraction (polymer - colloid 1)
    db lambdac1c1; // range for atraction (col1 - col1)
   
    db unnorm;
    db Z; //" Partition function" 
    db F_ent; // "entropic term Free energy"

    db colbulk1; // bulk fluid density colloids
    db chem1; // excess chem potential
    db col1_init_cut; // for initialising the col1 density 
    int Nc1; // # of colloids

    float conver; // Convergence tracker (mean field) aka polymer
    float conver_col1; // same for colloid 1
    float conver_col2; // same for colloid 2

    int tether; // tethering array value

    vec field, density; //vectors for mean field and the density profile (polymers)
    vec coldensity1; // colloid density
    vec Gv1; // to manipulate greens function at one instance
    vec Gv2; // " "
    vec n0, n1, n2, n3, nv1, nv2; // FMT polymer weighted densities

    vec cn0, cn1, cn2, cn3, cnv1, cnv2; // FMT weighted densities
    vec dphi0, dphi1, dphi2, dphi3, dphiv1, dphiv2; // derivative of free energy density
    vec cdphi0, cdphi1, cdphi2, cdphi3, cdphiv1, cdphiv2; // colloid derivative of free energy density

    vec c; // To hold the one body direct correlation function and the weights for the polymers
    vec cc; // one body... for colloids
    vec V; //external potential
    mat G1; // these will hold all G's for all s and z
    mat G2;

};






/*

 * 
 * Depracated code
 *

//Used for EIGEN LU solver (slower but keep just in case)

for(int i=0;i<Nz;i++)
 {
   for(int j=0;j<Nz;j++)
     {
       if(i==j)
         {
           T(i,j)=1.+(2.*H)+ds*field(i);
           if(i==0 || j==0)
             {
               T(i,j+1)=-H;
             }
           else if(i==Nz-1 || j==Nz-1)
             {
               T(i,j-1)=-H;T(i,j)=1.+(2.*H)+ds*field(i);
             }
           else
             {
               T(i,j-1)=-H;T(i,j+1)=-H;
             }
         }
       else if(fabs(i-j)>1)
         {T(i,j)=0;}
       else;
     }    }


for(int i=0; i<Nz;i++) // Neccessary for the Crank Nicolson scheme
     {
       if(i==0)
         {
           r(i)= (1.0 - 2.0*H - ds*field(i))*Gv1(i) + H*(Gv1(i+1));
         }
       else if(i>0 && i<Nz-1)
         {
           r(i)= (1.0 - 2.0*H - ds*field(i))*Gv1(i) + H*(Gv1(i+1) + Gv1(i-1));
         }
       else
       {
         r(i)= (1.0 - 2.0*H - ds*field(i))*Gv1(i) + H*(Gv1(i-1));					  
       }
     }



 */



#endif
