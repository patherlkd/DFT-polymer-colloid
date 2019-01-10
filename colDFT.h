/*

POLYMER COLLOID DFT

Luke Kristopher Davis:  ucapkda@ucl.ac.uk

 */

#ifndef DFT_H
#define DFT_H


class DFT {
   
public:
    DFT(int);


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

    db get_dz();

    // General DFT functions

    void init_field(db);
    void init_coldensity1();

    virtual void solveGs();
    virtual void comp_dens();
    void comp_POT();
    void update_mf(int);
    void update_col1();
    void CS(int);
    void norm();
    db comp_free_energy();
    db simp(int);
    db attractive(int);
    db comp_att_term(int zplace,vec& densityj,db eppij,db ri,db rj,db lambdaij); // New general attractive term in free energy
    
    // pair potentials
    db dinos_potential(db Z,db eppij,db dij, db lambdaij);
    db LEA(db);

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

    // Churns whole scheme
    virtual void evolve();
    void plotcheck(int); // plots important things to check whats going on (debugging)

private:
    unsigned int Np; // number of polymers
    unsigned int Nm; // number of monomers per polymer
    unsigned int Ns; // number of s points
    unsigned int Nz; // number of spatial points
    unsigned int iter; // iteration tracker
    unsigned int potential_mode; // to use Dinos or my shorter ranged potential

    // for plotting purposes

    gnuplot gnucol, gnudens, gnuw0, gnuw1, gnuw2, gnuw3, gnuwv1, gnuwv2, gnumf;
    gnuplot gnun0, gnun1, gnun2, gnun3, gnunv1, gnunv2;
    gnuplot gnup0, gnup1, gnup2, gnup3, gnupv1, gnupv2;
    gnuplot gnuc, gnug1, gnug2, gnuV, gnuF;
    gnuplot gnuLEA;

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
    db rc1; // radius of colloid
    
    
    db gamma; // convergence criterion
    db hs; // one body direct correlation function
    db dt; // for steepest descent (pol)
    db DT; // for colloids
    db wall_strength; // energetic strength of hard wall
    db att; // attactive term contribution
    db eps; // strength for gaussian attractive
    db lambda; //range for gaussian attractive
    db sigma; // the same as particle diameter
    db cut; // cutoff for attraction
    db conv_fact;
    db unnorm;
    db Z; //" Partition function" 
    db F_ent; // "entropic term Free energy"

    db colbulk; // bulk fluid density colloids
    db chem; // excess chem potential
    int Nc; // # of colloids

    float conver; // Convergence tracker (mean field) aka polymer
    float conver_col1; // same for colloid 1
    float conver_col2; // same for colloid 2

    int tether; // tethering array value
    int deb; // debugging parameter for gnuplot.h
    vec field, density; //vectors for mean field and the density profile (polymers)
    vec coldensity1; // colloid density
    vec Gv1; // to manipulate greens function at one instance
    vec Gv2; // " "
    vec n0, n1, n2, n3, nv1, nv2; // FMT polymer weighted densities

    vec cn0, cn1, cn2, cn3, cnv1, cnv2; // FMT weighted densities
    vec dphi0, dphi1, dphi2, dphi3, dphiv1, dphiv2; // derivative of free energy density
    vec cdphi0, cdphi1, cdphi2, cdphi3, cdphiv1, cdphiv2; // colloid derivative of free energy density

    vec c; // To hold the one body direct correlation function and the weights
    vec cc; // one body... for colloids
    vec V; //external potential
    mat G1; // these will hold all G's for all s and z
    mat G2;


};



/*

 * 
 * Depracated code
 *
db DFT::correct(db R, db x) {
     if(eq(x,R,0.001))
      {return frac(3.0,8.0);}
    else if(eq(x,R - dz,0.001))
      {return frac(7.0,6.0);}
    else if(eq(x,R -2.0*dz,0.001))
      {return frac(23.0,24.0);}
    else 
    {return 1.0;}
    return 1.0;
}

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
