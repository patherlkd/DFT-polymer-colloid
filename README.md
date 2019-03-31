# DFT-polymer-colloid
Density Functional Theory Code for a polymer film + colloids

# dftpolymercolloid2 Instructions:

```
Density functional theory (DFT) for polymers + colloids in a film
Usage: ./dftpolymercolloid2 [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --config TEXT=./config.dft  configuration file for dftpolymercolloid


Simulation parameters:
  -p,--potential_mode INT REQUIRED
                               Attractive potential to use for all particles (0 = long ranged gaussian, 1 = short ranged gaussian) 
  -c,--ctol,--convergence_tolerance FLOAT REQUIRED
                               Convergence tolerance for DFT simulation.
  -t,--dt,--timestep FLOAT REQUIRED
                               timestep for the numerical algorithm to advance system to equilibrium (for polymers)
  --DT,--coltimestep FLOAT REQUIRED
                              timestep for the numerical algorithm but for colloids
  --Nz,--spatial_points INT REQUIRED
                              Number of spatial points
  --polymers_off              Switch polymers off
  --topwall_off_c1            Switch top wall off for colloid 1
  --botwall_off_c1            Switch bottom wall off for colloid 1
  --topwall_off_c2            Switch top wall off for colloid 2
  --botwall_off_c2            Switch bottom wall off for colloid 2


Output files:
  --poly_dens_file TEXT        filename for polymer density
  --col1_dens_file TEXT        filename for colloid density 1
  --col2_dens_file TEXT        filename for colloid density 2
  --meanfield_file TEXT        filename for meanfield
  --external_pot_file TEXT REQUIRED
                               filename for external potential
  --system_out_file TEXT REQUIRED
                               filename for system output


Simulation box parameters:
  -H,--height FLOAT REQUIRED   Height (z) of the simulation box
  -A,--area FLOAT REQUIRED     area (x-y) of the simulation box
  --wall_strength FLOAT REQUIRED
                              Strength of wall repulsive wall interactions


Particle interactions:
  --epp FLOAT Needs: --lambdapp
                               Polymer - polymer cohesion strength 
  --epc1 FLOAT Needs: --lambdapc1
                               Polymer - colloid type 1 cohesion strength 
  --epc2 FLOAT Needs: --lambdapc2
                               Polymer - colloid type 2 cohesion strength 
  --ec1c1 FLOAT Needs: --lambdac1c1
                               colloid 1 - colloid 1 cohesion strength 
  --ec2c2 FLOAT Needs: --lambdac2c2
                               colloid 2 - colloid 2 cohesion strength 
  --ec1c2 FLOAT Needs: --lambdac1c2
                               colloid 1 - colloid 2 cohesion strength 
  --lambdapp FLOAT Needs: --epp
                               Polymer - polymer cohesion range 
  --lambdapc1 FLOAT Needs: --epc1
                               Polymer - colloid type 1 cohesion range
  --lambdapc2 FLOAT Needs: --epc2
                               Polymer - colloid type 2 cohesion range
  --lambdac1c1 FLOAT Needs: --ec1c1
                               colloid 1 - colloid type 1 cohesion range
  --lambdac2c2 FLOAT Needs: --ec2c2
                               colloid 2 - colloid type 2 cohesion range
  --lambdac1c2 FLOAT Needs: --ec1c2
                               colloid 1 - colloid type 2 cohesion range


Polymer parameters:
  --poly_bondlength FLOAT Needs: --npolymers --nbeads
                              Bondlength of beads in a polymer
  --poly_diameter FLOAT Needs: --poly_bondlength --npolymers --nbeads
                              diameter of beads in a polymer
  --npolymers INT             Number of polymers in system
  -N,--nbeads INT Needs: --npolymers
                              Number of beads in a polymer
  --Ns INT Needs: --npolymers --nbeads
                              Number of points for solving polymer diffusion equation
  --tether INT Needs: --npolymers --nbeads
                              z tethering index for polymer grafting


Colloid 1 parameters:
  --chem1 FLOAT               Colloid (1) excess chemical potential
  --Nc1 INT                   Colloid (1) number of colloids
  --rc1 FLOAT                 Colloids (1) radius
  --col1_init_cut FLOAT       Colloids (1) height which below density is zero (to begin with)


Colloid 2 parameters:
  --chem2 FLOAT               Colloid (2) excess chemical potential
  --Nc2 INT                   Colloid (2) number of colloids
  --rc2 FLOAT                 Colloids (2) radius
  --col2_init_cut FLOAT       Colloids (2) height which below density is zero (to begin with)
```
