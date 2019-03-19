#include "useful.h"
#include "colDFT.h"
#include "pstream.h"
#include <ctime>
#include "CLI11.hpp"
#include <iostream>

int main(int argc, char *argv[]) {

    CLI::App dftapp{"Density functional theory (DFT) for polymers + colloids in a film"};

    db epp = 0.0, epc1 = 0.0, ec1c1 = 0.0;
    db lambdapp = 0.0, lambdapc1 = 0.0, lambdac1c1 = 0.0;
    db gamma = 0.0, dt = 0.0, DT = 0.0;
    db height = 0.0, area = 0.0, wall_strength = 0.0;

    int potential_mode = 0; // use Dino potential by default
    int Nz = 0;

    std::string poly_dens_filename = "";
    std::string col1_dens_filename = "";
    std::string meanfield_filename = "";
    std::string external_pot_filename = "";
    std::string system_out_filename = "";

    dftapp.set_config("--config", "./config.dft", "configuration file for dftpolymercolloid", true);

    CLI::Option* opt_potentialmode = dftapp.add_option("-p,--potential_mode", potential_mode, " Attractive potential to use for all particles (0 = long ranged gaussian, 1 = short ranged gaussian) ");
    CLI::Option* opt_convertol = dftapp.add_option("-c,--ctol,--convergence_tolerance", gamma, " Convergence tolerance for DFT simulation.");
    CLI::Option* opt_timestep = dftapp.add_option("-t,--dt,--timestep", dt, " timestep for the numerical algorithm to advance system to equilibrium (for polymers)");
    CLI::Option* opt_timestep_col = dftapp.add_option("--DT,--coltimestep", DT, "timestep for the numerical algorithm but for colloids");
    CLI::Option* opt_spacepoints = dftapp.add_option("--Nz,--spatial_points", Nz, "Number of spatial points");

    opt_potentialmode->required()->group("Simulation parameters");
    opt_convertol->required()->group("Simulation parameters");
    opt_timestep->required()->group("Simulation parameters");
    opt_timestep_col->required()->group("Simulation parameters");
    opt_spacepoints->required()->group("Simulation parameters");

    CLI::Option* opt_polydensfilename = dftapp.add_option("--poly_dens_file", poly_dens_filename, " filename for polymer density");
    CLI::Option* opt_col1densfilename = dftapp.add_option("--col1_dens_file", col1_dens_filename, " filename for colloid density");
    CLI::Option* opt_meanfieldfilename = dftapp.add_option("--meanfield_file", meanfield_filename, " filename for meanfield");
    CLI::Option* opt_externalpotfilename = dftapp.add_option("--external_pot_file", external_pot_filename, " filename for external potential");
    CLI::Option* opt_systemoutfilename = dftapp.add_option("--system_out_file", system_out_filename, " filename for system output");

    opt_polydensfilename->group("Output files");
    opt_col1densfilename->group("Output files");
    opt_meanfieldfilename->group("Output files");
    opt_externalpotfilename->group("Output files")->required();
    opt_systemoutfilename->group("Output files")->required();


    CLI::Option* opt_box_height = dftapp.add_option("-H,--height", height, " Height (z) of the simulation box");
    CLI::Option* opt_box_area = dftapp.add_option("-A,--area", area, " area (x-y) of the simulation box");
    CLI::Option* opt_wall_strength = dftapp.add_option("--wall_strength", wall_strength, "Strength of wall repulsive wall interactions");

    opt_box_height->group("Simulation box parameters")->required();
    opt_box_area->group("Simulation box parameters")->required();
    opt_wall_strength->group("Simulation box parameters")->required();

    CLI::Option* opt_epp = dftapp.add_option("--epp", epp, " Polymer - polymer cohesion strength ");
    CLI::Option* opt_epc1 = dftapp.add_option("--epc1", epc1, " Polymer - colloid type 1 cohesion strength ");
    CLI::Option* opt_ec1c1 = dftapp.add_option("--ec1c1", ec1c1, " colloid 1 - colloid 1 cohesion strength ");
    CLI::Option* opt_lambdapp = dftapp.add_option("--lambdapp", lambdapp, " Polymer - polymer cohesion range ");
    CLI::Option* opt_lambdapc1 = dftapp.add_option("--lambdapc1", lambdapc1, " Polymer - colloid type 1 cohesion range");
    CLI::Option* opt_lambdac1c1 = dftapp.add_option("--lambdac1c1", lambdac1c1, " colloid 1 - colloid type 1 cohesion range");

    opt_epp->needs(opt_lambdapp);
    opt_lambdapp->needs(opt_epp);
    opt_epc1->needs(opt_lambdapc1);
    opt_lambdapc1->needs(opt_epc1);
    opt_ec1c1->needs(opt_lambdac1c1);
    opt_lambdac1c1->needs(opt_ec1c1);

    opt_epp->group("Particle interactions");
    opt_epc1->group("Particle interactions");
    opt_lambdapp->group("Particle interactions");
    opt_lambdapc1->group("Particle interactions");
    opt_lambdac1c1->group("Particle interactions");
    opt_ec1c1->group("Particle interactions");

    db poly_bondlength = 0.0;
    db poly_diameter = 0.0;
    int npoly = 0;
    int nbeads = 0;
    int Ns = 0;
    int tether = 0;

    CLI::Option* opt_poly_bondlength = dftapp.add_option("--poly_bondlength", poly_bondlength, "Bondlength of beads in a polymer");
    CLI::Option* opt_poly_diameter = dftapp.add_option("--poly_diameter", poly_diameter, "diameter of beads in a polymer");
    CLI::Option* opt_poly_npoly = dftapp.add_option("--npolymers", npoly, "Number of polymers in system");
    CLI::Option* opt_poly_nbeads = dftapp.add_option("-N,--nbeads", nbeads, "Number of beads in a polymer");
    CLI::Option* opt_poly_ns = dftapp.add_option("--Ns", Ns, "Number of points for solving polymer diffusion equation");
    CLI::Option* opt_poly_tether = dftapp.add_option("--tether", tether, "z tethering index for polymer grafting");

    opt_poly_npoly->group("Polymer parameters");
    opt_poly_bondlength->group("Polymer parameters")->needs(opt_poly_npoly)->needs(opt_poly_nbeads);
    opt_poly_diameter->group("Polymer parameters")->needs(opt_poly_bondlength)->needs(opt_poly_npoly)->needs(opt_poly_nbeads);
    opt_poly_nbeads->group("Polymer parameters")->needs(opt_poly_npoly);
    opt_poly_ns->group("Polymer parameters")->needs(opt_poly_npoly)->needs(opt_poly_nbeads);
    opt_poly_tether->group("Polymer parameters")->needs(opt_poly_npoly)->needs(opt_poly_nbeads);


    db chem1 = 0.0;
    db col1_rad = 0.0;
    db col1_init_cut = 0.0;
    int ncolloids1 = 0;

    CLI::Option* opt_chem1 = dftapp.add_option("--chem1", chem1, "Colloid (1) excess chemical potential");
    CLI::Option* opt_ncolloids1 = dftapp.add_option("--Nc1", ncolloids1, "Colloid (1) number of colloids");
    CLI::Option* opt_col1_rad = dftapp.add_option("--rc1", col1_rad, "Colloids (1) radius");
    CLI::Option* opt_col1_init_cut = dftapp.add_option("--col1_init_cut", col1_init_cut, "Colloids (1) height which below density is zero (to begin with)");


    opt_chem1->group("Colloid parameters");
    opt_ncolloids1->group("Colloid parameters");
    opt_col1_rad->group("Colloid parameters");
    opt_col1_init_cut->group("Colloid parameters");

    CLI11_PARSE(dftapp, argc, argv);

    std::cout << dftapp.config_to_str(true, false) << "\n";

    DFT sim;

    sim.set_poly_dens_filename(poly_dens_filename);
    sim.set_col1_dens_filename(col1_dens_filename);
    sim.set_meanfield_filename(meanfield_filename);
    sim.set_system_out_filename(system_out_filename);
    sim.set_external_pot_filename(external_pot_filename);

    sim.set_potential_mode(potential_mode);
    sim.set_gamma(gamma);
    sim.set_dt(dt);
    sim.set_DT(DT);
    sim.set_dz(Nz, height);
    sim.set_A(area);
    sim.set_wall_strength(wall_strength);



    sim.set_epp(epp);
    sim.set_epc1(epc1);
    sim.set_ec1c1(ec1c1);

    sim.set_lambdapp(lambdapp);
    sim.set_lambdapc1(lambdapc1);
    sim.set_lambdac1c1(lambdac1c1);

    sim.set_Np(npoly);
    sim.set_Nm(nbeads);
    sim.set_dia(poly_diameter);
    sim.set_ds(Ns);
    sim.set_tether(tether);

    sim.set_chem1(chem1);
    sim.set_ncolloids1(ncolloids1);
    sim.set_rc1(col1_rad);
    sim.set_colbulk1();
    sim.set_col1_init_cut(col1_init_cut);

    sim.set_D(0.16666 * sqr_d(poly_diameter));
    sim.set_H_solver(); // VITAL for solving the polymer density

    //sim.test_dinos_potential(Nz, sim.get_dz(), epp, poly_diameter, lambdapp);
    //sim.test_lukes_potential(Nz, sim.get_dz(), epp, poly_diameter, lambdapp);

    time_t START = time(NULL);
    sim.setup();
    sim.evolve();
    time_t END = time(NULL);

    double duration = difftime(END, START);
    std::cout << "FINISHED in " << duration << " seconds " << '\n';

    return 0;
}
