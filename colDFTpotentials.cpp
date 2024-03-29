/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



#include "colDFT.h"
#include "math.h"
#include "useful.h"
#include <string>

void DFT::comp_POT() {

    db z = 0.0;

    for (int i = 0; i < (int) (Nz * 0.5); i++) {
        z = (db) i*dz;

        if (z > dia && z < (h - dia))
            V(i) = 0.0;
        else
            V(i) = wall_strength * (-1.0 * (1.0 / exp(0.5)) - exp((dia - h) / (2.0 * dia)) + exp(-z / (2.0 * dia)) + exp((z - h) / (2.0 * dia)));
        //   V(i) = wall_strength * (exp((-z) / (2.0 * dia)) + exp((z - h) / (2.0 * dia)) - exp(-dia / dia));
    }

    for (int i = (int) (Nz * 0.5); i < Nz; i++) {
        z = (db) i*dz;

        if (z > dia && z < (h - dia))
            V(i) = 0.0;
        else
            V(i) = V(Nz - i - 1);
    }


}

void DFT::comp_POT_c1() {


    db z = 0.0;

    db dc1 = rc1 * 2.0;

    for (int i = 0; i < (int) (Nz * 0.5); i++) {
        z = (db) i*dz;
        if (botwall_off_c1) {
            Vc1(i) = 0.0;
        } else {
             if (z > dc1 && z < (h - dc1))
                Vc1(i) = 0.0;
            else
                Vc1(i) = wall_strength * (-1.0 * (1.0 / exp(0.5)) - exp((dc1 - h) / (2.0 * dc1)) + exp(-z / (2.0 * dc1)) + exp((z - h) / (2.0 * dc1)));
            //   V(i) = wall_strength * (exp((-z) / (2.0 * dia)) + exp((z - h) / (2.0 * dia)) - exp(-dia / dia));
        }
    }

    for (int i = (int) (Nz * 0.5); i < Nz; i++) {
        z = (db) i*dz;
        if (topwall_off_c1) {
            Vc1(i) = 0.0;
        } else {
            if (z > dc1 && z < (h - dc1))
                Vc1(i) = 0.0;
            else
                Vc1(i) = Vc1(Nz - i - 1);
        }
    }


}


void DFT::comp_POT_c2() {


    db z = 0.0;

    db dc2 = rc2 * 2.0;

    for (int i = 0; i < (int) (Nz * 0.5); i++) {
        z = (db) i*dz;
        if (botwall_off_c2) {
            Vc2(i) = 0.0;
        } else {
             if (z > dc2 && z < (h - dc2))
                Vc2(i) = 0.0;
            else
                Vc2(i) = wall_strength * (-1.0 * (1.0 / exp(0.5)) - exp((dc2 - h) / (2.0 * dc2)) + exp(-z / (2.0 * dc2)) + exp((z - h) / (2.0 * dc2)));
            //   V(i) = wall_strength * (exp((-z) / (2.0 * dia)) + exp((z - h) / (2.0 * dia)) - exp(-dia / dia));
        }
    }

    for (int i = (int) (Nz * 0.5); i < Nz; i++) {
        z = (db) i*dz;
        if (topwall_off_c2) {
            Vc2(i) = 0.0;
        } else {
            if (z > dc2 && z < (h - dc2))
                Vc2(i) = 0.0;
            else
                Vc2(i) = Vc2(Nz - i - 1);
        }
    }


}

db DFT::comp_att_term(int zplace, vec &densityj, db eppij, db ri, db rj, db lambdaij) {

    db delz = 0.0, sum = 0.0, dij = ri + rj;
    db pot_below_d = 0.0; // potential energy beyond hard surface contact

    if (potential_mode == 0) { // Dinos potential

        pot_below_d = dinos_potential(0.0, eppij, dij, lambdaij); // compute this once

        for (int iz = 0; iz < Nz; iz++) {

            delz = fabs((db) (iz - zplace) * dz);

            if (delz >= dij) {
                sum += simp(iz) * densityj(iz) * dinos_potential(delz, eppij, dij, lambdaij) * dz;
            } else {
                sum += simp(iz) * densityj(iz) * pot_below_d*dz;
            }

        }


    } else if (potential_mode == 1) { // Lukes shorter ranged potential

        pot_below_d = lukes_potential(dij, eppij, dij, lambdaij); // compute this once

        for (int iz = 0; iz < Nz; iz++) {

            delz = fabs((db) (iz - zplace) * dz);

            if (delz > 2 * dij) {
                continue;
            } else if (delz >= dij) {
                sum += simp(iz) * densityj(iz) * lukes_potential(delz, eppij, dij, lambdaij) * dz;
            } else {
                sum += simp(iz) * densityj(iz) * pot_below_d*dz;
            }

        }

    }

    return sum;

}

db DFT::dinos_potential(db Z, db eppij, db dij, db lambdaij) {

    if (eppij == 0.0) {
        return 0.0;
    }

    return -2.0 * eppij * lambdaij * (dij + lambdaij) * pi + 2.0 * eppij * pi * (lambdaij * (dij + lambdaij) -
            exp((dij - Z) / lambdaij) * lambdaij * (lambdaij + Z)) * Htheta(-dij + Z);

}

db DFT::lukes_potential(db Z, db eppij, db dij, db lambdaij) {

    db dij2 = dij*dij;
    db dij3 = dij2*dij;
    db lambdaij2 = lambdaij*lambdaij;
    db lambdaij3 = lambdaij2*lambdaij;
    db Z2 = Z*Z;
    db Z3 = Z2*Z;

    if (eppij == 0.0) {
        return 0.0;
    }

    return (2 * -eppij * pi * (4 * dij3 / 3. + 2 * dij2 * lambdaij + lambdaij3 - (lambdaij * Z2 / 2.0) + Z3 / 3.0 -
            exp((2 * dij - Z) / lambdaij) * lambdaij2 * (lambdaij + Z) + dij * (2 * lambdaij2 - Z2))) /
            (dij + lambdaij - exp(dij / lambdaij) * lambdaij);
}

void DFT::test_dinos_potential(unsigned int Nz, db dz, db eppij, db dij, db lambdaij, std::string filename) {

    DFT::system_out_file << "Testing dinos potential\n";

    std::ofstream test;

    test.open(filename);

    for (int i = 0; i < Nz; i++) {
        db z = (db) i * dz;
        test << z << "\t" << DFT::dinos_potential(z, eppij, dij, lambdaij) << "\n";
    }

    test.close();
}

void DFT::test_lukes_potential(unsigned int Nz, db dz, db eppij, db dij, db lambdaij, std::string filename) {

    DFT::system_out_file << "Testing lukes potential\n";

    std::ofstream test;

    test.open(filename);

    db z = 0.0;
    db atdia = DFT::lukes_potential(dij, eppij, dij, lambdaij);
    db pot = 0.0;

    for (int i = 0; i < Nz; i++) {
        z = (db) i * dz;


        if (z <= dij) {
            pot = atdia;
        } else if (z <= 2 * dij) {
            pot = DFT::lukes_potential(z, eppij, dij, lambdaij);
        } else {
            pot = 0.0;
        }

        test << z << "\t" << pot << "\n";
    }

    test.close();
}
