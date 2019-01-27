/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



#include "colDFT.h"
#include "math.h"
#include "useful.h"

void DFT::comp_POT() {

    db z;

    for (int i = 0; i < Nz; i++) {
        z = (db) i*dz;

        if (z > dia && z < (h - dia))
            V(i) = 0.0;
        else
            V(i) = wall_strength * (exp((-z) / (2.0 * dia)) + exp((z - h) / (2.0 * dia)) - exp(-dia / dia));
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

        pot_below_d = lukes_potential(0.0, eppij, dij, lambdaij); // compute this once

        for (int iz = 0; iz < Nz; iz++) {

            delz = fabs((db) (iz - zplace) * dz);

            if (delz >= dij) {
                sum += simp(iz) * densityj(iz) * lukes_potential(delz, eppij, dij, lambdaij) * dz;
            } else {
                sum += simp(iz) * densityj(iz) * pot_below_d*dz;
            }

        }

    }

    return sum;

}

db DFT::dinos_potential(db Z, db eppij, db dij, db lambdaij) {

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

    return (2 * eppij * pi * (4 * dij3 / 3. - +2 * dij2 * lambdaij + lambdaij3 - (lambdaij * Z2) / 2.0 + Z3 / 3.0 -
            exp((2 * dij - Z) / lambdaij) * lambdaij2 * (lambdaij + Z) + dij * (2 * lambdaij2 - Z2))) /
            (dij + lambdaij - exp(dij / lambdaij) * lambdaij);
}

/*
db DFT::attractive(int place) {
    db pot = 0, delz = 0, sum = 0, term = 0;

    for (int l = 0; l < Nz; l++) {
        delz = fabs((db) l * dz - (db) place * dz);
        if (delz >= sigma && delz <= cut)// cutoff disabled for Dinos potential  ATM! !!!!!!
        {

            //pot=-eps*2.0*pi*lambda*(lambda+delz)*exp(-(delz - sigma)/lambda); // Dino potential
            pot = (eps / conv_fact) * LEA(delz);

        } else if (delz < sigma) {
            pot = (eps / conv_fact) * LEA(sigma);

        } else;
        sum += simp(l) * pot * density(l) * dz;
    }

    return sum;

}

db DFT::LEA(db Z) {
    db Z2 = Z*Z, ans;

    ans = 8.45561 - 3.46718 * Z2 + 1.0138 * Z2 * Z + exp(-1.31579 * Z)*(-9.8651 - 12.9804 * Z);

    return ans;

}


 
 * From DFT.h lemons might be more recent
 * 
 db DFT::attractive(int place) {
    if (eps == 0.0)
        return 0.0;

    db pot = 0.0, delz = 0, sum = 0, term = 0;
    int n = Nz;
    for (int l = 0; l < Nz; l++) {
        delz = fabs((db) l * dz - (db) place * dz);

        if (delz <= cut) {
            pot = LEA((-eps / conv_fact), lambda, delz);
        } else;

        sum += simp(l) * pot * density(l) * dz;
    }

    return sum;

}

db DFT::LEA(db epp, db lambda, db Z)// Potential integrated over the plane                                                      
{
    db Z2 = Z*Z, ans;
    db e1 = exp(Z / lambda);
    db e2 = exp(-3.0 - Z / lambda);
    db e3 = exp(1 + Z / lambda);
    db lambda2 = lambda*lambda;
    db lambda3 = lambda2*lambda;

    //ans = 8.45561 - 3.46718*Z2 +  1.0138*Z2*Z+exp(-1.31579*Z)*(-9.8651-12.9804*Z);                                            
    ans = 0;
    ans = (1.0 / (3.0 * lambda)) * e2 * epp * pi * (-1.0 + Htheta(Z - 2.0 * lambda));
    ans *= ((31.0 - 12 * euler) * e1 * lambda3 + (12.0 * e3 * lambda3 - 6.0 *euler * euler * lambda2 * (Z + lambda) + e1 * ((2.0 * Z2 * Z) - 9.0 * Z2 * lambda + 7.0 * lambda3)) * Htheta(Z - lambda));
    return ans;
}
 */

