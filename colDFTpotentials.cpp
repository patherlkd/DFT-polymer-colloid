/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



#include "colDFT.h"

void DFT::comp_POT() {

    db z;

    for (int i = 0; i < Nz; i++) {
        z = (db) i*dz;

        if (z > dia && z < (h - dia))
            V(i) = 0.0;
        else
            V(i) = wall_strength*(exp((-z) / (2.0 * dia)) + exp((z - h) / (2.0 * dia)) - exp(-dia / dia));

        gnuV.send2file(z, V(i));

    }
    // gnuV.plot('l');

}


db DFT::attractive(int place) {
    db pot = 0, delz = 0, sum = 0, term = 0;
    int n = Nz;
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
