/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"


void DFT::comp_FMT_pol() {



    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_pol(); // Compute the weighted densities
    // comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 
    comp_dphi_pol(); // White bear
    //Zero_vec(c,Nz);


    for (int i = 0; i < Nz; i++) // Convolution in real space because lol. 
    {
        W0(i) = 0;
        W1(i) = 0;
        W2(i) = 0;
        W3(i) = 0;
        Wv1(i) = 0;
        Wv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            W0(i) += simp(j)*(dphi0(j)) * w0(r, (db) i * dz - (db) j * dz) * dz;
            W1(i) += simp(j)*(dphi1(j)) * w1(r, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(dphi2(j)) * w2(r, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(dphi3(j)) * w3(r, (db) i * dz - (db) j * dz) * dz;
            Wv1(i) += simp(j)*(dphiv1(j))*-wv1(r, ((db) i * dz - (db) j * dz)) * dz; // N.B negative sign since here z - z' -> z' -z
            Wv2(i) += simp(j)*(dphiv2(j))*-wv2(r, ((db) i * dz - (db) j * dz)) * dz; // also here
        }

    }

    c = (W0 + W1 + W2 + W3 + Wv1 + Wv2); // this is the thing to put in the update of mean field

    //   test << W2 <<endl;


}

void DFT::comp_FMT_col() {

    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_col(); // Compute the weighted densities
    //comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 
    comp_dphi_col(); // White bear
    //Zero_vec(c,Nz);


    for (int i = 0; i < Nz; i++) // Convolution in real space because lol. 
    {
        W0(i) = 0;
        W1(i) = 0;
        W2(i) = 0;
        W3(i) = 0;
        Wv1(i) = 0;
        Wv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            W0(i) += simp(j)*(cdphi0(j)) * w0(rc, (db) i * dz - (db) j * dz) * dz;
            W1(i) += simp(j)*(cdphi1(j)) * w1(rc, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(cdphi2(j)) * w2(rc, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(cdphi3(j)) * w3(rc, (db) i * dz - (db) j * dz) * dz;
            Wv1(i) += simp(j)*(cdphiv1(j))*-wv1(rc, ((db) i * dz - (db) j * dz)) * dz; // N.B negative sign since here z - z' -> z' -z
            Wv2(i) += simp(j)*(cdphiv2(j))*-wv2(rc, ((db) i * dz - (db) j * dz)) * dz; // also here
        }

    }

    cc = (W0 + W1 + W2 + W3 + Wv1 + Wv2); // this is the thing to put in the update of mean field

    //   test << W2 <<endl;

}



