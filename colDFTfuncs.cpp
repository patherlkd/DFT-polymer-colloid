/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "useful.h"
#include <iostream>

void DFT::comp_n_pol() {


    for (int i = 0; i < Nz; i++) {
        static db den = 0;
        n0(i) = 0;
        n1(i) = 0;
        n2(i) = 0;
        n3(i) = 0;
        nv1(i) = 0;
        nv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            den = density(j);

            n0(i) += simp(j) * den * w0(r, (db) i * dz - (db) j * dz) * dz;
            n1(i) += simp(j) * den * w1(r, (db) i * dz - (db) j * dz) * dz;
            n2(i) += simp(j) * den * w2(r, (db) i * dz - (db) j * dz) * dz;
            n3(i) += simp(j) * den * w3(r, (db) i * dz - (db) j * dz) * dz;
            nv1(i) += simp(j) * den * wv1(r, (db) i * dz - (db) j * dz) * dz;
            nv2(i) += simp(j) * den * wv2(r, (db) i * dz - (db) j * dz) * dz;

        }

    }

}

void DFT::comp_dphi_pol_v2() {

    db N = (db) Nm;


    db nf0, nf1, nf2, nf3, nfv1, nfv2;

    db np0, np1, np2, np3, npv1, npv2;

    db Rf = rc1, Rp = r;

    for (int i = 0; i < Nz; i++) {

        nf0 = cn0(i);
        nf1 = cn1(i);
        nf2 = cn2(i);
        nf3 = cn3(i);
        nfv1 = cnv1(i);
        nfv2 = cnv2(i);

        np0 = n0(i);
        np1 = n1(i);
        np2 = n2(i);
        np3 = n3(i);
        npv1 = nv1(i);
        npv2 = nv2(i);

        

    }

}

void DFT::comp_FMT_pol() {



    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_pol(); // Compute the weighted densities
    // comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 

    // comp_dphi_pol(); // White bear old version
    comp_dphi_pol_v2();



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

void DFT::comp_n_col1() {

    for (int i = 0; i < Nz; i++) {
        static db den = 0;
        cn0(i) = 0;
        cn1(i) = 0;
        cn2(i) = 0;
        cn3(i) = 0;
        cnv1(i) = 0;
        cnv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            den = coldensity1(j);

            cn0(i) += simp(j) * den * w0(rc1, (db) i * dz - (db) j * dz) * dz;
            cn1(i) += simp(j) * den * w1(rc1, (db) i * dz - (db) j * dz) * dz;
            cn2(i) += simp(j) * den * w2(rc1, (db) i * dz - (db) j * dz) * dz;
            cn3(i) += simp(j) * den * w3(rc1, (db) i * dz - (db) j * dz) * dz;
            cnv1(i) += simp(j) * den * wv1(rc1, (db) i * dz - (db) j * dz) * dz;
            cnv2(i) += simp(j) * den * wv2(rc1, (db) i * dz - (db) j * dz) * dz;

        }

    }

}

void DFT::comp_n_col2() {

    for (int i = 0; i < Nz; i++) {
        static db den = 0;
        ccn0(i) = 0;
        ccn1(i) = 0;
        ccn2(i) = 0;
        ccn3(i) = 0;
        ccnv1(i) = 0;
        ccnv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            den = coldensity2(j);

            ccn0(i) += simp(j) * den * w0(rc2, (db) i * dz - (db) j * dz) * dz;
            ccn1(i) += simp(j) * den * w1(rc2, (db) i * dz - (db) j * dz) * dz;
            ccn2(i) += simp(j) * den * w2(rc2, (db) i * dz - (db) j * dz) * dz;
            ccn3(i) += simp(j) * den * w3(rc2, (db) i * dz - (db) j * dz) * dz;
            ccnv1(i) += simp(j) * den * wv1(rc2, (db) i * dz - (db) j * dz) * dz;
            ccnv2(i) += simp(j) * den * wv2(rc2, (db) i * dz - (db) j * dz) * dz;

        }

    }

}



void DFT::comp_dphi_col1_v2() {

    db N = (db) Nm;


    db nf0, nf1, nf2, nf3, nfv1, nfv2;

    db np0, np1, np2, np3, npv1, npv2;

    db Rf = rc1, Rp = r;



    for (int i = 0; i < Nz; i++) {

        nf0 = cn0(i);
        nf1 = cn1(i);
        nf2 = cn2(i);
        nf3 = cn3(i);
        nfv1 = cnv1(i);
        nfv2 = cnv2(i);

        np0 = n0(i);
        np1 = n1(i);
        np2 = n2(i);
        np3 = n3(i);
        npv1 = nv1(i);
        npv2 = nv2(i);



        cdphi0(i) = -log(1.0 - nf3 - np3);



        cdphi1(i) = (nf2 + np2) / (1.0 - nf3 - np3);


        cdphi2(i) = (nf1 + np1) / (1 - nf3 - np3) +
                ((2 * (nf2 + np2) - 3 * Power(nfv2 + npv2, 2))*
                (nf3 + np3 + Power(1 - nf3 - np3, 2) *
                log(1 - nf3 - np3))) /
                (36. * Power(1 - nf3 - np3, 2)*(nf3 + np3) * pi); // HS term

        cdphi3(i) = -((-nf0 - np0) / (1 - nf3 - np3)) +
                ((nf1 + np1)*(nf2 + np2) -
                (nfv1 + npv1)*(nfv2 + npv2)) / Power(1 - nf3 - np3, 2)
                + ((Power(nf2 + np2, 2) -
                3 * (nf2 + np2) * Power(nfv2 + npv2, 2))*
                (nf3 + np3 - 2 * (1 - nf3 - np3) * log(1 - nf3 - np3))) /
                (36. * Power(1 - nf3 - np3, 2)*(nf3 + np3) * pi) -
                ((Power(nf2 + np2, 2) -
                3 * (nf2 + np2) * Power(nfv2 + npv2, 2))*
                (nf3 + np3 + Power(1 - nf3 - np3, 2) *
                log(1 - nf3 - np3))) /
                (36. * Power(1 - nf3 - np3, 2) * Power(nf3 + np3, 2) * pi) +
                ((Power(nf2 + np2, 2) -
                3 * (nf2 + np2) * Power(nfv2 + npv2, 2))*
                (nf3 + np3 + Power(1 - nf3 - np3, 2) *
                log(1 - nf3 - np3))) /
                (18. * Power(1 - nf3 - np3, 3)*(nf3 + np3) * pi); // HS



        cdphiv1(i) = (-nfv2 - npv2) / (1. - nf3 - np3);

        cdphiv2(i) = (-nfv1 - npv1) / (1 - nf3 - np3) -
                ((nf2 + np2)*(nfv2 + npv2)*
                (nf3 + np3 + Power(1 - nf3 - np3, 2) *
                log(1 - nf3 - np3))) /
                (6. * Power(1 - nf3 - np3, 2)*(nf3 + np3) * pi);

        if (!DFT::polymers_off) {

            cdphi2(i) += ((1 - N) * np0 * (1 - Power(npv2, 2) / Power(np2, 2))*
                    ((2 * (1 - Power(nfv2 + npv2, 2)) * Power(Rf, 2) *
                    Power(Rp, 2)) /
                    (9. * Power(1 - nf3 - np3, 3) * Power(Rf + Rp, 2)) +
                    (2 * Power(nfv2 + npv2, 2) * Rf * Rp) /
                    (Power(nf2 + np2, 2) * Power(1 - nf3 - np3, 2)*
                    (Rf + Rp)) + ((1 -
                    Power(nfv2 + npv2, 2) / Power(nf2 + np2, 2)) * Rf * Rp
                    ) / (Power(1 - nf3 - np3, 2)*(Rf + Rp)))) /
                    (N * (1 / (1 - nf3 - np3) +
                    (2 * (nf2 + np2)*(1 - Power(nfv2 + npv2, 2)) *
                    Power(Rf, 2) * Power(Rp, 2)) /
                    (9. * Power(1 - nf3 - np3, 3) * Power(Rf + Rp, 2)) +
                    ((nf2 + np2)*(1 -
                    Power(nfv2 + npv2, 2) / Power(nf2 + np2, 2)) * Rf * Rp
                    ) / (Power(1 - nf3 - np3, 2)*(Rf + Rp)))); // CH term


            cdphi3(i) += ((1 - N) * np0 * (1 - Power(npv2, 2) / Power(np2, 2))*
                    (Power(1 - nf3 - np3, -2) +
                    (2 * (nf2 + np2)*(1 - Power(nfv2 + npv2, 2)) *
                    Power(Rf, 2) * Power(Rp, 2)) /
                    (3. * Power(1 - nf3 - np3, 4) * Power(Rf + Rp, 2)) +
                    (2 * (nf2 + np2)*(1 -
                    Power(nfv2 + npv2, 2) / Power(nf2 + np2, 2)) * Rf * Rp
                    ) / (Power(1 - nf3 - np3, 3)*(Rf + Rp)))) /
                    (N * (1 / (1 - nf3 - np3) +
                    (2 * (nf2 + np2)*(1 - Power(nfv2 + npv2, 2)) *
                    Power(Rf, 2) * Power(Rp, 2)) /
                    (9. * Power(1 - nf3 - np3, 3) * Power(Rf + Rp, 2)) +
                    ((nf2 + np2)*(1 -
                    Power(nfv2 + npv2, 2) / Power(nf2 + np2, 2)) * Rf * Rp
                    ) / (Power(1 - nf3 - np3, 2)*(Rf + Rp)))); // CH term


            cdphiv2(i) += ((1 - N) * np0 * (1 - Power(npv2, 2) / Power(np2, 2))*
                    ((-4 * (nf2 + np2)*(nfv2 + npv2) * Power(Rf, 2) *
                    Power(Rp, 2)) /
                    (9. * Power(1 - nf3 - np3, 3) * Power(Rf + Rp, 2)) -
                    (2 * (nfv2 + npv2) * Rf * Rp) /
                    ((nf2 + np2) * Power(1 - nf3 - np3, 2)*(Rf + Rp)))) /
                    (N * (1 / (1 - nf3 - np3) +
                    (2 * (nf2 + np2)*(1 - Power(nfv2 + npv2, 2)) *
                    Power(Rf, 2) * Power(Rp, 2)) /
                    (9. * Power(1 - nf3 - np3, 3) * Power(Rf + Rp, 2)) +
                    ((nf2 + np2)*(1 -
                    Power(nfv2 + npv2, 2) / Power(nf2 + np2, 2)) * Rf * Rp
                    ) / (Power(1 - nf3 - np3, 2)*(Rf + Rp)))); // CH term
        }

    }


}

void DFT::comp_FMT_col1() {

    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_col1(); // Compute the weighted densities
    //comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 
    // comp_dphi_col1(); // White bear
    //Zero_vec(c,Nz);

    comp_dphi_col1_v2();

    for (int i = 0; i < Nz; i++) // Convolution in real space because lol. 
    {
        W0(i) = 0;
        W1(i) = 0;
        W2(i) = 0;
        W3(i) = 0;
        Wv1(i) = 0;
        Wv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            W0(i) += simp(j)*(cdphi0(j)) * w0(rc1, (db) i * dz - (db) j * dz) * dz;
            W1(i) += simp(j)*(cdphi1(j)) * w1(rc1, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(cdphi2(j)) * w2(rc1, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(cdphi3(j)) * w3(rc1, (db) i * dz - (db) j * dz) * dz;
            Wv1(i) += simp(j)*(cdphiv1(j))*-wv1(rc1, ((db) i * dz - (db) j * dz)) * dz; // N.B negative sign since here z - z' -> z' -z
            Wv2(i) += simp(j)*(cdphiv2(j))*-wv2(rc1, ((db) i * dz - (db) j * dz)) * dz; // also here
        }

    }

    cc = (W0 + W1 + W2 + W3 + Wv1 + Wv2); // this is the thing to put in the update of mean field

    //   test << W2 <<endl;

}


void DFT::comp_FMT_col2() {

    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_col2(); // Compute the weighted densities
    //comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 
    // comp_dphi_col1(); // White bear
    //Zero_vec(c,Nz);

    comp_dphi_col1_v2();

    for (int i = 0; i < Nz; i++) // Convolution in real space because lol. 
    {
        W0(i) = 0;
        W1(i) = 0;
        W2(i) = 0;
        W3(i) = 0;
        Wv1(i) = 0;
        Wv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            W0(i) += simp(j)*(ccdphi0(j)) * w0(rc2, (db) i * dz - (db) j * dz) * dz;
            W1(i) += simp(j)*(ccdphi1(j)) * w1(rc2, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(ccdphi2(j)) * w2(rc2, (db) i * dz - (db) j * dz) * dz;
            W3(i) += simp(j)*(ccdphi3(j)) * w3(rc2, (db) i * dz - (db) j * dz) * dz;
            Wv1(i) += simp(j)*(ccdphiv1(j))*-wv1(rc2, ((db) i * dz - (db) j * dz)) * dz; // N.B negative sign since here z - z' -> z' -z
            Wv2(i) += simp(j)*(ccdphiv2(j))*-wv2(rc2, ((db) i * dz - (db) j * dz)) * dz; // also here
        }

    }

    ccc = (W0 + W1 + W2 + W3 + Wv1 + Wv2); // this is the thing to put in the update of mean field

    //   test << W2 <<endl;

}



