/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "useful.h"
#include <iostream>

using namespace std;

DFT::DFT() {

}
 

void DFT::setup(){

    Gv1.resize(Nz);
    Gv2.resize(Nz);
    field.resize(Nz);


    n0.resize(Nz);
    n1.resize(Nz);
    n2.resize(Nz);
    n3.resize(Nz);
    nv1.resize(Nz);
    nv2.resize(Nz);


    cn0.resize(Nz);
    cn1.resize(Nz);
    cn2.resize(Nz);
    cn3.resize(Nz);
    cnv1.resize(Nz);
    cnv2.resize(Nz);

    dphi0.resize(Nz);
    dphi1.resize(Nz);
    dphi2.resize(Nz);
    dphi3.resize(Nz);
    dphiv1.resize(Nz);
    dphiv2.resize(Nz);

    cdphi0.resize(Nz);
    cdphi1.resize(Nz);
    cdphi2.resize(Nz);
    cdphi3.resize(Nz);
    cdphiv1.resize(Nz);
    cdphiv2.resize(Nz);

    c.resize(Nz);
    cc.resize(Nz);
    V.resize(Nz);


    density.resize(Nz);
    coldensity1.resize(Nz);
    G1.resize(Nz, Ns);
    G2.resize(Nz, Ns);

    Zero_vec(Gv1, Nz);
    Zero_vec(Gv2, Nz);

    Zero_vec(field, Nz);
    Zero_vec(density, Nz);
    Zero_vec(coldensity1, Nz);

    Zero_vec(n0, Nz);
    Zero_vec(n1, Nz);
    Zero_vec(n2, Nz);
    Zero_vec(n3, Nz);
    Zero_vec(nv1, Nz);
    Zero_vec(nv2, Nz);


    Zero_vec(cn0, Nz);
    Zero_vec(cn1, Nz);
    Zero_vec(cn2, Nz);
    Zero_vec(cn3, Nz);
    Zero_vec(cnv1, Nz);
    Zero_vec(cnv2, Nz);

    Zero_vec(c, Nz);
    Zero_vec(cc, Nz);

    Zero_mat(G1, Nz, Ns);
    Zero_mat(G2, Nz, Ns);

}

DFT::~DFT() {
    DFT::system_out_file.close();
}

void DFT::evolve() {

    iter = 0;
    conver = 1;

    init_field(0.80); // Initialize mean field
    init_coldensity1(20); // Initialize colloid density (bulk?)
    comp_POT(); // Compute external potential. 

    export_data();
    
    while (conver > gamma) // Do while un converged. Well durr. 
    {
        solveGs(); // Solve for the propagators for polymer density 
        comp_dens(); // Compute density(Z)
        if (conver == 1) {
            DFT::comp_n_pol();
        }
        
        if(conver ==1){
            DFT::comp_n_col1();
        }
        //  norm();
        update_mf(); //Update mean field now. Argument is 1 for FMT hard sphere stuff. DO NOT set to 0. Plez.
        update_col1(); // Update the colloid density

        //export_data();

        DFT::system_out_file << "MF convergence: " << conver << endl;
        DFT::system_out_file << "Col density 1 convergence: " << conver_col1 << endl;
        //        cout << "Col density 2 convergence: " << conver_col2 << endl;

        iter++;
    }

   // export_data();
    
}

/*
 Export essential output data from the simulation
 */
void DFT::export_data() {

    DFT::poly_dens_file.open(DFT::poly_dens_filename);
    DFT::col1_dens_file.open(DFT::col1_dens_filename);
    DFT::meanfield_file.open(DFT::meanfield_filename);
    DFT::external_pot_file.open(DFT::external_pot_filename);

    for (int i = 0; i < Nz; i++) {
        DFT::poly_dens_file << (db) i * dz << "\t" << DFT::density(i) << endl;
        DFT::col1_dens_file << (db) i * dz << "\t" << DFT::coldensity1(i) << endl;
        DFT::meanfield_file << (db) i * dz << "\t" << DFT::field(i) << endl;
        DFT::external_pot_file << (db) i * dz << "\t" << DFT::V(i) << endl;
    }

    DFT::poly_dens_file.close();
    DFT::col1_dens_file.close();
    DFT::meanfield_file.close();
    DFT::external_pot_file.close();

}

void DFT::update_col1() {

    db ARG = 0, max = 0, old_d = 0, diff = 0;

    comp_FMT_col1();

    for (int i = 0; i < Nz; i++) {
        old_d = coldensity1(i);
        ARG = chem1 + cc(i) - V(i) - DFT::comp_att_term(i, density, epc1, r, rc1, lambdapc1); // double check if cc(i) is negative
        coldensity1(i) = (1.0 - dt) * coldensity1(i) + dt * colbulk1 * exp(ARG);
        diff = fabs(old_d - coldensity1(i));
        if (diff > max)
            max = diff;

    }

    conver_col1 = max;

}

void DFT::update_mf() {
    hs = 0;

    att = 0;
    db max = 0;
    db old_mf = 0, diff = 0;

    comp_FMT_pol();


    for (int i = 0; i < Nz; i++) {
        old_mf = field(i);

        field(i) = old_mf + dt * (-old_mf + c(i) + V(i) + DFT::comp_att_term(i, density, epp, r, r, lambdapp));

        diff = fabs(old_mf - field(i));
        if (diff > max)
            max = diff;

    }

    conver = max;
}

void DFT::norm() {
    
    vec d1;
    d1.resize(Nz);
    db unnorm = 0, norm = 0;
    for (int i = 0; i < Nz; i++) {
        d1(i) = coldensity1(i);
        unnorm += d1(i) * dz*A;
    }
    for (int i = 0; i < Nz; i++) {

        d1(i) = Nc1 * d1(i) / unnorm;
        coldensity1(i) = d1(i);
        // d2 << (db)i*dz << "\t" << d1(i)<<endl;
        norm += d1(i) * dz*A;
    }

    //   cout << "TOTAL N: " << norm << endl;
    //d2.close();
}

void DFT::init_coldensity1(unsigned int cut) // Initialize the colloid density
{
    for (int i = 0; i < Nz; i++) {

        if (i <= cut || i == Nz - 1) {
            coldensity1(i) = 0.0;
        } else if (i > cut)
            coldensity1(i) = colbulk1 * exp(V(i));
    }
}

void DFT::init_field(db a) {

    for (int i = 0; i < Nz; i++) {
        field(i) = a;
    }
}

db DFT::correct(db R, db x) {
    if (eq(x, R, 0.001)) {
        return frac(3.0, 8.0);
    } else if (eq(x, R - dz, 0.001)) {
        return frac(7.0, 6.0);
    } else if (eq(x, R - 2.0 * dz, 0.001)) {
        return frac(23.0, 24.0);
    } else {
        return 1.0;
    }
    return 1.0;
}

