/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "iostream"
#include "string"
#include "fstream"

db DFT::comp_free_energy() {

    ifstream math_att("attract_energy_term.txt");
    db sum = 0, F, F_field, F_exc, F_att, F_pot;
    string line;
    int z_int;
    db Ival;
    vec att_;
    att_.resize(Nz);
    Zero_vec(att_, Nz);
    // IDEAL TERM
    Z = ((db) Np * (db) Nm) / unnorm;
    F_ent = -(db) Np * log(Z);
    // Field TERM
    for (int i = 0; i < Nz; i++) {
        sum += field(i) * density(i) * dz;
    }
    F_field = -sum;
    sum = 0;
    // White Bear terms
    for (unsigned int i = 0; i < Nz; i++) {
        sum += (WBF(i) + CHF(i)) * dz;
    }
    F_exc = sum;
    sum = 0;
    // External potential
    for (int i = 0; i < Nz; i++) {
        sum += V(i) * density(i) * dz;
    }
    F_pot = sum;
    sum = 0;
    // Attractive terms
    int cnt = 0;
    while (getline(math_att, line)) {

        stringstream ss(line);

        ss >> z_int >> Ival;
        if (cnt < Nz)
            att_(cnt) = Ival;
        cnt++;
    }
    cout << "out of getline" << endl;
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nz; j++) {
            sum += density(i) * density(j) * eps * att_(abs(i - j)) * dz*dz;
        }
    }
    F_att = sum;
    sum = 0;

    F = F_ent + F_pot + F_att + F_field + F_exc;

    cout << "Free energy= " << F << endl;

    math_att.close();
    return F;
}