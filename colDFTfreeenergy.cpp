/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "iostream"
#include "string"
#include "fstream"
#include "useful.h"

using namespace std;

db DFT::WBF(unsigned int z) {

    db t1, t2, t3, logn3;
    logn3 = log(1.0 - n3(z));
    t1 = -n0(z) * logn3;
    t2 = frac(n1(z) * n2(z) - nv1(z) * nv2(z), 1.0 - n3(z));
    t3 = (n2(z) * n2(z) * n2(z) - 3.0 * n2(z) * nv2(z) * nv2(z));
    t3 *= frac(n3(z)+(1.0 - n3(z))*(1.0 - n3(z)) * logn3, 36.0 * pi * n3(z) * n3(z)*(1.0 - n3(z)*(1.0 - n3(z))));

    return t1 + t2 + t3;
}

db DFT::CHF(unsigned int z) {

    db f_, logy, eta;
    eta = 1.0 - frac(nv2(z) * nv2(z), n2(z) * n2(z));
    logy = frac(1.0, (1.0 - n3(z)));
    logy += frac(n2(z) * r*eta, 2.0 * (1.0 - n3(z))*(1.0 - n3(z)));
    logy += frac(n2(z) * n2(z) * r * r*eta, 18.0 * (1.0 - n3(z))*(1.0 - n3(z))*(1.0 - n3(z)));
    logy = log(logy);
    f_ = frac(1.0 - (db) Nm, (db) Nm) * n0(z) * eta*logy;
    return f_;
}

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