/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include "colDFT.h"
#include "fstream"
#include "useful.cpp"

using namespace std;

void DFT::comp_dens() {
    unnorm = 0;
    bool nan_check(false);
    ofstream dens("NSP1_data/easter_data/packing256_0.0.txt");

    for (int i = 0; i < Nz; i++) {
        density(i) = 0;
        if (i >= 0) {
            {
                for (int j = 0; j < Ns; j++) {
                    density(i) += simp(i) * G2(i, (Ns - 1) - j) * G1(i, j) * ds * dz*A; // Convolute the Greens

                }
            }
        }
        if (isnan((float) density(i)))
            nan_check = true;
        unnorm += simp(i) * density(i) * dz*A; // Compute unnormalized density integral

    }

    db norm = 0, norm1 = 0;
    for (int i = 0; i < Nz; i++) {

        density(i) = (density(i)*(db) Nm * (db) Np) / unnorm;
        if (!nan_check)
            dens << ((db) i * dz) << '\t' << density(i) << endl;
        else {
            cout << "ERROR: Nan in the density. Aborted." << endl;
        }
        norm1 += simp(i) * density(i) * dz;
        norm += simp(i) * density(i) * dz*A;
    }
    cout << "# of monomers =" << norm << endl;
    //cout << "# of monomers(dz) =" << norm1<< endl;

    dens.close();
}

////  SOLVE GS ---------------------------------------------------------------------

void DFT::solveGs() {

    vec r(Nz), a(Nz), b(Nz), c(Nz - 1), Gsol(Nz);
    mat T(Nz, Nz);
    a(0) = 0;
    Zero_vec(Gsol, Nz);
    Zero_vec(r, Nz);

    /// Solving for G1 (Free propagator)  ======================================
    for (int i = 0; i < Nz; i++) {
        if (i > tether && i < Nz - 1)
            Gv1(i) = 1.0; // unity initial conditions
        else
            Gv1(i) = 0.0;

        b(i) = (1.0 + 2.0 * H);
        if (i < Nz - 1) {
            c(i) = -H;
        } else;
        if (i > 0) {
            a(i) = -H;
        } else;

        if (i > 0 && i < Nz - 1 && abs(a(i)) + abs(c(i)) >= abs(b(i))) {
            cout << "a + c > b???" << endl;
        }

    }



    for (int s = 0; s < Ns; s++) {

        for (int i = 0; i < Nz; i++) // Neccessary for the Crank Nicolson scheme
        {
            if (i == 0) {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv1(i) + H * (Gv1(i + 1));
                //r(i)= (1.0 -ds*field(i))*Gv1(i);
            } else if (i > 0 && i < Nz - 1) {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv1(i) + H * (Gv1(i - 1) + Gv1(i + 1));
                //r(i)= (1.0 -ds*field(i))*Gv1(i);
            } else {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv1(i) + H * (Gv1(i - 1));
                //r(i)= (1.0 -ds*field(i))*Gv1(i);
            }
        }



        // Gsol = T.llt().solve(r); //(EIGEN)
        tridag(a, b, c, r, Gsol); // Actually solves the matrix equation.

        Gv1 = Gsol;
        G1.col(s) = Gsol;
    }

    ///   Solving for G2 (Tethered propagator)  ========================================
    for (int i = 0; i < Nz; i++) {
        if (i == tether)
            Gv2(i) = 1.0;
        else
            Gv2(i) = 0.0;
    }


    for (int s = 0; s < Ns; s++) {

        for (int i = 0; i < Nz; i++) // Neccessary for the Crank Nicolson scheme
        {
            if (i == 0) {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv2(i) + H * (Gv2(i + 1));
                //r(i)= (1.0 -ds*field(i))*Gv2(i);
            } else if (i > 0 && i < Nz - 1) {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv2(i) + H * (Gv2(i - 1) + Gv2(i + 1));
                //r(i)= (1.0 -ds*field(i))*Gv2(i);
            } else {
                r(i) = (1.0 - 2.0 * H - ds * field(i)) * Gv2(i) + H * (Gv2(i - 1));
                //r(i)= (1.0 -ds*field(i))*Gv2(i);
            }
        }
        // Gsol = T.llt().solve(r); //(EIGEN)
        tridag(a, b, c, r, Gsol); // Actually solves the matrix equation.
        Gv2 = Gsol;
        // modulate(Gsol,Nz);
        //Gv2=T*r;
        G2.col(s) = Gv2;
    }


    if (deb == 1) {
        for (int i = 0; i < Nz; i++) {
            for (int j = 0; j < Ns; j++) {
                gnug1.send3file((db) (i) * dz, (db) (j) * ds, G1(i, j));
                gnug2.send3file((db) (i) * dz, (db) (j) * ds, G2(i, j));
            }
        }
    }

}


db DFT::simp(int a) // trapezoidal atm
{
    if (a == 0 || a == Nz - 1)
        return 1.0 * 0.5;
    else
        return 1.0;
}