/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"

using namespace std;


DFT::DFT(int d) {


    gnumf.makefile("NSP1_data/easter_data/TESTmf256_0.0.txt");
    // gnumf.makefile("NSP1_data/easter_data/DFT2048_ideal_mf.txt");
    gnudens.makefile("NSP1_data/easter_data/packing256_0.0.txt"); // Density profile file
    gnucol.makefile("NSP1_data/easter_data/col_dens.txt");
    gnuV.makefile("NSP1_data/Vtest.txt");
    //gnuF.makefile("NSP1_data/easter_data/FE__0.0.txt");// Free energy file
    //gnuLEA.makefile("LEA_term_0.2.txt");
    deb = d;
    if (deb == 1) {
        gnuc.makefile("FMT_field_contribution.txt");
        gnug1.makefile("uniform_greens_function.txt");
        gnug2.makefile("tethered_greens_function.txt");
    } else if (deb == 2) {
        gnuw0.makefile("w0.txt");
        gnuw1.makefile("w1.txt");
        gnuw2.makefile("w2.txt");
        gnuw3.makefile("w3.txt");
        gnuwv1.makefile("wv1.txt");
        gnuwv2.makefile("wv2.txt");
    } else if (deb == 3) {
        gnun0.makefile("n0.txt");
        gnun1.makefile("n1.txt");
        gnun2.makefile("n2.txt");
        gnun3.makefile("n3.txt");
        gnunv1.makefile("nv1.txt");
        gnunv2.makefile("nv2.txt");
    } else if (deb == 4) {
        gnup0.makefile("p0.txt");
        gnup1.makefile("p1.txt");
        gnup2.makefile("p2.txt");
        gnup3.makefile("p3.txt");
        gnupv1.makefile("pv1.txt");
        gnupv2.makefile("pv2.txt");
    }

    Nc = 1000;

    h = 100.0;
    e = 20.0;
    set_b(1.0);
    set_dia(1.0);
    set_Np(25);
    set_Nm(300);
    set_D(0.16666 * sqr_d(b));
    set_ds(8000);
    set_dz(512, h);
    set_A(757.575757576);
    // set_A(500.0);
    set_gamma(0.0001);
    set_dt(0.01);
    DT = 0.07;

    conv_fact = 0.26424111765;
    r = 0.5 * dia;
    rc = r + r;
    chem = 5;
    colbulk = (db) Nc / (A * h);
    sigma = 0.76;
    cut = 2 * dia;
    eps = 0.0;
    lambda = 0.76;



    graft = ((db) Np) / A;
    tether = 4;

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
    coldensity.resize(Nz);
    G1.resize(Nz, Ns);
    G2.resize(Nz, Ns);

    Zero_vec(Gv1, Nz);
    Zero_vec(Gv2, Nz);

    Zero_vec(field, Nz);
    Zero_vec(density, Nz);
    Zero_vec(coldensity, Nz);



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

    H = (ds * D) / (2 * dz * dz);
   cout << "H: " << H << endl;
    cout << "dz: " << dz << endl;
}


void DFT::evolve() {

    iter = 0;
    conver = 1;
    init_field(0.80); // Initialize mean field
    init_coldensity(); // Initialize colloid density (bulk?)
    comp_POT(); // Compute external potential. 



    while (conver > gamma) // Do while un converged. Well durr. 
    {
        solveGs(); // Solve for the propagators
        comp_dens(); // Compute density(Z)
        if (conver == 1) {
            comp_n_col();
            comp_n_pol();
        }
        //  norm();
        update_mf(1); //Update mean field now. Argument is 1 for FMT hard sphere stuff. DO NOT set to 0. Plez.
        update_col(); // Update the colloid density

        if (iter % 1 == 0) {
            cout << "MF convergence: " << conver << endl;
            cout << "Col density convergence: " << conver_col << endl;
        }


        //     gnudens.plot('l'); // This plots the density real time (see gnuplot.h (A useful creation I made.. but needs improving) for how) 
        //     gnucol.plot('p');
        iter++;
    }


}

void DFT::update_col() {
    ofstream d2("NSP1_data/easter_data/col_dens.txt");
    db ARG = 0, max = 0, old_d = 0, diff = 0;

    comp_FMT_col();

    for (int i = 0; i < Nz; i++) {
        old_d = coldensity(i);
        ARG = chem + cc(i) - V(i);
        coldensity(i) = (1.0 - DT) * coldensity(i) + dt * colbulk * exp(ARG);
        diff = fabs(old_d - coldensity(i));
        if (diff > max)
            max = diff;

        gnucol.send2file((db) i*dz, coldensity(i));
        d2 << (db) i * dz << "\t" << coldensity(i) << endl;
    }

    conver_col = max;
    d2.close();
}

void DFT::update_mf(int ch) {
    hs = 0;

    att = 0;
    db max = 0, ave = 0;
    db old_mf = 0, diff = 0;
    if (ch == 1) {
        comp_FMT_pol();
    }

    for (int i = 0; i < Nz; i++) {
        old_mf = field(i);
        if (ch == 0)
            CS(i);

        field(i) = old_mf + dt * (-old_mf + c(i) + V(i));

        gnumf.send2file((db) i*dz, field(i));
        if (deb == 1)
            gnuc.send2file((db) i * dz, c(i));
        diff = fabs(old_mf - field(i));
        if (diff > max)
            max = diff;


        //    ave+=fabs(old_mf-field(i));   
    }

    conver = max;
}

void DFT::norm() {
    //ofstream d2("NSP1_data/easter_data/col_dens.txt");
    vec d1;
    d1.resize(Nz);
    db unnorm = 0, norm = 0;
    for (int i = 0; i < Nz; i++) {
        d1(i) = coldensity(i);
        unnorm += d1(i) * dz*A;
    }
    for (int i = 0; i < Nz; i++) {

        d1(i) = Nc * d1(i) / unnorm;
        coldensity(i) = d1(i);
        // d2 << (db)i*dz << "\t" << d1(i)<<endl;
        norm += d1(i) * dz*A;
    }

    cout << "TOTAL N: " << norm << endl;
    //d2.close();
}

void DFT::init_coldensity() // Initialize the colloid density
{
    for (int i = 0; i < Nz; i++) {

        if (i <= 20 || i == Nz - 1) {
            coldensity(i) = 0.0;
        } else if (i > 20)
            coldensity(i) = colbulk * exp(V(i));
    }
}

void DFT::init_field(db a) {
    
    for(int i = 0;i < Nz; i++){
        field(i) = 0.0;
    }
    
  /*  ifstream imf("NSP1_data/easter_data/DFT512_ideal_mf.txt");
    rini(57389);
    db Z_, F_;
    string line;
    /*for(int i=0;i < field.rows(); i++)
      {
        if( i*dz >=0  && i*dz <(Nz)*dz)
          field(i) = a*exp(-(db)i*dz);
          }
    int i = 0;
    while (getline(imf, line)) {

        stringstream ss(line);

        ss >> Z_ >> F_;
        if (i < Nz)
            field(i) = F_;
        i++;
    }
    imf.close();*/
}