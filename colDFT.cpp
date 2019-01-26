/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"

using namespace std;

DFT::DFT(int deb) {

    this->deb = deb;

    gnumf.makefile("NSP1_data/easter_data/TESTmf256_0.0.txt");
    // gnumf.makefile("NSP1_data/easter_data/DFT2048_ideal_mf.txt");
    gnudens.makefile("NSP1_data/easter_data/packing256_0.0.txt"); // Density profile file
    gnucol.makefile("NSP1_data/easter_data/col_dens.txt");
    gnuV.makefile("NSP1_data/Vtest.txt");
    //gnuF.makefile("NSP1_data/easter_data/FE__0.0.txt");// Free energy file
    //gnuLEA.makefile("LEA_term_0.2.txt");


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

    H = (ds * D) / (2 * dz * dz);
    cout << "H: " << H << endl;
    cout << "dz: " << dz << endl;
}

void DFT::evolve() {

    iter = 0;
    conver = 1;
    init_field(0.80); // Initialize mean field
    init_coldensity1(); // Initialize colloid density (bulk?)
    comp_POT(); // Compute external potential. 



    while (conver > gamma) // Do while un converged. Well durr. 
    {
        solveGs(); // Solve for the propagators for polymer density 
        comp_dens(); // Compute density(Z)
        if (conver == 1) {
            comp_n_col1();
            comp_n_pol();
        }
        //  norm();
        update_mf(1); //Update mean field now. Argument is 1 for FMT hard sphere stuff. DO NOT set to 0. Plez.
        update_col1(); // Update the colloid density


        cout << "MF convergence: " << conver << endl;
        cout << "Col density 1 convergence: " << conver_col1 << endl;
        //        cout << "Col density 2 convergence: " << conver_col2 << endl;

        iter++;
    }


}

void DFT::update_col1() {
    ofstream d2("NSP1_data/easter_data/col_dens.txt");
    db ARG = 0, max = 0, old_d = 0, diff = 0;

    comp_FMT_col1();

    for (int i = 0; i < Nz; i++) {
        old_d = coldensity1(i);
        ARG = chem + cc(i) - V(i);
        coldensity1(i) = (1.0 - DT) * coldensity1(i) + dt * colbulk * exp(ARG);
        diff = fabs(old_d - coldensity1(i));
        if (diff > max)
            max = diff;

        gnucol.send2file((db) i*dz, coldensity1(i));
        d2 << (db) i * dz << "\t" << coldensity1(i) << endl;
    }

    conver_col1 = max;
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
        d1(i) = coldensity1(i);
        unnorm += d1(i) * dz*A;
    }
    for (int i = 0; i < Nz; i++) {

        d1(i) = Nc * d1(i) / unnorm;
        coldensity1(i) = d1(i);
        // d2 << (db)i*dz << "\t" << d1(i)<<endl;
        norm += d1(i) * dz*A;
    }

    //   cout << "TOTAL N: " << norm << endl;
    //d2.close();
}

void DFT::init_coldensity1() // Initialize the colloid density
{
    for (int i = 0; i < Nz; i++) {

        if (i <= 20 || i == Nz - 1) {
            coldensity1(i) = 0.0;
        } else if (i > 20)
            coldensity1(i) = colbulk * exp(V(i));
    }
}

void DFT::init_field(db a) {

    for (int i = 0; i < Nz; i++) {
        field(i) = 0.0;
    }
}

void DFT::plotcheck(int op) {
    if (op == 0) {
        gnudens.plot('l');
        //gnumf.plot('p');
    } else if (op == 2) {
        gnuw0.xrange(0.0, 3.0);
        gnuw1.xrange(0.0, 3.0);
        gnuw2.xrange(0.0, 3.0);
        gnuw3.xrange(0.0, 3.0);
        gnuwv1.xrange(0.0, 3.0);
        gnuwv2.xrange(0.0, 3.0);
        gnuw0.plot('l');
        gnuw1.plot('l');
        gnuw2.plot('l');
        gnuw3.plot('l');
        gnuwv1.plot('l');
        gnuwv2.plot('l');
    } else if (op == 3) {
        gnun0.plot('p');
        gnun1.plot('p');
        gnun2.plot('p');
        gnun3.plot('p');
        gnunv1.plot('p');
        gnunv2.plot('p');
    } else if (op == 4) {
        gnup0.plot('p');
        gnup1.plot('p');
        gnup2.plot('p');
        gnup3.plot('p');
        gnupv1.plot('p');
        gnupv2.plot('p');
    } else if (op == 1) {
        gnuc.plot('l');

        gnug1.splot(13), gnug2.splot(11);
    } else
        cout << "ERROR plot debugging op has to = 0 or 1!\n" << endl;

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