/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "useful.h"

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

        if (deb == 2) {
            gnuw0.send2file((db) i*dz, w0(r, (db) i * dz));
            gnuw1.send2file((db) i*dz, w1(r, (db) i * dz));
            gnuw2.send2file((db) i*dz, w2(r, (db) i * dz));
            gnuw3.send2file((db) i*dz, w3(r, (db) i * dz));
            gnuwv1.send2file((db) i*dz, wv1(r, (db) i * dz));
            gnuwv2.send2file((db) i*dz, wv2(r, (db) i * dz));
        }

        if (deb == 3) {
            gnun0.send2file((i)*(dz), n0(i));
            gnun1.send2file((i)*(dz), n1(i));
            gnun2.send2file((i)*(dz), n2(i));
            gnun3.send2file((i)*(dz), n3(i));
            gnunv1.send2file((i)*(dz), nv1(i));
            gnunv2.send2file((i)*(dz), nv2(i));
        }

    }

}

void DFT::comp_dphi_pol() //d(WB)/dn
{
    db N = (db) Nm;

    db a, b, c, d, e;
    db A, B, C, D;
    db loga, logc;
    db N0, N1, N2, N3, Nv1, Nv2;
    db R;

    R = frac(r*rc, rc + r);

    for (int i = 0; i < Nz; i++) {

        N0 = cn0(i) + n0(i);
        N1 = cn1(i) + n1(i);
        N2 = cn2(i) + n2(i);
        N3 = cn3(i) + n3(i);
        Nv1 = cnv1(i) + nv1(i);
        Nv2 = cnv2(i) + nv2(i);

        a = 1.0 - N3;
        loga = log(a);
        b = 1.0 - frac(sq(nv2(i)), sq(n2(i)));
        B = 1.0 - frac(sq(Nv2), sq(N2));
        c = recip(a) + frac(2.0 * N2 * (1 - sq(Nv2) * sq(R)), 9.0 * cu(a)) + frac(N2 * B*R, sq(a));
        logc = log(c);
        d = frac(1.0 - (db) Nm, (db) Nm);
        e = recip(sq(a)) + frac(2.0 * N2 * (1 - sq(Nv2) * sq(R)), 3.0 * sq(sq(a))) + frac(2.0 * N2 * B*R, cu(a));

        dphi0(i) = -loga + d * b*logc;

        dphi1(i) = frac(N2, a);

        dphi2(i) = frac(N1, a);
        dphi2(i) += n0(i) * d * frac(b * (frac(2.0 * (1 - sq(Nv2)) * R, 9.0 * cu(a)) + frac(2.0 * sq(Nv2) * R, sq(N2) * sq(a)) + frac(B*R, sq(a))), c);
        dphi2(i) += frac((2.0 * N2 - 3.0 * sq(Nv2))*(N3 + sq(a) * loga), 36.0 * sq(a) * N3 * pi);
        dphi2(i) += frac(2.0 * d * n0(i) * sq(nv2(i)) * logc, cu(n2(i)));


        dphi3(i) = frac(n0(i), a);
        dphi3(i) += frac(N1 * N2 - Nv1*Nv2, sq(a));
        dphi3(i) += d * n0(i) * frac(b*e, c);
        dphi3(i) += frac((sq(N2) - 3.0 * N2 * sq(Nv2))*(N3 - 2 * a * loga), 36.0 * sq(a) * N3 * pi);
        dphi3(i) -= frac((sq(N2) - 3.0 * N2 * sq(Nv2))*(N3 + 2 * sq(a) * loga), 36.0 * sq(a) * sq(N3) * pi);
        dphi3(i) += frac((sq(N2) - 3.0 * N2 * sq(Nv2))*(N3 + 2 * sq(a) * loga), 18.0 * cu(a)*(N3) * pi);


        dphiv1(i) = -frac(Nv2, a);

        dphiv2(i) = -frac(Nv1, a) + d * n0(i) * frac(b * (-(frac(4.0 * N2 * Nv2 * sq(R), 9.0 * cu(a)) + frac(2.0 * Nv2*R, N2 * sq(a)))), c);
        dphiv2(i) -= frac(N2 * Nv2 * (N3 + sq(a) * loga), 6.0 * sq(a) * N3 * pi);
        dphiv2(i) -= frac(2.0 * d * n0(i) * nv2(i) * logc, sq(n2(i)));


        if (deb == 4) {
            gnup0.send2file((i)*(dz), dphi0(i));
            gnup1.send2file((i)*(dz), dphi1(i));
            gnup2.send2file((i)*(dz), dphi2(i));
            gnup3.send2file((i)*(dz), dphi3(i));
            gnupv1.send2file((i)*(dz), dphiv1(i));
            gnupv2.send2file((i)*(dz), dphiv2(i));
        }
    }
}

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


void DFT::comp_n_col() {

    for (int i = 0; i < Nz; i++) {
        static db den = 0;
        cn0(i) = 0;
        cn1(i) = 0;
        cn2(i) = 0;
        cn3(i) = 0;
        cnv1(i) = 0;
        cnv2(i) = 0;
        for (int j = 0; j < Nz; j++) {
            den = coldensity(j);

            cn0(i) += simp(j) * den * w0(rc, (db) i * dz - (db) j * dz) * dz;
            cn1(i) += simp(j) * den * w1(rc, (db) i * dz - (db) j * dz) * dz;
            cn2(i) += simp(j) * den * w2(rc, (db) i * dz - (db) j * dz) * dz;
            cn3(i) += simp(j) * den * w3(rc, (db) i * dz - (db) j * dz) * dz;
            cnv1(i) += simp(j) * den * wv1(rc, (db) i * dz - (db) j * dz) * dz;
            cnv2(i) += simp(j) * den * wv2(rc, (db) i * dz - (db) j * dz) * dz;

        }

    }

}


void DFT::comp_dphi_col() // d(WB)/dn
{
    db N = (db) Nm;

    db a, b, c, d, e;
    db A, B, C, D;
    db loga, logc;
    db N0, N1, N2, N3, Nv1, Nv2;
    db R;

    R = frac(r*rc, rc + r);

    for (int i = 0; i < Nz; i++) {
        N0 = cn0(i) + n0(i);
        N1 = cn1(i) + n1(i);
        N2 = cn2(i) + n2(i);
        N3 = cn3(i) + n3(i);
        Nv1 = cnv1(i) + nv1(i);
        Nv2 = cnv2(i) + nv2(i);

        a = 1.0 - N3;
        loga = log(a);
        b = 1.0 - frac(sq(nv2(i)), sq(n2(i)));
        B = 1.0 - frac(sq(Nv2), sq(N2));
        c = recip(a) + frac(2.0 * N2 * (1 - sq(Nv2) * sq(R)), 9.0 * cu(a)) + frac(N2 * B*R, sq(a));
        logc = log(c);
        d = frac(1.0 - (db) Nm, (db) Nm);
        e = recip(sq(a)) + frac(2.0 * N2 * (1 - sq(Nv2) * sq(R)), 3.0 * sq(sq(a))) + frac(2.0 * N2 * B*R, cu(a));

        cdphi0(i) = -loga;

        cdphi1(i) = frac(N2, a);

        cdphi2(i) = dphi0(i) - frac(2.0 * d * n0(i) * sq(nv2(i)) * logc, cu(n2(i)));

        cdphi3(i) = dphi3(i);

        cdphiv1(i) = dphiv1(i);

        cdphiv2(i) = dphiv2(i) + frac(2.0 * d * n0(i) * nv2(i) * logc, sq(n2(i)));


    }
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





