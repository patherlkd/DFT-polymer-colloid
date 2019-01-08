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

            cn0(i) += simp(j) * den * w0(rc, (db) i * dz - (db) j * dz) * dz;
            cn1(i) += simp(j) * den * w1(rc, (db) i * dz - (db) j * dz) * dz;
            cn2(i) += simp(j) * den * w2(rc, (db) i * dz - (db) j * dz) * dz;
            cn3(i) += simp(j) * den * w3(rc, (db) i * dz - (db) j * dz) * dz;
            cnv1(i) += simp(j) * den * wv1(rc, (db) i * dz - (db) j * dz) * dz;
            cnv2(i) += simp(j) * den * wv2(rc, (db) i * dz - (db) j * dz) * dz;

        }

    }

}


void DFT::comp_dphi_col1() // d(WB)/dn
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


void DFT::comp_FMT_col1() {

    vec W0(Nz), W1(Nz), W2(Nz), W3(Nz), Wv1(Nz), Wv2(Nz);
    comp_n_col1(); // Compute the weighted densities
    //comp_n_tot(); // Compute the total weighted densities
    //RF(); (comment comp_dphi() out and uncomment this to use the previous Rosenfeld functional.) 
    comp_dphi_col1(); // White bear
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

void DFT::CS(int place) {

    db n; // packing fraction
    db pi = 3.141592653589793238;
    db t1 = 0, t2 = 0, t3 = 0, ts1 = 0, ts2 = 0, ts3 = 0; // split computation into smaller terms
    n = (density(place)*(cube_d(dia)) * pi) / 6.0;
    t1 = (4.0 * n - 3 * sqr_d(n)) / sqr_d(1.0 - n);

    ts1 = (4.0 - 6.0 * n) / sqr_d(1.0 - n);
    ts2 = 2.0 * (4.0 * n - 3.0 * sqr_d(n) / cube_d(1.0 - n));
    ts3 = -2.0 * (-(1.0 / (2 * cube_d(1.0 - n))) + (3.0 / 2.0)*(2.0 - n) / sqr_d(sqr_d(1.0 - n))) * cube_d(1.0 - n)*(1.0 - (1.0 / (db) Nm));
    ts3 *= (1.0 / (2.0 - n));

    t2 = n * (ts1 + ts2 + ts3);

    t3 = -(1.0 - (1.0 / (db) Nm)) * log((2.0 - n) / (2.0 * cube_d(1.0 - n)));

    hs = t1 + t2 + t3;

}

void DFT::RF() {
    db N = (db) Nm;
    for (int i = 0; i < Nz; i++) {
        db a, b, c, d, e, f, g, h, k, z, A, B, loga, logc;
        db sqn2, sqnv2, a4, a3, a2, dno, cun3, cun2, sqn3, r2;
        r2 = r*r;
        a = 1.0 - n3(i);
        a2 = a*a;
        a3 = a * a*a;
        a4 = a2*a2;
        loga = log(a);
        sqn2 = n2(i) * n2(i);
        sqnv2 = nv2(i) * nv2(i);
        sqn3 = n3(i) * n3(i);
        cun2 = n2(i) * n2(i) * n2(i);
        cun3 = n3(i) * n3(i) * n3(i);
        b = 1.0 - frac(sqnv2, sqn2);
        d = frac(1.0 - N, N);
        c = (1.0 / a) + ((n2(i) * b * r) / 2 * a2) + ((sqn2 * b * r2) / 18.0 * a3);
        e = ((sqnv2 * r) / sqn2 * a2);
        e += (b * r / 2.0 * a2);
        e += (sqnv2 * r2) / (9.0 * n2(i) * a3);
        e += (n2(i) * b * r2) / (9.0 * a3);
        logc = log(c);
        f = frac(1.0, a2);
        f += (n2(i) * b * r) / a3;
        f += (sqn2 * b * r2) / 6.0 * a4;
        g = -frac(nv2(i) * r, n2(i) * a2) - frac(nv2(i) * r2, 9.0 * a3);
        h = n3(i)+ (a2 * loga);
        z = n3(i) -(2.0 * a * loga);
        k = n2(i) * sqnv2;
        dno = d * n0(i);
        A = (sqn2 - sqnv2);
        B = (cun2 - 3.0 * n2(i) * sqnv2);

        dphi0(i) = -loga + (d * b * logc);

        dphi1(i) = n2(i) / a;

        dphi2(i) = n1(i) / a;
        dphi2(i) += d * n0(i) * b * (e / c);
        dphi2(i) += (A) / (8.0 * a2 * pi);
        dphi2(i) += 2.0 * d * n0(i) * sqnv2 * (logc / cun2);

        dphi3(i) = (n0(i) / a)+((n1(i) * n2(i) - nv1(i) * nv2(i)) / a2);
        dphi3(i) += d * n0(i) * b * (f / c);
        dphi3(i) += (B * (1.0 / 3.0)) / (4.0 * a3 * pi);

        dphiv1(i) = -(nv2(i) / a);

        dphiv2(i) = -(nv1(i) / a);
        dphiv2(i) += d * n0(i) * b * (g / c);
        dphiv2(i) -= (n2(i) * nv2(i)) / (4.0 * a2 * pi);
        dphiv2(i) -= 2.0 * d * n0(i) * nv2(i)*(logc / sqn2);


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


