/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"


using namespace std;

db DFT::w0(db R, db x) {
    db w = 0;

    w = w2(R, x) / (R * R * pi * 4.0);
    return w;
}

db DFT::w1(db R, db x) {
    db w = 0;
    w = w2(R, x) / (R * pi * 4.0);
    return w;
}

db DFT::w2(db R, db x) {
    db w = 0;
    w = correct(R, x)*2.0 * (pi * R) * heaviside(R, x);

    return w;
}

db DFT::w3(db R, db x) {
    db w = 0;
    w = correct(R, x) * pi * (R * R - x * x) * heaviside(R, x);
    return w;
}

db DFT::wv1(db R, db x) {
    db w = 0;
    w = wv2(R, x) / (4.0 * pi * R);
    return w;
}

db DFT::wv2(db R, db x) {
    db w = 0;
    w = correct(R, x)*2.0 * pi * x * heaviside(R, x);
    return w;
}


