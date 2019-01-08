/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "iostream"

using namespace std;

void DFT::set_dia(db dia) {
    this->dia = dia;
    this->r = dia*0.5;
}

void DFT::set_Np(unsigned int N) {
    this->Np = N;
}

void DFT::set_Nm(unsigned int N) {
    this->Nm = N;
}

void DFT::set_D(db coef) {
    this->D = coef;
}

void DFT::set_b(db B) {
    this->b = B;
}

void DFT::set_ds(unsigned int N) {
    this->Ns = N;
    this->ds = ((db) Nm / (db) Ns);
    cout << "ds" << ds << endl;
}

void DFT::set_dz(unsigned int N, db Z) {
    Nz = N;
    dz = frac(Z, (db) Nz - 1);
}

void DFT::set_A(db a) {
    A = a;
}

void DFT::set_gamma(db g) {
    gamma = g;
}

void DFT::set_dt(db t) {
    dt = t;
}