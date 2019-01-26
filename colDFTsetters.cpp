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
    this->r = dia * 0.5;
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


void DFT::set_poly_dens_filename(std::string s){
    DFT::poly_dens_filename=s;
}
void DFT::set_col1_dens_filename(std::string s){
    DFT::col1_dens_filename=s;
}
void DFT::set_meanfield_filename(std::string s){
    DFT::meanfield_filename=s;
}
void DFT::set_external_pot_filename(std::string s){
    DFT::external_pot_filename=s;
}
void DFT::set_system_out_filename(std::string s){
    DFT::system_out_filename=s;
}