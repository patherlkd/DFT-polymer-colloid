/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "colDFT.h"
#include "iostream"

using namespace std;

void DFT::set_potential_mode(unsigned int pm) {
    this->potential_mode = pm;
}

void DFT::set_epp(db epp) {
    this->epp = epp;
}

void DFT::set_lambdapp(db lpp) {
    this->lambdapp = lpp;
}

void DFT::set_lambdapc1(db lpc) {
    this->lambdapc1 = lpc;
}

void DFT::set_epc1(db epc1) {
    this->epc1 = epc1;
}

void DFT::set_ec1c1(db ep){
    this->ec1c1 = ep; 
}

void DFT::set_lambdac1c1(db lc1c1){
    this->lambdac1c1 = lc1c1;
}

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

void DFT::set_wall_strength(db ws) {
    this->wall_strength = ws;
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
    
   // DFT::system_out_file << "Ns: "<<Ns<<"\n";
   // DFT::system_out_file << "ds: "<<ds<<"\n";
}

void DFT::set_dz(unsigned int N, db Z) {

    this->Nz = N;
    this->dz = frac(Z, (db) Nz - 1.0);
    this->h = Z;
}

void DFT::set_A(db a) {
    this->A = a;
  //  DFT::system_out_file << "A= " << A << "\n";
}

void DFT::set_polymers_off(){
    this->polymers_off = true;
}

void DFT::set_topwall_off_c1(){
    this->topwall_off_c1 = true;
}

void DFT::set_gamma(db g) {
    this->gamma = g;
}

void DFT::set_dt(db t) {
    this->dt = t;
}

void DFT::set_DT(db DT) {
    this->DT = DT;
}

void DFT::set_tether(unsigned int t) {
    this->tether = t;
}

void DFT::set_H_solver() {

    this->H = (ds * D) / (2 * dz * dz); // used in numerically solving for the greens function

}

void DFT::set_chem1(db ch) {
    this->chem1 = ch;
}

void DFT::set_ncolloids1(unsigned int nc) {
    this->Nc1 = nc;
}

void DFT::set_rc1(db rc) {
    this->rc1 = rc;
}

void DFT::set_colbulk1(){
    this->colbulk1 = (db)Nc1/(A*h);
}

void DFT::set_col1_init_cut(db c){
    this->col1_init_cut = c;
}

void DFT::set_poly_dens_filename(std::string s) {
    DFT::poly_dens_filename = s;
}

void DFT::set_col1_dens_filename(std::string s) {
    DFT::col1_dens_filename = s;
}

void DFT::set_meanfield_filename(std::string s) {
    DFT::meanfield_filename = s;
}

void DFT::set_external_pot_filename(std::string s) {
    DFT::external_pot_filename = s;
}

void DFT::set_system_out_filename(std::string s) {
    DFT::system_out_filename = s;

    DFT::system_out_file.open(DFT::system_out_filename);

}