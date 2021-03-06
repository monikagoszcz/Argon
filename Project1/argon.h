/*
 * argon.h
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#ifndef _argon_h_
#define _argon_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "coor.h"

struct Atom
{
    Coor r;
    Coor p;
    Coor E;
    Coor F;

};

struct Parameters
{
    double n;
    double m;
    double e;
    double R;
    double f;
    double L;
    double a;
    double T_0;
    double tau;
    double S_0;
    double S_d;
    int S_out;
    int S_xyz;
	double N;
};

struct State
{
    double V;
    double P;
    double H;
    double T;
    std::vector<Atom> Atoms;
};

Parameters getParameters(std::string inputFileName);

void setInitialState(const Parameters &Parameters, State & state);
void setInitialLocation(const Parameters &Parameters, std::vector<Atom> &Atoms);
void setInitialEnergies(const Parameters &Parameters, std::vector<Atom> &Atoms);
void setInitialMomentum(const Parameters &Parameters, std::vector<Atom> &Atoms);

void setPotentialForcesAndPressure(const Parameters &Parameters, State & state);

double calcPotentialS(const Parameters &Parameters, double ri);
double calcPotentialP(const Parameters &parameters, double rRatioPow6);
Coor calcForcesP(const Parameters &parameters, Coor ri, Coor rj, double rRatioPow6);
Coor calcForcesS(const Parameters &Parameters, Coor ri, double mri);

void simulate(const Parameters &Parameters, State &state, std::ofstream &outputFileXYZ, std::ofstream &outputFileChar);

void updateState(const Parameters &Parameters, State & state);
void setEnergyAndTemperature(const Parameters &Parameters, State & state);

double getUniRandom();
int getPlusOrMinus();

void outputXYZ(std::vector<Atom> &Atoms, std::ofstream &outputFile);
void outputChar(State &state, std::ofstream &outputFile, double &t);
void outputMomentum(std::vector<Atom> &Atoms);

#endif 

