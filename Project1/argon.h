/*
 * argon.h
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#ifndef ARGON_H_
#define ARGON_H_


#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

struct coor
{
	double x;
	double y;
	double z;
};

struct atom
{
	coor r;
	coor p;
	coor E;
	coor F;
	
};

struct parameters
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
};

struct state
{
	double V;
	double P;
	double H;
	double T;
	vector<atom> atoms; 
};

parameters getParameters(string inputFileName);

void setInitialState(parameters &parameters, state & state);
void setInitialLocation(parameters &parameters, vector<atom> & atoms);
void setInitialEnergies(parameters &parameters, vector<atom> & atoms);
void setInitialMomentum(parameters &parameters, vector<atom> & atoms);
void setPotentialForcesAndPressure(parameters &parameters, state & state);
double calcPotentialS(parameters &parameters, double ri);
double calcPotentialP(parameters &parameters, double rij);
coor calcForcesP(parameters &parameters, coor ri, coor rj);
coor calcForcesS(parameters &parameters, coor ri);

void simulate(parameters &parameters, state & state, ofstream & outputFileXYZ, ofstream & outputFileChar);

void updateState(parameters &parameters, state & state);
void setEnergyAndTemperature(parameters &parameters, state & state);

double getUniRandom();
int getPlusOrMinus();
double calcVectorModulus(coor vector);
coor subtractVectors(coor v1, coor v2);
coor addVectors(coor v1, coor v2);
coor calcOppositeVector(coor v);

void outputXYZ(vector<atom> & atoms, ofstream & outputFile);
void outputChar(state & state, ofstream & outputFile, double & t);
void outputMomentum(vector<atom> & atoms);

#endif /* ARGON_H_ */
