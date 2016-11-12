/*
 * argon.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "argon.h"
#include "coor.h"

#include <cassert>
#include <math.h>
#include <cstdlib>
#include <cstring>

#define PI (3.141592653589793)

using namespace std;

#ifndef _WIN32 /*Enables compilation not only on Windows*/
int fopen_s(FILE **f, const char *name, const char *mode) {
    int ret = 0;
    assert(f);
    *f = fopen(name, mode);
    if (!*f)
        ret = errno;
    return ret;
}
#define fscanf_s fscanf
#endif

/* physical constants */
namespace {

const double kB = 8.31e-3;

}

Parameters getParameters(std::string inputFileName)
{
    Parameters param = {};
	std::ifstream inputFile;
	inputFile.open(inputFileName);
 
		inputFile   >> param.n >> param.m 
					>> param.e >> param.R 
					>> param.f >> param.L 
					>> param.a >> param.T_0 
					>> param.tau >> param.S_0
					>> param.S_d >> param.S_out
					>> param.S_xyz;

	param.N = param.n * param.n * param.n;

    return param;
}

void setInitialState(Parameters &parameters, State &state)
{
    setInitialLocation(parameters, state.Atoms);
    setInitialEnergies(parameters, state.Atoms);
    setInitialMomentum(parameters, state.Atoms);
    setPotentialForcesAndPressure(parameters, state);
}

void setInitialEnergies(Parameters &parameters, vector<Atom> & Atoms)
{
    for(auto & Atom : Atoms)
    {
        Atom.E ={- 0.5 * kB * parameters.T_0 * log(getUniRandom()) ,
                     - 0.5 * kB * parameters.T_0 * log(getUniRandom()),
                     - 0.5 * kB * parameters.T_0 * log(getUniRandom())};
    }
}

void setInitialLocation(Parameters &parameters, vector<Atom> &Atoms)
{
    Coor b[3] = {{parameters.a, 0, 0},
                 {parameters.a / 2.0, parameters.a * sqrt(3.0) / 2.0, 0},
                 {parameters.a / 2.0, parameters.a * sqrt(3.0) / 6.0, parameters.a * sqrt(2.0 / 3.0)}};

    Atoms.resize(static_cast<int>(parameters.N));
    int i = 0;

    for(int i0 = 0; i0 < parameters.n; i0++)
        for(int i1 = 0; i1 < parameters.n; i1++)
            for(int i2 = 0; i2 < parameters.n; i2++, i++)
            {
                Atoms[i].r = {(i0 - (parameters.n-1)/2.) * b[0].x
                              + (i1 - (parameters.n-1)/2.) * b[1].x
                              + (i2 - (parameters.n-1)/2.) * b[2].x,
                              (i0 - (parameters.n-1)/2.) * b[0].y
                              + (i1 - (parameters.n-1)/2.) * b[1].y
                              + (i2 - (parameters.n-1)/2.) * b[2].y,
                              (i0 - (parameters.n-1)/2.) * b[0].z
                              + (i1 - (parameters.n-1)/2.) * b[1].z
                              + (i2 - (parameters.n-1)/2.) * b[2].z};
            }

}

void setInitialMomentum(Parameters &parameters, vector<Atom> &Atoms)
{
    Coor Psum = {0, 0, 0};

	for (auto & Atom : Atoms)
    {
        Atom.p = {
            getPlusOrMinus() * sqrt(2.*parameters.m*Atom.E.x),
            getPlusOrMinus() * sqrt(2.*parameters.m*Atom.E.y),
            getPlusOrMinus() * sqrt(2.*parameters.m*Atom.E.z)};

        Psum.x += Atom.p.x;
        Psum.y += Atom.p.y;
        Psum.z += Atom.p.z;
    }

	for (auto & Atom : Atoms)
    {
        Atom.p = {
            Atom.p.x - Psum.x / parameters.N,
            Atom.p.y - Psum.y / parameters.N,
            Atom.p.z - Psum.z / parameters.N};
    }
}


void setPotentialForcesAndPressure(Parameters &parameters, State &state)
{

    state.V = 0;
    state.T = 0;
    state.P = 0;

	for (auto & Atom : state.Atoms)
    {
        Atom.F = {0,0,0};
    }

    for(int i = 0; i < parameters.N; i++)
    {
        state.V += calcPotentialS(parameters, calcVectorModulus(state.Atoms[i].r)); //potencjal od scianek (10)

        Coor Fis = calcForcesS(parameters, state.Atoms[i].r); //sily odpychania od scianek (14)
        state.Atoms[i].F = addVectors(state.Atoms[i].F, Fis ); //akumulacja do F
        state.P += calcVectorModulus(Fis) / 4. / PI / parameters.L / parameters.L; //akumulacja cisnienia chwilowego (15)

        for(int j = 0; j < i; j++)
        {
            double rij = calcVectorModulus( subtractVectors(state.Atoms[i].r, state.Atoms[j].r));
            state.V += calcPotentialP(parameters, rij); //obliczanie potencjalu par (9) i akumulacja do V
            Coor Fip = calcForcesP(parameters, state.Atoms[i].r, state.Atoms[j].r); //obliczanie sil miedzyAtomowych
			state.Atoms[i].F = addVectors(state.Atoms[i].F, Fip); //  akumulacja do Fi
            state.Atoms[j].F = addVectors(state.Atoms[j].F, calcOppositeVector(Fip)); // akuulacja do Fj
        }

    }
}

double calcPotentialS(Parameters &parameters, double ri)
{
    if(ri <= parameters.L)
        return 0;

    return (0.5 * parameters.f * (ri - parameters.L) * (ri - parameters.L));
}

double calcPotentialP(Parameters &parameters, double rij)
{
    double rRatio = parameters.R / rij;

    return parameters.e * (pow( rRatio , 12 ) - 2 * pow(rRatio, 6) );
}

Coor calcForcesP(Parameters &parameters, Coor ri, Coor rj)
{
    Coor riMinusrj = subtractVectors(ri, rj);
    double rij = calcVectorModulus( riMinusrj );
    double rRatio = parameters.R/rij;

    double A = 12 * parameters.e * (pow( rRatio , 12 ) - pow(rRatio, 6) ) / (rij * rij);

    Coor Fij = {A * riMinusrj.x,
                A * riMinusrj.y,
                A * riMinusrj.z};

    return Fij;
}

Coor calcForcesS(Parameters &parameters, Coor ri)
{
    double mri = calcVectorModulus(ri);
    Coor fS = {0, 0, 0};

    if(mri > parameters.L)
    {
        double A  = parameters.f * (parameters.L - mri)  / mri;
        fS = {ri.x * A, ri.y * A, ri.z * A};
    }

    return fS;
}

void simulate(Parameters &parameters, State &state, ofstream & outputFileXYZ, ofstream & outputFileChar)
{
    double t = 0;

    for(int s = 1; s <= (parameters.S_0 + parameters.S_d); s++)
    {
        t += parameters.tau;
        updateState(parameters, state);

        if(!(s % parameters.S_out))
            outputChar(state, outputFileChar, t);

        if( !(s % parameters.S_xyz))
            outputXYZ(state.Atoms, outputFileXYZ);
    }

}

void updateState(Parameters &parameters, State &state)
{

	for (auto & Atom : state.Atoms)
    {
        Atom.p =
        {
            Atom.p.x + 0.5 * Atom.F.x * parameters.tau,
            Atom.p.y + 0.5 * Atom.F.y * parameters.tau,
            Atom.p.z + 0.5 * Atom.F.z * parameters.tau
        };

        Atom.r =
        {
            Atom.r.x + Atom.p.x * parameters.tau / parameters.m,
            Atom.r.y + Atom.p.y * parameters.tau / parameters.m,
            Atom.r.z + Atom.p.z * parameters.tau / parameters.m
        };
    }

    setPotentialForcesAndPressure(parameters, state);

	for (auto & Atom : state.Atoms)
    {
        Atom.p =
        {
            Atom.p.x  + 0.5 * Atom.F.x * parameters.tau,
            Atom.p.y  + 0.5 * Atom.F.y * parameters.tau,
            Atom.p.z  + 0.5 * Atom.F.z * parameters.tau
        };
    }

    setEnergyAndTemperature(parameters, state);
}

void setEnergyAndTemperature(Parameters &parameters, State &state)
{
    state.H = 0;

    for(int i = 0; i < parameters.N; i++)
    {
        state.H += calcVectorModulus(state.Atoms[i].p) *
                calcVectorModulus(state.Atoms[i].p) / 2. / parameters.m;
    }

    state.T = state.H * 2. / 3. / parameters.N / kB;
    state.H += state.V;
}

double getUniRandom()
{
    return (double)rand() / (RAND_MAX);
}

int getPlusOrMinus()
{
    return (rand() % 2) ? 1 : -1;
}

void outputXYZ(vector<Atom> & Atoms, ofstream &outputFile)
{
    outputFile << Atoms.size() << endl << "comment"<< endl;

	for (auto & Atom : Atoms)
    {
        outputFile << "Ar " << Atom.r.x << " " << Atom.r.y << " " << Atom.r.z << endl;
    }
}

void outputMomentum(vector<Atom> &Atoms)
{
    ofstream outputFile;
    outputFile.open("outputMomentum.txt");

	for (auto & Atom : Atoms)
    {
        outputFile << Atom.p.x << "\t" << Atom.p.y << "\t"
                   << Atom.p.z << endl;
    }

    outputFile.close();
}

void outputChar(State & state, ofstream &outputFile, double &t)
{
    outputFile << t << "\t" << state.H << "\t" << state.V << "\t"
               << state.T << "\t" << state.P << endl;
}

