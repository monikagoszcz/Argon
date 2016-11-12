/*
 * main.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "argon.h"
#include "coor.h"

int main(int /* argc */, char* /* argv */[])
{
    State state;
    Parameters parameters = getParameters("data.txt");

    setInitialState(parameters, state);

    std::ofstream outputFileXYZ;
    outputFileXYZ.open("output.xyz");

	std::ofstream outputFileChar;
    outputFileChar.open("output.txt");

    outputFileChar << "t" << "\t" << "H" << "\t" << "V" << "\t" << "T" << "\t " << "P" << std::endl;

    outputXYZ(state.Atoms, outputFileXYZ);

    outputMomentum(state.Atoms);

    simulate(parameters, state, outputFileXYZ, outputFileChar);

    return 0;
}


