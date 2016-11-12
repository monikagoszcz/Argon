#include "coor.h"

double calcVectorModulus(Coor vector)
{
	return sqrt(vector.x * vector.x +
		vector.y * vector.y +
		vector.z * vector.z);
}

Coor subtractVectors(Coor v1, Coor v2)
{
	return{ v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

Coor addVectors(Coor v1, Coor v2)
{
	return{ v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

Coor calcOppositeVector(Coor v)
{
	return{ -v.x, -v.y, -v.z };
}