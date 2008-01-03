/* gjkD - An implementation of the Gilbert-Johnson-Keerthi algorithm
 * for the collision detection of convex objects, written in D.
 * Copyright (C) 2007-2008 Mason A. Green
 *
 * This file is part of gjkD
 *
 * gjkD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gjkD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
module collision;

import tango.util.collection.ArraySeq;

import math;

bool gjk(Vector[] p1, Vector[] p2, inout Vector[] Bsimplex, inout Vector plot, inout int simpIndex)
{
	ArraySeq!(Vector) simplex;
	simplex = new ArraySeq!(Vector);

    Vector maxD = p1[0] - p2[2];
    simplex.append(maxD);

	float dMag = maxD.magnitude;
    plot = maxD;

	foreach(inout Vector v; Bsimplex) v = maxD;

	bool collide = false;
	int i = 0;

	while(i < 6)  // There can be as many as six, three on each polytope.
	{
        simplex.append(support(p1, p2, maxD.neg()));

        Bsimplex[i+1] = simplex.tail();

		if((dMag*dMag - maxD*simplex.tail()) < 1) return false;

		if(i == 0)  							                        // Line Test
		{
			Vector ab = simplex.get(1) - simplex.get(0);
			Vector origin;
            float t = ((origin-simplex.get(0)) * ab);
			float denom = ab * ab;
			if(t >= denom)
			{
				maxD = simplex.get(1);
			}
			else
			{
				t /= denom;
				maxD = simplex.get(0) + ab*t;
			}
		}
		else maxD = closestPointTriangle(simplex, i+1, collide); 	    // Triangle Test

		dMag = maxD.magnitude;                                          // Seperation distance between polytopes
		plot = maxD;
		simpIndex = ++i;
		if((dMag <= 0.3) || (collide == true) ) return true;
	}
	return false;
}

Vector closestPointTriangle(inout ArraySeq!(Vector) simplex, int i, inout bool collide)
{
	Vector origin;
	Vector ab = simplex.get(i-1) - simplex.get(i);
	Vector ac = simplex.get(i-2) - simplex.get(i);
	Vector ao = origin - simplex.get(i);

	float d1 = ab*ao;
	float d2 = ac*ao;
	if(d1 <= 0 && d2 <= 0)
	{
		return simplex.get(i); 			    // Origin in vertex region outside A
	}

	Vector bo = origin - simplex.get(i-1);
	float d3 = ab*bo;
	float d4 = ac*bo;

	float vc = d1*d4 - d3*d2;
	if(vc <= 0 && d1 >= 0 && d3 <= 0)		// Origin in edge region outside AB
	{
		float v = d1/(d1 - d3);
		return (simplex.get(i) + ab*v);
	}

	Vector co = origin - simplex.get(i-2);
	float d5 = ab * co;
	float d6 = ac * co;
	float vb = d5*d2 - d1*d6;
	if(vb <= 0 && d2 >= 0 && d6 <= 0)		// Origin in edge region ouside AC
	{
		float w = d2/(d2-d6);
		return (simplex.get(i) + ac*w);
	}

    collide = true;                         // Origin inside face region. Collision!!!

    float va = d3*d6 - d5*d4;
	float denom = 1 / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;

	return (simplex.get(i) + ab*v + ac*w);
}

Vector support(Vector[] p1, Vector[] p2, Vector maxD)
{
	Vector minkowskiSum;
	Vector p1Simp = p1[0];
	Vector p2Simp = p2[0];

	maxD.normalize();

	// To find extreme vertexes I employ a brute force method.  This is sufficient for small 2D polygons.
	// Dobkin Kirkpatrick (DK), BSP, or hill climbing algorithms should be investigated for larger polygons and/or 3D.

	for(int i = 1; i < p1.length; i++)
		if(maxD*p1[i] > maxD*p1Simp) p1Simp = p1[i];

	maxD = maxD*-1;

	for(int i = 1; i < p2.length; i++)
		if(maxD*p2[i] > maxD*p2Simp) p2Simp = p2[i];

	minkowskiSum = p1Simp - p2Simp;
	return minkowskiSum;
}
