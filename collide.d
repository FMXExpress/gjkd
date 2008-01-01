/* gjkD - An implementation of the Gilbert-Johnson-Keerthi algorithm 
 * for the collision detection of convex objects, written in D.
 * Copyright (C) 2007-2008 Mason A. Green (mason.green@gmail.com)
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

import math;
import tango.io.Stdout;

bool gjk(Vector[] p1, Vector[] p2, inout Vector[] Bsimplex, inout Vector plot, inout int simpIndex)
{
					
	Vector[] simplex;
	Vector origin, a, d;

	bool collide = false;
	int i,b = 0;
	float dMag; 

	simplex.length = p1.length*p2.length;
	simplex[0] = d = p1[0] - p2[2];
	dMag = d.magnitude;

	foreach(inout Vector v; Bsimplex) v = d;
	
	while(true) 
	{
		a = support(p1, p2, d.neg());
		plot = d;
					
		if((dMag*dMag - a*d) < 1) break;
			
		simplex[i+1] = Bsimplex[i+1] = a;

		if(i == 0)  							// Line Test
		{
				
			Vector ab = simplex[1] - simplex[0]; 
			float t = ((origin-simplex[0]) * ab);
			float denom = ab * ab;
				
			if(t >= denom)
			{
				d = simplex[1];
				b = 1;
			}
			else
			{
				t /= denom;
				d = simplex[0] + ab*t;
			}
				
		}
		else	d = closestPointTriangle(simplex, i+1, collide, b); 	// Triangle Test
		
		dMag = d.magnitude;
		plot = d;
		simpIndex = ++i;
		if((dMag <= 0.1) || (collide == true) ) return true; 	// Rounding tolerance
		
	} 

	return false;
}

Vector closestPointTriangle(Vector[] simplex, int i, inout bool collide, inout int b)
{
	Vector origin;
	Vector ab = simplex[i-1] - simplex[i];
	Vector ac = simplex[i-2] - simplex[i];
	Vector ao = origin - simplex[i];

	float d1 = ab*ao;
	float d2 = ac*ao;
	if(d1 <= 0 && d2 <= 0)
	{
		b = i;
		return simplex[i]; 			// Origin in vertex region outside A
	}

	Vector bo = origin - simplex[i-1];
	float d3 = ab*bo;
	float d4 = ac*bo;
	
	float vc = d1*d4 - d3*d2;
	if(vc <= 0 && d1 >= 0 && d3 <= 0)		// Origin in edge region outside AB
	{
		float v = d1/(d1 - d3);
		b = i-1;
		return (simplex[i] + ab*v);
	}

	Vector co = origin - simplex[i-2];
	float d5 = ab * co;
	float d6 = ac * co;
	float vb = d5*d2 - d1*d6;
	if(vb <= 0 && d2 >= 0 && d6 <= 0)		// Origin in edge region ouside AC	
	{
		float w = d2/(d2-d6);
		b = i-2;
		return simplex[i] + ac*w;
	}
	
	collide = true;

	float va = d3*d6 - d5*d4;			// Origin inside face region. Collision!!!
	float denom = 1 / (va + vb + vc);	
	float v = vb * denom;
	float w = vc * denom;

	return simplex[i] + ab*v + ac*w;
		
}

Vector support(Vector[] p1, Vector[] p2, Vector d)
{
		
	Vector minkowskiSum;
	Vector p1Simp = p1[0];
	Vector p2Simp = p2[0];

	d.normalize();

	// To find extreme vertexes I employ a brute force method.  This is sufficient for small 2D polygons.
	// DK, BSP, or hill climbing should be investigated for larger polygons and/or 3D. 

	for(int i = 1; i < p1.length; i++)		
		if(d*p1[i] > d*p1Simp) p1Simp = p1[i];

	d = d*-1;

	for(int i = 1; i < p2.length; i++)		
		if(d*p2[i] > d*p2Simp) p2Simp = p2[i];

	minkowskiSum = p1Simp - p2Simp;
	return minkowskiSum;
		
}


