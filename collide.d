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

bool gjk(Vector[] p1, Vector[] p2, inout ArraySeq!(Vector) simplex, inout Vector plot,
         inout ArraySeq!(Vector) rb1Simplex, inout ArraySeq!(Vector) rb2Simplex, inout int edgeFlag)
{
	rb1Simplex.append(p1[0]);
	rb2Simplex.append(p2[0]);
    Vector maxD = rb1Simplex.tail() - rb2Simplex.tail();
    simplex.append(maxD);
	float dMag = maxD.magnitude;

	bool collide = false;

	while(true)  // There can be as many as six, three on each polytope.
	{
        simplex.append(support(p1, p2, rb1Simplex, rb2Simplex, maxD.neg()));

        if((maxD*maxD - maxD*simplex.tail()) < EPSILON)
		{
            rb1Simplex.removeTail();    // Remove extra vertex in Simplex
            rb2Simplex.removeTail();
            simplex.removeTail();
            plot = maxD;
            return false;
		}

		if(simplex.size() == 2)  							                        // Line Test
		{
			Vector ab = simplex.get(1) - simplex.get(0);
			Vector origin;
            float t = ((origin-simplex.get(0)) * ab);
            float denom = ab*ab;
            if(t >= denom)
            {
                maxD = simplex.get(1);
                edgeFlag = 1;
            }
                else
            {
                t /= denom;
                maxD = simplex.get(0) + ab*t;
                edgeFlag = 0;
            }
		}
		else { maxD = pointTriangle(simplex, collide, edgeFlag); }	    // Triangle Test

        dMag = maxD.magnitude;
		if(dMag <= 0.3 || collide == true)
		{
		    plot = maxD;
            return true;
		}
	}
	return false;
}

private Vector pointTriangle(inout ArraySeq!(Vector) simplex, inout bool collide, inout int edgeFlag)
{
	Vector origin;
	int i = simplex.size()-1;

	Vector ab = simplex.get(i-1) - simplex.tail();
	Vector ac = simplex.get(i-2) - simplex.tail();;
	Vector ao = origin - simplex.tail();

	float d1 = ab*ao;
	float d2 = ac*ao;
	if(d1 <= 0.0f && d2 <= 0.0f)
	{
	    edgeFlag = i;
	    return simplex.tail();			             // Origin in vertex region outside A
	}

    Vector bo = origin - simplex.get(i-1);           // Origin in vertex region outside B
 	float d3 = ab*bo;
	float d4 = ac*bo;
	float vc = d1*d4 - d3*d2;
	if(vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)		// Origin in edge region outside AB
	{
	    edgeFlag = i-1;
		float v = d1/(d1 - d3);
		return (simplex.tail() + ab*v);
	}

    Vector co = origin - simplex.get(i-2);
    float d5 = ab*co;
    float d6 = ac*co;
	float vb = d5*d2 - d1*d6;
	if(vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)	    // Origin in edge region ouside AC
	{
		float w = d2/(d2-d6);
        edgeFlag = i-2;
		return (simplex.tail() + ac*w);
	}

    collide = true;                                 // Origin inside face region. Penetration!!!

    float va = d3*d6 - d5*d4;
	float denom = 1f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;

	return (simplex.tail() + ab*v + ac*w);
}

private Vector support(Vector[] p1, Vector[] p2, inout ArraySeq!(Vector) rb1Simplex, inout ArraySeq!(Vector) rb2Simplex, Vector maxD)
{
	Vector p1Simp = p1[0];
	Vector p2Simp = p2[0];

	// To find extreme vertexes I employ a brute force method.  This is sufficient for small 2D polygons.
	// Dobkin Kirkpatrick (DK), BSP, or hill climbing algorithms should be investigated for larger polygons and/or 3D.

 	for(int i = 1; i < p1.length; i++)
		if(maxD*p1[i] > maxD*p1Simp) p1Simp = p1[i];

	maxD = maxD*-1;

	for(int i = 1; i < p2.length; i++)
		if(maxD*p2[i] > maxD*p2Simp) p2Simp = p2[i];

    rb1Simplex.append(p1Simp);              // Maintain a list for each polytope
    rb2Simplex.append(p2Simp);

	Vector minkowskiSum = p1Simp - p2Simp;
	return minkowskiSum;
}

void segmentSegment(Vector p1, Vector q1, Vector p2, Vector q2, inout Vector c1, inout Vector c2)
{
    Vector d1 = q1 - p1;
    Vector d2 = q2 - p2;
    Vector r = p1 - p2;
    float a = d1*d1;
    float e = d2*d2;
    float f = d2*r;
    float s, t;

    // Check if either or both segments degenerate into a vertex
    if (a <= EPSILON && e <= EPSILON)                           // Vertex-Vertex
    {
        // Both segments are a vertex
        s = t = 0.0f;
        c1 = p1;
        c2 = p2;
        return;
    }

    if (a <= EPSILON)                                           // Vertex-Edge
    {
        // First segment is a vertex
        s = 0.0f;
        t = f / e; // s = 0 => t = (b*s + f) / e = f / e
        t = clip(t, 0.0f, 1.0f);
    }
    else
    {
        float c = d1*r;
        if (e <= EPSILON)                                       // Vertex-Edge
        {
            // Second segment is a vertex
            t = 0.0f;
            s = clip(-c / a, 0.0f, 1.0f); // t = 0 => s = (b*t - c) / a = -c / a
        }
        else                                                    // Edge-Edge
        {
            // The general case
            float b = d1*d2;
            float denom = a*e-b*b;
            c = d1*r;

            if (denom != 0.0f)                                  // Segments not parallel
            {
                s = clip((b*f - c*e) / denom, 0.0f, 1.0f);
            }
            else                                                // Segments paralllel
            {
                s = 0.0f;
            }

            t = (b*s + f) / e;

            float tnm = b*s + f;
            if (tnm < 0.0f)
            {
                t = 0.0f;
                s = clip(-c / a, 0.0f, 1.0f);
            }
            else if (tnm > e)
            {
                t = 1.0f;
                s = clip((b - c) / a, 0.0f, 1.0f);
            }
            else
            {
                t = tnm / e;
            }
        }
    }

    c1 = p1 + (d1 * s);
    c2 = p2 + (d2 * t);

    return;
}

private float clip(float n, float min, float max)
{
    if (n < min) return min;
    if (n > max) return max;
    return n;
}
