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
 * along with gjkD.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
module collision;

import tango.util.collection.ArraySeq;
import tango.math.Math;

import math;

// The Gilbert-Johnson-Keerthi algorithm
bool gjk(Vector[] p1, Vector[] p2, inout ArraySeq!(Vector) sAB, inout ArraySeq!(Vector) sA,
         inout ArraySeq!(Vector) sB, inout Entry e)
{
    sA.append(p1[0]);
    sB.append(p2[0]);
    sAB.append(sA.tail() - sB.tail());

    e = constructEntry(sAB.tail, sAB.tail, sA.tail, sB.tail, sA.tail, sB.tail);

    bool penetrate;
    int failsafe;

    while (failsafe++ < 10) // Don't want to get caught in an infinite loop!
    {
        sAB.append(support(p1, p2, sA, sB, e.v.neg));

        if ((e.v*e.v - e.v*sAB.tail()) <= EPSILON) return false;

        // Line Test
        if (sAB.size == 2)
        {
            e = constructEntry(sAB.get(0), sAB.get(1), sA.get(0), sB.get(0), sA.get(1), sB.get(1));
            if(e.key <= EPSILON) penetrate = true;
        }
        // Triangle Test
        else e = pointTriangle(sAB, sA, sB, penetrate);

        if (penetrate == true) return true;
    }
    return false;
}

private Entry pointTriangle(ArraySeq!(Vector) sAB, ArraySeq!(Vector) sA, ArraySeq!(Vector) sB, inout bool penetrate)
{
    int i = sAB.size()-1;
    Vector ab = sAB.get(i-1) - sAB.tail;
    Vector ac = sAB.get(i-2) - sAB.tail;
    Vector ao = sAB.tail.neg;

    double d1 = ab*ao;
    double d2 = ac*ao;
    // Origin in vertex region outside A
    if (d1 <= 0.0f && d2 <= 0.0f)
    {
        return constructEntry(sAB.get(i), sAB.tail, sA.get(i), sB.get(i), sA.tail,sB.tail);
    }

    Vector bo = sAB.get(i-1).neg;
    double d3 = ab*bo;
    double d4 = ac*bo;
    double vc = d1*d4 - d3*d2;
    // Origin in edge region outside AB
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
    {
        return constructEntry(sAB.get(i-1), sAB.tail, sA.get(i-1), sB.get(i-1), sA.tail, sB.tail);
    }

    Vector co = sAB.get(i-2).neg;
    double d5 = ab*co;
    double d6 = ac*co;
    double vb = d5*d2 - d1*d6;
    // Origin in edge region ouside AC
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
    {
        return constructEntry(sAB.get(i-2), sAB.tail, sA.get(i-2), sB.get(i-2), sA.tail, sB.tail);
    }
    // Origin inside face region. Penetration!!!
    penetrate = true;
    Entry fooBar;
    return fooBar;
}

private Vector support(Vector[] p1, Vector[] p2, inout ArraySeq!(Vector) sA, inout ArraySeq!(Vector) sB, Vector maxD)
{
    Vector p1Simp = p1[0];
    Vector p2Simp = p2[0];

    // To find extreme vertexes I employ a brute force method.  This is sufficient for small 2D polygons.
    // Dobkin Kirkpatrick (DK), BSP, or hill climbing algorithms should be investigated for larger polygons and/or 3D.

    for (int i = 1; i < p1.length; i++)
        if (maxD*p1[i] > maxD*p1Simp) p1Simp = p1[i];

    maxD = maxD*-1;

    for (int i = 1; i < p2.length; i++)
        if (maxD*p2[i] > maxD*p2Simp) p2Simp = p2[i];

    sA.append(p1Simp);              // Maintain a list for each polytope
    sB.append(p2Simp);

    Vector minkowskiSum = p1Simp - p2Simp;
    return minkowskiSum;
}

private bool closest_is_internal(Entry e )
{
    return (e.s > 0 && e.t > 0);
}

private Entry constructEntry(Vector A, Vector B, Vector p0, Vector q0, Vector p1, Vector q1)
{
    Entry e;

    e.y0 = A;
    e.y1 = B;

    e.p0 = p0;
    e.p1 = p1;

    e.q0 = q0;
    e.q1 = q1;

    Vector ab = B - A;
    double t = A.neg*ab;

    if ( t <= 0.0f)
    {
        t = 0.0f;
        e.v = A;
    }
    else
    {
        double denom = ab*ab;
        if (t >= denom)
        {
            e.v = B;
            t = 1.0f;
        }
        else
        {
            t /= denom;
            e.v = A + t * ab;
        }
    }

    e.s = 1-t;
    e.t = t;
    e.key = e.v.magnitude;

    return e;
}

// Expanding Polytope Algorithim (EPA)

Entry epa(Vector[] p1, Vector[] p2, inout ArraySeq!(Vector) sAB, inout ArraySeq!(Vector) sA,
           inout ArraySeq!(Vector) sB)
{
    // We want the final Simplex containing the origin
    while (sAB.size > 3)
    {
        sAB.removeHead;
        sA.removeHead;
        sB.removeHead;
    }

    auto Qheap = new ArraySeq!(Entry);

    for ( int i = 0; i < sAB.size; i++ )
    {
        int next_i = (i+1) % sAB.size;

        Entry e = constructEntry(sAB.get(i),sAB.get(next_i),sA.get(i),sB.get(i), sA.get(next_i),sB.get(next_i));

        if (closest_is_internal(e))
        {
            Qheap.append(e);
        }
    }

    bool close_enough = false;
    int failSafe = 0;
    Entry e;

    do
    {
        int j = 0;
        for (int i = 0; i < Qheap.size; i++)
        {
            if (Qheap.get(i).key < Qheap.get(j).key) j = i;
        }

        e = Qheap.get(j);
        Qheap.removeAt(j);

        Vector v = e.v;

        sAB.append(support(p1, p2, sA, sB, v));

        double vl = v.magnitude;

        double dot = v*sAB.tail/vl;
        close_enough = (dot - vl) <= EPSILON;

        if (close_enough == false)
        {
            Entry e1 = constructEntry(e.y0, sAB.tail, e.p0, e.q0, sA.tail, sB.tail);

            if (closest_is_internal(e1))
            {
                Qheap.append(e1);
            }

            Entry e2 = constructEntry(e.y1, sAB.tail, e.p1, e.q1, sA.tail, sB.tail);

            if (closest_is_internal(e2))
            {
                Qheap.append(e2);
            }
        }
    }
    while (close_enough == false && failSafe++ < 100);

    return e;
}
