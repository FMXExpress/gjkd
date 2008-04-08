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
module collide;

import tango.util.collection.LinkSeq;
import tango.math.Math;

import math;
import gjkSys;

const SIMPLEX_EPSILON = 0.01f;

// The Gilbert-Johnson-Keerthi algorithm

bool gjk(RigidBody rBody1, RigidBody rBody2, inout Vector[] sAB, inout Vector[] sA, inout Vector[] sB, inout Entry e)
{
    sA  ~= rBody1.coldStartGjk(rBody2.pos);
    sB  ~= rBody2.coldStartGjk(rBody1.pos);
    sAB ~= (sA[0] - sB[0]);

    e = constructEntry(sAB[0], sAB[0], sA[0], sB[0], sA[0], sB[0]);

    bool penetrate;
    int  failsafe;

    while (failsafe++ < 20)       /// Don't want to get caught in an infinite loop!
    {
        Vector v1 = rBody1.support(-e.v);
        Vector v2 = rBody2.support(e.v);

        sA ~= v1; sB ~= v2;
        sAB ~= v1 - v2;

        if (e.v*e.v - e.v*sAB[sAB.length - 1] < SIMPLEX_EPSILON) return false;

        /// Line Test
        if (sAB.length == 2)
        {
            e = constructEntry(sAB[0], sAB[1], sA[0], sB[0], sA[1], sB[1]);
        }
        /// Triangle Test
        else e = pointTriangle(sAB, sA, sB, penetrate);

        if (penetrate) return true;
    }
    // This should never happen... ;-)
    return false;
}

///
private Entry pointTriangle(Vector[] sAB, Vector[] sA, Vector[] sB, inout bool penetrate)
{
    int      i  = sAB.length - 1;
    Vector ab = sAB[i - 1] - sAB[sAB.length - 1];
    Vector ac = sAB[i - 2] - sAB[sAB.length - 1];
    Vector ao = -sAB[sAB.length - 1];

    float   d1 = ab*ao;
    float   d2 = ac*ao;

    /// Origin in vertex region outside A
    if (d1 <= 0.0f && d2 <= 0.0f)
        return constructEntry(sAB[i], sAB[sAB.length - 1], sA[i], sB[i], sA[sA.length - 1], sB[sB.length - 1]);

    Vector bo = -sAB[i - 1];
    float   d3 = ab*bo;
    float   d4 = ac*bo;
    float   vc = d1 * d4 - d3 * d2;
    /// Origin in edge region outside AB
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
        return constructEntry(sAB[i - 1], sAB[sAB.length - 1], sA[i - 1], sB[i - 1], sA[sA.length - 1], sB[sB.length - 1]);

    Vector co = -sAB[i - 2];
    float   d5 = ab*co;
    float   d6 = ac*co;
    float   vb = d5 * d2 - d1 * d6;
    /// Origin in edge region ouside AC
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
        return constructEntry(sAB[i - 2], sAB[sAB.length - 1], sA[i - 2], sB[i - 2], sA[sA.length - 1], sB[sB.length - 1]);
    /// Origin inside face region. Penetration!!!
    penetrate = true;
    Entry fooBar;
    return fooBar;
}

///
private bool closest_is_internal(Entry e)
{
    return e.s > 0 && e.t > 0;
}

///
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
    float   t  = -A*ab;

    if (t <= 0.0f)
    {
        t   = 0.0f;
        e.v = A;
    }
    else
    {
        float denom = ab*ab;
        if (t >= denom)
        {
            e.v = B;
            t   = 1.0f;
        }
        else
        {
            t  /= denom;
            e.v = A + t * ab;
        }
    }

    e.s     = 1 - t;
    e.t     = t;
    e.key = e.v.magnitude;
    return e;
}

/// Expanding Polytope Algorithim (EPA)
Entry epa(RigidBody rBody1, RigidBody rBody2, inout Vector[] sAB, inout Vector[] sA, inout Vector[] sB)
{
    /// Line
    if (sAB.length == 2)
    {
        /// This only works when you have a shallow penetration. Need to optimize this for deep penetrations as well
        if (sAB[0].magnitude < sAB[sAB.length - 1].magnitude)
            return constructEntry(sAB[0], sAB[0], sA[0], sB[0], sA[0], sB[0]);
        else
            return constructEntry(sAB[sAB.length - 1], sAB[sAB.length - 1], sA[sA.length - 1], sB[sB.length - 1], sA[sA.length - 1], sB[sB.length - 1]);
    }

    /// We want the final Simplex containing the origin
    if (sAB.length > 3)
    {
        Vector[] A = sA[sA.length - 3 .. sA.length];
        Vector[] B = sB[sB.length - 3 .. sB.length];
        Vector[] C = sAB[sAB.length - 3 .. sAB.length];
        sA = A; sB = B; sAB = C;
    }

    auto Qheap = new LinkSeq!(Entry);

    for (int i = 0; i < sAB.length; i++)
    {
        int   next_i = (i + 1) % sAB.length;
        Entry e      = constructEntry(sAB[i], sAB[next_i], sA[i], sB[i], sA[next_i], sB[next_i]);

        if (closest_is_internal(e))
            Qheap.append(e);
    }

    bool  close_enough = false;
    int   failSafe     = 0;
    Entry e;

    do
    {
        int j = 0;
        for (int i = 0; i < Qheap.size; i++)
            if (Qheap.get(i).key < Qheap.get(j).key) j = i;

        e = Qheap.get(j);
        Qheap.removeAt(j);

        Vector v = e.v;

        Vector v1 = rBody1.support(e.v);
        Vector v2 = rBody2.support(-e.v);

        sA ~= v1; sB ~= v2;
        sAB ~= v1 - v2;

        float vl = v.magnitude;
        float dot = v*sAB[sAB.length - 1] / vl;

        close_enough = (dot - vl) <= SIMPLEX_EPSILON;

        if (!close_enough)
        {
            Entry e1 = constructEntry(e.y0, sAB[sAB.length - 1], e.p0, e.q0, sA[sA.length - 1], sB[sB.length - 1]);

            if (closest_is_internal(e1))
                Qheap.append(e1);


            Entry e2 = constructEntry(e.y1, sAB[sAB.length - 1], e.p1, e.q1, sA[sA.length - 1], sB[sB.length - 1]);
            if (closest_is_internal(e2))
                Qheap.append(e2);
        }
    } while (!close_enough && failSafe++ < 10);
    return e;
}

