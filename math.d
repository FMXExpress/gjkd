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
module math;

import tango.math.Math;
import tango.util.collection.ArraySeq;

const EPSILON = float.epsilon;

struct Vector
{
    float x = 0f, y = 0f;

    static Vector opCall(float ax, float ay)
    {
        Vector u;
        u.x = ax;
        u.y = ay;
        return u;
    }

    void normalize()
    {
        float m = magnitude();
        x *= 1.0f/m;
        y *= 1.0f/m;
    }

    Vector opAdd(Vector u)
    {
        return Vector(x + u.x, y + u.y);
    }

    void opAddAssign(Vector u)
    {
        x += u.x;
        y += u.y;
    }

    void opSubAssign(Vector u)
    {
        x -= u.x;
        y -= u.y;
    }

    Vector opSub(Vector u)
    {
        return Vector(x - u.x, y - u.y);
    }

    float opMul(Vector u)			// Vector Dot Product
    {
        return(x*u.x + y*u.y);
    }

    Vector opMul(float s)			// Scaler Multiplication
    {
        return Vector(x*s, y*s);
    }

    float opXor(Vector v)
    {
        return x * v.y - y * v.x;
    }

    Vector opDiv(float s)
    {
        return Vector(x/s, y/s);
    }

    // negation
    Vector opNeg()
    {
        return Vector(-x, -y);
    }

    float magnitude()
    {
        return sqrt(x*x + y*y);
    }
}

struct Entry
{

    Vector y0, y1;
    Vector p0, p1;
    Vector q0, q1;
    Vector v;

    float key;
    float t;
    float s;
}
