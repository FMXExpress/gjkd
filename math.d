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

const EPSILON = 0.00008f;

struct Vector
{
    double x = 0f, y = 0f, z = 0f;

    static Vector opCall(double ax, double ay, double az)
    {
        Vector u;
        u.x = ax;
        u.y = ay;
        u.z = az;
        return u;
    }

    void normalize()
    {
        double m = magnitude();
        x *= 1.0f/m;
        y *= 1.0f/m;
        z *= 1.0f/m;
    }

    Vector opAdd(Vector u)
    {
        return Vector(x + u.x, y + u.y, z + u.z);
    }

    void opAddAssign(Vector u)
    {
        x += u.x;
        y += u.y;
        z += u.z;
    }

    void opSubAssign(Vector u)
    {
        x -= u.x;
        y -= u.y;
        z -= u.z;
    }

    Vector opSub(Vector u)
    {
        return Vector(x - u.x, y - u.y, z - u.z);
    }

    real opMul(Vector u)			// Vector Dot Product
    {
        return(x*u.x + y*u.y + z*u.z);
    }

    Vector opMul(double s)			// Scaler Multiplication
    {
        return Vector(x*s, y*s, z*s);
    }

    Vector opXor(Vector u)
    {
        return Vector(y*u.z - z*u.y, -x*u.z + z*u.x, x*u.y - y*u.x);
    }

    Vector opDiv(double s)
    {
        return Vector(x/s, y/s, z/s);
    }

    Vector neg()
    {
        return Vector(-x, -y, -z);
    }

    real magnitude()
    {
        return sqrt(x*x + y*y + z*z);
    }
}

struct Entry
{

    Vector y0, y1;
    Vector p0, p1;
    Vector q0, q1;
    Vector v;

    double key;
    double t;
    double s;
}
