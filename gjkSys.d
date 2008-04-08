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
module gjkSys;

import tango.math.Math;
import tango.core.Array;
import tango.util.collection.ArraySeq;

import math;
import collide;
import chainHull;

const SCALE = 5;            // Poltgon scale factor
const CIRCLE_SEGS = 25;

class RigidSys
{
    RigidBody[] rb;
    Vector[] mink, minkHull;

    Vector range;
    Vector cp1,cp2;

    int shape1 = 1;			// Polygon #1 shape
    int shape2 = 5;			// Polygon #2 shape
    float[][] minkSum;
    bool penetrate;

    this(long MAXRB)
    {
        rb.length = MAXRB;
        rb[0] = new RigidBody(shape1);
        rb[1] = new RigidBody(shape2);

        rb[0].pos = Vector(40.0f, 20.0f);
        rb[0].vel = Vector(0.0f,0.0f);
        rb[0].omega = 0.01;

        rb[1].pos = Vector(40.0f, 40.0f);
        rb[1].vel = Vector(0.1f,0.0f);
        rb[1].omega = 0.01;
    }

    void update()								                                // Update Universe
    {
        float dt = 1.0f/60.0f;

        rb[0].update(dt);
        rb[1].update(dt);

        // Narrow Phase Collision Detection

        Vector[] sA;    // Rigid Body 1 support map
        Vector[] sB;    // Rigid Body 2 support map
        Vector[] sAB;   // CSO;
        Entry e;        // Stores barycentric coordinates

        penetrate = gjk(rb[0], rb[1], sAB, sA, sB, e);

        if (penetrate)
            e = epa(rb[0], rb[1], sAB, sA, sB);

        cp1 = e.s * e.p0 + e.t * e.p1;
        cp2 = e.s * e.q0 + e.t * e.q1;
        range = e.v;

        /* cp1 and cp2 represent the closest points on each polytope.  If you normalize the range,
           you have the contact normal, which sets you up nicely for collision response */

        minkDiff();
    }

    void spawn(int hull)							                            // Change Polygon Shape
    {
        if (hull == 1) rb[0].shape(shape1);
        else rb[1].shape(shape2);
    }

    private void minkDiff()								                    // Calculate Minkowski Difference for display
    {
        int scale = rb[1].vertex.length*rb[0].vertex.length;
        minkSum = new float[][](scale,2);
        mink.length = minkHull.length = scale;

        int i = 0;
        foreach(rb1; rb[1].vertex)
            foreach(rb2; rb[0].vertex)
            {
                minkSum[i][0] = rb2.x - rb1.x;
                minkSum[i++][1] = rb2.y - rb1.y;
            }

        sort(minkSum);

        i = 0;
        foreach(inout m; mink)
        {
            m.x = minkSum[i][0];
            m.y = minkSum[i++][1];
        }

        foreach(inout v; minkHull)
        {
            v.x = 0;    // Clear Vector
            v.y = 0;
        }

        chainHull_2D(mink,minkHull);					            // Find Minkowski Hull
    }

}

private class RigidBody
{

    Vector[] V;
    Vector[] vertex;

    // State variables
    Vector pos;					// Position of center of mass
    float q;					// Rotation position

    // Derived quantities (auxiliary variables)
    Vector vel;					// linear velocity
    float omega;				// angular velocity

    int type;

    float radius;

    this(int s)
    {
        type = s;
        shape(type);
        q = 0.0001f;
        transform();
    }

    void shape(int hull)
    {
        type = hull;
        switch (hull)
        {
        case 1:		// Triangle
        {
            V = null;
            vertex = null;
            V ~= Vector(0,1);
            V ~= Vector(1,-1);
            V ~= Vector(-1,-1);
            vertex.length = V.length;
            break;
        }
        case 2:		// Quad
        {
            V = null;
            vertex = null;
            V ~= Vector(1,1);
            V ~= Vector(1,-1);
            V ~= Vector(-1,-1);
            V ~= Vector(-1,1);
            vertex.length = V.length;
            break;

        }
        case 3:		// Pentagon
        {
            V = null;
            vertex = null;
            V ~= Vector(1,1);
            V ~= Vector(2,0);
            V ~= Vector(0,-2);
            V ~= Vector(-2,0);
            V ~= Vector(-1,1);
            vertex.length = V.length;
            break;
        }
        case 4:		// Hexagon
        {
            V = null;
            vertex = null;
            V ~= Vector(1,1);
            V ~= Vector(1.5,0);
            V ~= Vector(0.5,-3);
            V ~= Vector(-0.5,-3);
            V ~= Vector(-1.5,0);
            V ~= Vector(-1,1);
            vertex.length = V.length;
            break;
        }
        case 5:		// Circle
        {
            radius = 1.5;
            V = null;
            vertex = null;

            int segs = CIRCLE_SEGS;
            Vector c = pos;
            float r = radius;
            float coef = 2.0*PI/segs;

            for(int n = 0; n <= segs; n++)
            {
                float rads = n*coef;
                V ~= Vector(r*cos(rads), r*sin(rads));
            }
            vertex.length = V.length;

            break;
        }
        }
    }

    void update(float dt)
    {
        pos.x += vel.x*dt;
        pos.y += vel.y*dt;
        q += omega*dt;

        transform();
    }

    void transform()
    {
        // Update world coordinates
        float degrees = q * 180f/PI;			// convert Polar rotation to cartesian coordinates

        while (degrees > 360f) degrees -= 360f;
        while (degrees < -360f) degrees += 360f;

        float cd = cos(degrees);
        float sd = sin(degrees);

        foreach(int i, v; V)
        {
            vertex[i].x = pos.x + SCALE*(v.x*cd + v.y*sd);
            vertex[i].y = pos.y + SCALE*(-v.x*sd + v.y*cd);
        }
    }

    Vector support(Vector n)
    {
        Vector r;
        if(type == 5)
        {
            r = radius * n.getNormal * SCALE;
            r = r + pos;
        }
        else
        {
            int i = vertex.length-1;
            r = vertex[i--];
            while (i>=0)
            {
                if ( (vertex[i] - r) * n >= 0 )
                {
                    r = vertex[i];
                }
                i--;
            }
        }
        return r;
    }

    Vector coldStartGjk(Vector p)
    {
        if(type == 5)
        {
            return support(p);
        }
        else
        {
            return vertex[0];
        }
    }

    Vector getCenter()
    {
        return pos;
    }
}
