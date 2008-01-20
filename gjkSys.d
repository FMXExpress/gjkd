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

class RigidSys
{
    RigidBody[] rb;
    Vector[] mink, minkHull;

    Vector range;
    Vector cp1,cp2;

    int shape1 = 1;			// Polygon #1 shape
    int shape2 = 1;			// Polygon #2 shape
    int stateSize = 12;
    double[][] minkSum;
    double[] x0, xEnd;
    bool penetrate;

    this(long MAXRB)
    {

        rb.length = MAXRB;
        rb[0] = new RigidBody(shape1);
        rb[1] = new RigidBody(shape2);

        minkSum = new double[][]((shape1+2)*(shape2+2),2);
        mink.length = minkHull.length = (shape1+2)*(shape2+2);

        rb[0].pos.x = 40f;
        rb[0].pos.y = 20f;
        rb[0].vel.x = 0f;
        rb[0].vel.y = 0f;

        rb[1].pos.x = 45f;
        rb[1].pos.y = 50f;
        rb[1].vel.x = -5f;
        rb[1].vel.y = -5f;

        x0.length = stateSize * rb.length;
        xEnd.length = stateSize * rb.length;
        bodiesToArray(xEnd);
    }

    void update()								                                // Update Universe
    {
        x0 = xEnd;

        rk4(0.01f);
        arrayToBodies(xEnd);

        foreach(inout RigidBody b; rb)
        {
            if (abs(b.omega.z) > .2f) b.omega.z *= .75f;                       // Limit rotation speed
            b.update();
        }

        // Narrow Phase Collision Detection

        auto sA = new ArraySeq!(Vector);    // Rigid Body 1 support map
        auto sB = new ArraySeq!(Vector);    // Rigid Body 2 support map
        auto sAB = new ArraySeq!(Vector);   // CSO
        Entry e;                            // Stores barycentric coordinates

        penetrate = gjk(rb[0].vertex[0], rb[1].vertex[0], sAB, sA, sB, e);

        if (penetrate == true)
            e = epa(rb[0].vertex[0], rb[1].vertex[0], sAB, sA, sB);

        cp1 = e.s * e.p0 + e.t * e.p1;
        cp2 = e.s * e.q0 + e.t * e.q1;
        range = e.v;

        /* cp1 and cp2 represent the closest points on each polytope.  If you normalize the range,
           you have the contact normal, which sets you up nicely for collision response */

        minkDiff();
    }

    void spawn(int hull)							                            // Change Polygon Shape
    {
        if (hull == 1)
        {

            rb[0].shape(shape1);
            minkSum = new double[][]((shape1+2)*(shape2+2),2);
            mink.length = minkHull.length = (shape1+2)*(shape2+2);
        }
        else
        {
            rb[1].shape(shape2);
            minkSum = new double[][]((shape1+2)*(shape2+2),2);
            mink.length = minkHull.length = (shape1+2)*(shape2+2);
        }
    }

    private void minkDiff()								                    // Calculate Minkowski Difference for display
    {
        int i = 0;
        for (int j; j < rb[0].vertex[0].length; j++)
            for (int k; k < rb[1].vertex[0].length; k++)
            {
                minkSum[i][0] = rb[0].vertex[0][j].x - rb[1].vertex[0][k].x;
                minkSum[i++][1] = rb[0].vertex[0][j].y - rb[1].vertex[0][k].y;
            }

        sort(minkSum);

        i = 0;
        foreach(inout Vector m; mink)
        {
            m.x = minkSum[i][0];
            m.y = minkSum[i++][1];
        }

        foreach(inout Vector v; minkHull)
        {
            v.x = 0;    // Clear Vector
            v.y = 0;
        }

        chainHull_2D(mink,minkHull);					            // Find Minkowski Hull
    }

    private void rk4(double h)				        // Runge Kutta 4th order ODE Solver
    {
        double[] inp, k1, k2, k3, k4;
        double hh, h6;

        inp.length = k1.length = k2.length = k3.length = k4.length = rb.length * stateSize;

        double x = 1f/6f;
        dxdt(h, x0, k1);			// Step 1
        foreach(int i, inout double p; inp)
        p = x0[i] + k1[i]*h*0.5f;
        dxdt(h, inp, k2);			// Step 2
        foreach(int i, inout double p; inp)
        p = x0[i] + k2[i]*h*0.5f;
        dxdt(h, inp, k3);			// Step 3
        foreach(int i, inout double p; inp)
        p = x0[i] + k3[i]*h;
        dxdt(h, inp, k4);			// Step 4
        foreach(int i, inout double p; xEnd)
        p = x0[i] + (k1[i] + 2f*k2[i] + 2f*k3[i] + k4[i])*h*x;

    }

    private void dxdt(double t, double[] x, inout double[] xdot)
    {
        arrayToBodies(x);

        foreach(int i, RigidBody b; rb)
        {
            ddtStateToArray(rb[i], xdot, i*stateSize);
        }
    }

    private void ddtStateToArray(RigidBody b, inout double[] xdot, int i)
    {
        xdot[i++] = b.vel.x;
        xdot[i++] = b.vel.y;
        xdot[i++] = b.vel.z;

        xdot[i++] = b.omega.x/2.0f;
        xdot[i++] = b.omega.y/2.0f;
        xdot[i++] = b.omega.z/2.0f;

        xdot[i++] = b.force.x;
        xdot[i++] = b.force.y;
        xdot[i++] = b.force.z;

        xdot[i++] = b.torque.x;
        xdot[i++] = b.torque.y;
        xdot[i++] = b.torque.z;
    }

    private void bodiesToArray(inout double[] x)
    {
        int i = 0;
        foreach(RigidBody b; rb)
        {
            x[i++] = b.pos.x;
            x[i++] = b.pos.y;
            x[i++] = b.pos.z;

            x[i++] = b.q.x;
            x[i++] = b.q.y;
            x[i++] = b.q.z;

            x[i++] = b.P.x;
            x[i++] = b.P.y;
            x[i++] = b.P.z;

            x[i++] = b.L.x;
            x[i++] = b.L.y;
            x[i++] = b.L.z;
        }
    }

    private void arrayToBodies(double[] x)
    {
        int i = 0;
        foreach(inout RigidBody b; rb)
        {
            b.pos.x = x[i++];
            b.pos.y = x[i++];
            b.pos.z = x[i++];

            b.q.x = x[i++];
            b.q.y = x[i++];
            b.q.z = x[i++];

            b.P.x = x[i++];
            b.P.y = x[i++];
            b.P.z = x[i++];

            b.L.x = x[i++];
            b.L.y = x[i++];
            b.L.z = x[i++];
        }
    }
}

private class RigidBody
{
    // Constant variables
    double mass;
    double iBody;
    double iBodyInv;

    Vector[][] V;
    Vector[][] vertex;

    // State variables
    Vector pos;					// Position of center of mass
    Vector q;					// rotation
    Vector P;					// linear momentum
    Vector L;					// angular momentum

    // Derived quantities (auxiliary variables)
    Vector iInv;				// Inverse of current I
    Vector vel;					// linear velocity
    Vector omega;				// angular velocity

    // Computed quantities
    Vector force;
    Vector torque;
    Vector collisionPoint;

    real degrees;

    this(int s)
    {
        shape(s);
    }

    void shape(int hull)
    {
        switch (hull)
        {
        case 1:		// Triangle
        {
            vertex = new Vector[][](1,3);
            V = new Vector[][](1,3);

            V[0][0].x = 0.0f;
            V[0][0].y = 1.0f;
            V[0][1].x = 1.0f;
            V[0][1].y = -1.0f;
            V[0][2].x = -1.0f;
            V[0][2].y = -1.0f;
            break;
        }
        case 2:		// Quad
        {
            vertex = new Vector[][](1,4);
            V = new Vector[][](1,4);

            V[0][0].x = 1.0f;
            V[0][0].y = 1.0f;
            V[0][1].x = 1.0f;
            V[0][1].y = -1.0f;
            V[0][2].x = -1.0f;
            V[0][2].y = -1.0f;
            V[0][3].x = -1.0f;
            V[0][3].y = 1.0f;
            break;

        }
        case 3:		// Pentagon
        {
            vertex = new Vector[][](1,5);
            V = new Vector[][](1,5);

            V[0][0].x = 1.0f;
            V[0][0].y = 1.0f;
            V[0][1].x = 2.0f;
            V[0][1].y = 0.0f;
            V[0][2].x = 0.0f;
            V[0][2].y = -2.0f;
            V[0][3].x = -2.0f;
            V[0][3].y = 0.0f;
            V[0][4].x = -1.0f;
            V[0][4].y = 1.0f;
            break;
        }
        case 4:		// Hexagon
        {
            vertex = new Vector[][](1,6);
            V = new Vector[][](1,6);

            V[0][0].x = 1.0f;
            V[0][0].y = 1.0f;
            V[0][1].x = 1.5f;
            V[0][1].y = 0.0f;
            V[0][2].x = 0.5f;
            V[0][2].y = -3.0f;
            V[0][3].x = -0.5f;
            V[0][3].y = -3.0f;
            V[0][4].x = -1.5f;
            V[0][4].y = 0.0f;
            V[0][5].x = -1.0f;
            V[0][5].y = 1.0f;
            break;
        }
        }
    }

    void update()						    // Update world coordinates
    {
        degrees = q.z * 180f/PI;			// convert Polar rotation to cartesian coordinates

        while (degrees > 360f) degrees -= 360f;
        while (degrees < -360f) degrees += 360f;

        real cd = cos(degrees);
        real sd = sin(degrees);

        for (int i = 0; i < V[0].length; i++)
        {
            vertex[0][i].x = pos.x + SCALE*(V[0][i].x*cd + V[0][i].y*sd);
            vertex[0][i].y = pos.y + SCALE*(-V[0][i].x*sd + V[0][i].y*cd);
        }
    }
}

