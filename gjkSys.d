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

module gjkSys;

import tango.io.Stdout;
import tango.math.Random;
import tango.math.Math;
import tango.util.collection.HashSet;
import tango.core.Array;

import math;
import collide;
import chainHull;

const bool hash = false;		// Utilize Hash table
const int GRID_SIZE = 16;

const float X = 2.5f;			// X polygon scale factor
const float Y = 2.5f;			// Y polygon scale factor


class RigidSys
{
	RigidBody[] rb; 
	Vector[] mink, minkHull;
	Vector[] Bsimplex;
	Vector plot;

	int simpIndex = 0;
	int shape1 = 1;			// Polygon #1 shape
	int shape2 = 1;			// Polygon #2 shape
	int stateSize = 12;
	float[][] minkSum;
	float[] x0, xEnd;
	bool collisionState;
	
	this(long MAXRB)
	{
				
		rb.length = MAXRB;
		Bsimplex.length = 10;
				
		rb[0] = new RigidBody(shape1);
		rb[1] = new RigidBody(shape2);

		minkSum = new float[][]((shape1+2)*(shape2+2),2);
		mink.length = minkHull.length = (shape1+2)*(shape2+2);

		rb[0].pos.x = 40f;
		rb[0].pos.y = 20f;
		rb[0].vel.x = 0f;
		rb[0].vel.y = 0f;

		rb[1].pos.x = 45f;
		rb[1].pos.y = 30f;
		rb[1].vel.x = 5f;
		rb[1].vel.y = 5f;

	       	x0.length = stateSize * rb.length;
		xEnd.length = stateSize * rb.length;		
		bodiesToArray(xEnd);
				
	}
	
	void update()								// Update Universe
	{
		x0 = xEnd;
		
		rk4(0.01f);
		arrayToBodies(xEnd);

		foreach(inout RigidBody b; rb)
		{
			if(abs(b.omega.z) > .2f) b.omega.z *= .75f;
			b.update();
		}

		foreach(inout Vector v; Bsimplex) { v.x = 0; v.y = 0; }		// Clear Simplex
		collisionState = checkForCollisions();
		minkDiff();

	}
	
	void spawn(int hull)							// Change Polygon Shape
	{
		if(hull == 1) 
		{
		
			rb[0].shape(shape1);
			minkSum = new float[][]((shape1+2)*(shape2+2),2);
			mink.length = minkHull.length = (shape1+2)*(shape2+2);
		}
		else
		{
			rb[1].shape(shape2);
			minkSum = new float[][]((shape1+2)*(shape2+2),2);
			mink.length = minkHull.length = (shape1+2)*(shape2+2);
		}

	}

	void minkDiff()								// Calculate Minkowski Difference
	{	
		int i = 0;
		for(int j; j < rb[0].vertex[0].length; j++)
			for(int k; k < rb[1].vertex[0].length; k++)
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

		foreach(inout Vector v; minkHull) { v.x = 0; v.y = 0; }		// Clear Vector

		chainHull_2D(mink,minkHull);					// Find Minkowski Hull
	}
	
	void rk4(float h)				// Runge Kutta 4th order ODE Solver
	{
		float[] inp, k1, k2, k3, k4;
		float hh, h6;

		inp.length = k1.length = k2.length = k3.length = k4.length = rb.length * stateSize;

		float x = 1f/6f;
		dxdt(h, x0, k1);			// Step 1
		foreach(int i, inout float p; inp)
			p = x0[i] + k1[i]*h*0.5f;
		dxdt(h, inp, k2);			// Step 2
		foreach(int i, inout float p; inp)
			p = x0[i] + k2[i]*h*0.5f;
		dxdt(h, inp, k3);			// Step 3
		foreach(int i, inout float p; inp)
			p = x0[i] + k3[i]*h;
		dxdt(h, inp, k4);			// Step 4
		foreach(int i, inout float p; xEnd)
			p = x0[i] + (k1[i] + 2f*k2[i] + 2f*k3[i] + k4[i])*h*x;

	}

	void dxdt(float t, float[] x, inout float[] xdot)
	{
		arrayToBodies(x);
		
		foreach(int i, RigidBody b; rb)
		{
			ddtStateToArray(rb[i], xdot, i*stateSize);
		}
	}

	void ddtStateToArray(RigidBody b, inout float[] xdot, int i)
	{
		xdot[i++] = b.vel.x;
		xdot[i++] = b.vel.y;
		xdot[i++] = b.vel.z;

		xdot[i++] = b.omega.x;
		xdot[i++] = b.omega.y;
		xdot[i++] = b.omega.z;

		xdot[i++] = b.force.x;
		xdot[i++] = b.force.y;
		xdot[i++] = b.force.z;

		xdot[i++] = b.torque.x;
		xdot[i++] = b.torque.y;
		xdot[i++] = b.torque.z;
	}
	
	void bodiesToArray(inout float[] x)
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

	void arrayToBodies(float[] x)
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

	void computeForceAndTorque(float t, inout RigidBody b)
	{
		auto rand = new Random;
	 
		b.force.x = rand.next(10);
		b.force.y = rand.next(10);
		b.force.z = rand.next(10);
	
		b.torque.x = rand.next(100);
		b.torque.y = rand.next(100);
		b.torque.z = rand.next(100); 
	}

	bool checkForCollisions()
	{
			
		bool collide = false;
	
		if(hash == true)							// Hash Table for Collisions
		{
			alias HashSet!(int) set;
			set[GRID_SIZE] hashTable;
	
			for(int i = 0; i < GRID_SIZE; i++)
				hashTable[i] = new set;

			foreach(int i, RigidBody b; rb)
			{
			
				if(b.pos.x > 100) b.vel.x -= 5;
					else if(b.pos.x < 0) b.vel.x += 5;

				if(b.pos.y > 100) b.vel.y -= 5;
					else if(b.pos.y < 0) b.vel.y += 5;
		
				float aabb4 = floor((b.pos.x + X) * 0.04) + floor((b.pos.y + Y) * 0.04) * 4;
				float aabb1 = floor((b.pos.x - X) * 0.04) + floor((b.pos.y - Y) * 0.04) * 4;
			
				if((aabb4 < GRID_SIZE) && (aabb1 >= 0))
				{

					if(aabb4 == aabb1)				// Single Cell
					{
						hashTable[rndint(aabb4)].add(i);
					}
					else if(aabb4 == aabb1+1)			// Two Cells, adjacent
					{
						hashTable[rndint(aabb1)].add(i);
						hashTable[rndint(aabb4)].add(i);
					}
					else if(aabb4 == aabb1+4)			// Two Cells, over & under 
					{
						hashTable[rndint(aabb1)].add(i);
						hashTable[rndint(aabb4)].add(i);
					}
					else						// Four cells
					{
						float aabb3 = floor((b.pos.x - X) * 0.04) + floor((b.pos.y + Y) * 0.04) * 4;
						float aabb2 = floor((b.pos.x + X) * 0.04) + floor((b.pos.y - Y) * 0.04) * 4;

						hashTable[rndint(aabb1)].add(i);
						hashTable[rndint(aabb2)].add(i);
						hashTable[rndint(aabb3)].add(i);
						hashTable[rndint(aabb4)].add(i);
					}
				}

			}


			for(int i = 0; i < hashTable.length; i++)
			{
				if(hashTable[i].size() > 1)
				{

					int[] t;
					int x = 0;
				
					t.length = hashTable[i].size();
				
					foreach(value; hashTable[i]) t[x++] = value;
									
					for(int j = 0; j < t.length-1; j++)
						for(int k=j+1; k < t.length; k++) collide = testAABB(t[j],t[k]);			
				}
			}
		}
		else
		{
			foreach(int i, RigidBody b; rb)
			{
				if(b.vel.x > 50) b.vel.x = 50f;
				if(b.vel.y > 50) b.vel.y = 50f;

				if(b.pos.x > 100) b.vel.x -= 5f;
				else if(b.pos.x < 0) b.vel.x += 5f;

				if(b.pos.y > 100) b.vel.y -= 5f;
				else if(b.pos.y < 0) b.vel.y += 5f;
			}

			for(int j = 0; j < rb.length-1; j++)
				for(int k=j+1; k < rb.length; k++) collide = testAABB(j,k);
		}

		return collide;
	}

	bool testAABB(int b1, int b2)
	{
		return gjk(rb[b1].vertex[0], rb[b2].vertex[0], Bsimplex, plot, simpIndex);
	}

}

class RigidBody
{
	// Constant variables
	float mass;
	float iBody;
	float iBodyInv;
	
	Vector[][] V;
	Vector[][] vertex;

	// State variables
	Vector pos;					// Position of center of mass
	Vector q;					// rotation
	Vector P;					// linear momentum
	Vector L;					// angular momentum

	// Derived quantities (auxiliary variables)
	Vector iInv;					// Inverse of current I
	Vector vel;					// linear velocity
	Vector omega;					// angular velocity

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
		switch(hull)
		{
			case 1:		// Triangle
			{
				vertex = new Vector[][](1,3);
				V = new Vector[][](1,3);

				V[0][0].x = 0f;
				V[0][0].y = Y; 
				V[0][1].x = X;
				V[0][1].y = -Y;
				V[0][2].x = -X;
				V[0][2].y = -Y;
				break;
			}
			case 2:		// Quad
			{
				vertex = new Vector[][](1,4);
				V = new Vector[][](1,4);

				V[0][0].x = X;
				V[0][0].y = Y; 
				V[0][1].x = X;
				V[0][1].y = -Y;
				V[0][2].x = -X;
				V[0][2].y = -Y;
				V[0][3].x = -X;
				V[0][3].y = Y;
				break;

			}
			case 3:		// Pentagon
			{
				vertex = new Vector[][](1,5);
				V = new Vector[][](1,5);

				V[0][0].x = X;
				V[0][0].y = Y; 
				V[0][1].x = 2f*X;
				V[0][1].y = 0f;
				V[0][2].x = 0f;
				V[0][2].y = -2f*Y;
				V[0][3].x = -2f*X;
				V[0][3].y = 0f;
				V[0][4].x = -X;
				V[0][4].y = Y;
				break;
			}
			case 4:		// Hexagon
			{
				vertex = new Vector[][](1,6);
				V = new Vector[][](1,6);

				V[0][0].x = X;
				V[0][0].y = Y; 
				V[0][1].x = 1.5f*X;
				V[0][1].y = 0f;
				V[0][2].x = 0.5f*X;
				V[0][2].y = -3f*Y;
				V[0][3].x = -0.5f*X;
				V[0][3].y = -3f*Y;
				V[0][4].x = -1.5f*X;
				V[0][4].y = 0;
				V[0][5].x = -X;
				V[0][5].y = Y;
				break;
			}
		}
	}

	void update()						// Update coordinates
	{

		degrees = q.z * 180f/PI;			// convert Polar rotation to cartesian coordinates
		

		while(degrees > 360f) degrees -= 360f;
		while(degrees < -360f) degrees += 360f;
 
		real cd = cos(degrees);
		real sd = sin(degrees);
	
		for(int i = 0; i < V[0].length; i++)
		{
			vertex[0][i].x = pos.x + (V[0][i].x*cd + V[0][i].y*sd);
			vertex[0][i].y = pos.y + (-V[0][i].x*sd + V[0][i].y*cd);
		}
		
	}

}

