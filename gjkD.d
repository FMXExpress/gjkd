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
module gjkD;

import tango.stdc.stringz;

import derelict.opengl.gl;
import derelict.opengl.glu;
import derelict.sdl.sdl;

import gjkSys;
import math;

const char[] WINDOW_TITLE = "gjkD v0.2";

//The screen attributes
const int SCREEN_WIDTH = 600;
const int SCREEN_HEIGHT = 600;
const int SCREEN_BPP = 32;

RigidSys system;

bool running;				// The main loop flag

//Module constructor.
static this()
{
	DerelictGL.load();		// Load Derelict libraries
	DerelictGLU.load();
    DerelictSDL.load();

    if (SDL_Init(SDL_INIT_VIDEO) < 0)
        throw new Exception("Failed to initialize SDL: " ~ getSDLError());
}

// Module destructor
static ~this()
{
	SDL_Quit();
}

void main(char[][] args)
{
	bool fullScreen = false;
  	system = new RigidSys(2);

    if (args.length > 1) fullScreen = args[1] == "-fullscreen";

	createGLWindow(WINDOW_TITLE, SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, fullScreen);
	initGL();

    running = true;
    while (running)			    // Main Program Loop
	{
		processEvents();        // User input
		system.update();   	    // Update Rigid Bodies
        drawScene();  		    // Draw Scene

        SDL_GL_SwapBuffers();
        SDL_Delay(10);          // Pause
    }
}

void processEvents()
{
    SDL_Event event;
    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
            case SDL_KEYUP:
                keyReleased(event.key.keysym.sym);
                break;
            case SDL_QUIT:
                running = false;
                break;
            default:
                break;
        }
    }
}

void keyReleased(int key)               // Controls
{
    switch (key)
    {
        case SDLK_ESCAPE:
            running = false;
            break;
        case SDLK_RIGHT:
            system.rb[1].vel.x += 5;
            break;
        case SDLK_LEFT:
            system.rb[1].vel.x -= 5;
            break;
        case SDLK_UP:
            system.rb[1].vel.y += 5;
            break;
        case SDLK_DOWN:
            system.rb[1].vel.y -= 5;
            break;
        case SDLK_SPACE:
            system.rb[1].vel.x = 0;
            system.rb[1].vel.y = 0;
            system.rb[1].omega.z = 0;
            system.rb[1].pos.x = 50;
            system.rb[1].pos.y = 50;
            system.rb[1].q.z = 0;
            break;
        case SDLK_RSHIFT:
            system.rb[1].omega.z += 0.01;
            break;
        case SDLK_RETURN:
            system.rb[1].omega.z -= 0.01;
            break;
        case SDLK_d:
            system.rb[0].vel.x += 5;
            break;
        case SDLK_a:
            system.rb[0].vel.x -= 5;
            break;
        case SDLK_w:
            system.rb[0].vel.y += 5;
            break;
        case SDLK_s:
            system.rb[0].vel.y -= 5;
            break;
        case SDLK_e:
            system.rb[0].omega.z += 0.01;
            break;
        case SDLK_q:
            system.rb[0].omega.z -= 0.01;
            break;
        case SDLK_c:
            system.rb[0].vel.x = 0;
            system.rb[0].vel.y = 0;
            system.rb[0].omega.z = 0;
            system.rb[0].pos.x = 50;
            system.rb[0].pos.y = 50;
            system.rb[0].q.z = 0;
            break;
        case SDLK_LEFTBRACKET:
            if(system.shape1 == 4) system.shape1 = 1;
            else system.shape1++;
            system.spawn(1);
            break;
        case SDLK_RIGHTBRACKET:
            if(system.shape2 == 4) system.shape2 = 1;
            else system.shape2++;
            system.spawn(2);
            break;
        default:
            break;
    }
}

void initGL()
{
   	glLoadIdentity();

	glMatrixMode( GL_PROJECTION );
	gluOrtho2D(0,100,0,100);				            // Use 2d Coordinate system
	glMatrixMode( GL_MODELVIEW );
    glDisable(GL_DEPTH_TEST);

    glShadeModel(GL_SMOOTH);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f); 			// Black Background

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glLoadIdentity();
}

void drawScene()
{
	// Clear The Screen And The Depth Buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glColor3f(0,0,1);						                        // Blue
	glLineWidth(2);

	glBegin(GL_LINES);						                        // Draw center Bullseye
		glVertex2d(50,55);
		glVertex2d(50,45);
		glVertex2d(45,50);
		glVertex2d(55,50);
	glEnd();

	glColor3f(1,0,0);						                        // Red

	foreach(int i, RigidBody b; system.rb)
	{
		if(i == 1) glColor3f(0,1,0);				                // Green

		glBegin(GL_LINE_LOOP);					                    // Draw Polygon
			for(int j; j < b.V[0].length; j++)
				glVertex2d(b.vertex[0][j].x, b.vertex[0][j].y);
		glEnd();

		glLoadIdentity();
		glFlush();
	}

	glTranslatef(50,50,0);

	if(system.collisionState == false) glColor3f(1.0, 1.0, 1.0);	// White - Clear
	else glColor3f(0, 0, 1);					                    // Blue - Collision

	glBegin(GL_LINE_LOOP);						                    // Draw Minkowski Hull
		int k = 0;
		foreach(Vector m; system.minkHull)
			if(m.x != 0 && m.y != 0) { glVertex2d(m.x, m.y); k++;}
	glEnd();

	glLoadIdentity();
	glFlush();
	glTranslatef(50,50,0);

	glColor3f(1,1,0);						                        // Yellow
	glBegin(GL_LINE_LOOP);						                    // Draw Simplex
		for(int j=system.simpIndex;(j>=0)&&(j>system.simpIndex-3);j--)
			glVertex2d(system.Bsimplex[j].x, system.Bsimplex[j].y);
	glEnd();

	glColor3f(1,0,0);						                        // Red
	glPointSize(6);

	glBegin(GL_POINTS);						                        // Draw closest point to origin
		glVertex2d(system.plot.x, system.plot.y);
	glEnd();

	glFlush();
}

void createGLWindow(char[] title, int width, int height, int bits, bool fullScreen)
{
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 6);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    SDL_WM_SetCaption(toStringz(title), null);

    int mode = SDL_OPENGL;
    if (fullScreen) mode |= SDL_FULLSCREEN;

    if (SDL_SetVideoMode(width, height, bits, mode) is null)
    {
        throw new Exception("Failed to open OpenGL window: " ~ getSDLError());
    }
}

char[] getSDLError()
{
    return fromUtf8z(SDL_GetError());
}
