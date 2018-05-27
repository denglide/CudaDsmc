/*
The MIT License (MIT)

Copyright (c) 2013 Denis Gladkov

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <cuda_gl_interop.h>

#include "model.h"
#include "sim_settings.h"
#include "draw_helper.h"

#include "dsmc_base.h"
#include "dsmc_kernel.cuh"

uint	g_width;
uint	g_height;

bool	g_bUseBWScheme = false;

//TODO: move to frame object
float3	position = make_float3(0,0,-37);
float3	rotation = make_float3(45.0f,45.0f,0);

mouse_state_t		g_mouse_state;
settings_t			g_settings;
simulation_state_t	g_simState;

float3	g_streamVel = make_float3(0.0,0,0);

//TODO: move to benchmarking object
cudaEvent_t start, stop;
float	gpuTime = 0;

Model	model;
dsmc::DSMCSolver*	solver = 0;

void	reshape(int w, int h)
{
	g_width = w;
	g_height = h;

	ReshapePerspective(w, h);
}

void	display()
{
	if(g_bUseBWScheme)
		glClearColor(1,1,1,1);
	else
		glClearColor(0,0,0,0);

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();

	glColor3f(1,0,0);

	glTranslatef(position.x,position.y, position.z);

	glRotatef(rotation.x,1,0,0);
	glRotatef(rotation.y,0,1,0);

	glPointSize(1);

	solver->Render();

	if(g_simState.drawVectorField)
	{
		glColor3f(0,0,1);
		solver->RenderVectorField();
	}

	if(g_bUseBWScheme)
	{
		glColor3f(0,0,0);
		DrawBBox(g_settings.boundaries);
	}
	else
		DrawBBoxColored(g_settings.boundaries);

	if(g_settings.geometry)
		RenderModel(model);

	glColor3f(0,0,1);

	if(g_simState.drawGrid)
		DrawRegularGrid(model.GetGrid());

	if(g_simState.sampleConcentration || g_simState.drawVectorField)
		DrawSimulationGrid(g_settings.boundaries,g_settings.grid_dim, g_simState.slice);

	glLoadIdentity();

	if(g_bUseBWScheme)
		glColor3f(0,0,0);
	else
		glColor3f(1,0,0);

	glDisable(GL_DEPTH_TEST);

	setup2DView(800,600);

	DrawText2D(10, 10, "Press left mouse button and drag to move");

	DrawText2D(300, 10, "Direct Simulation Monte Carlo");

	if(g_simState.paused)
		DrawText2D(300, 20, "Simulation is paused");

	DrawText2D(10, 20, "Press right mouse button and drag to rotate");
	DrawText2D(10, 30, "Press +/- to zoom");
	DrawText2D(10, 40, "Press space to toggle pause");

	static	char s[100];
	sprintf(s,"FPS: %4.2f",g_simState.fps);

	DrawText2D(10, 50, s);

	sprintf(s,"Simulation time: %f", g_simState.simTime);

	DrawText2D(10, 60, s);

	sprintf(s,"Real time: %4.2f", g_simState.realTime);

	DrawText2D(10, 70, s);

	if(g_simState.sampleConcentration)
	{
		DrawText2D(700, 10, "Concentration");
		solver->RenderConcMap(650, 20, 128, 128);
	}

	reshape(g_width, g_height);

	glEnable(GL_DEPTH_TEST);

	glFlush();

	glutSwapBuffers();
}

void	mouse(int btn, int state, int x, int y)
{
	if(state == GLUT_DOWN)
	{
		if(btn == GLUT_LEFT_BUTTON)
			g_mouse_state.leftBtnPressed = true;
		if(btn == GLUT_RIGHT_BUTTON)
			g_mouse_state.rightBtnPressed = true;

		g_mouse_state.pos_x = x;
		g_mouse_state.pos_y = y;
	}else
	{
		if(btn == GLUT_LEFT_BUTTON)
			g_mouse_state.leftBtnPressed = false;
		if(btn == GLUT_RIGHT_BUTTON)
			g_mouse_state.rightBtnPressed = false;

		g_mouse_state.pos_x = -1;
		g_mouse_state.pos_y = -1;
	}

	if(btn == 3)
		position.z += 0.5f;

	if(btn == 4)
		position.z -= 0.5f;
}

void mouse_motion(int x,int y)
{
	if(!g_simState.demoMode)
	{
		if(g_mouse_state.leftBtnPressed)
		{
			position.x += (x-g_mouse_state.pos_x)/100.0;
			position.y += (g_mouse_state.pos_y-y)/100.0;
		}

		if(g_mouse_state.rightBtnPressed)
		{
			rotation.x += (g_mouse_state.pos_y-y)/10.0;
			rotation.y += (x-g_mouse_state.pos_x)/10.0;
		}

		g_mouse_state.pos_x = x;
		g_mouse_state.pos_y = y;
	}
}

void	update()
{
	checkCUDAError("update");

	static uint prevTime = glutGet(GLUT_ELAPSED_TIME);
	static float alpha = 0.0f;
	static bool	stopZoom = false;
	static uint frame = 0;
	static	float tt = 0.0f;

	uint	time2 = glutGet(GLUT_ELAPSED_TIME);

	float	delta = (time2 - prevTime)/1000.0;
	frame++;

	tt += delta;

	if(g_settings.maxRunTime > 0.0f && g_simState.simTime > g_settings.maxRunTime)
		throw	std::runtime_error("Running time is out");

	if (tt >= 1.0f)
	{
			g_simState.fps = frame/tt;
			tt = 0.0f;
			frame = 0;
	}

	if(g_simState.demoMode)
	{
		rotation.y -= delta*10;//0.5;
		if(rotation.x > 0)
			rotation.x -= delta*3;//0.1;

		if(position.z < -30 && !stopZoom)
		{
			position.z += delta*3;//0.1f;
		}else
		{
			stopZoom = true;
			position.z = -30 + 15*sin(alpha);
			alpha+=delta*0.2f;//0.01f;
		}
	}

	prevTime = time2;

	if(!g_simState.paused)
	{
		g_simState.advanceRealTime(delta);
		solver->RunSimulationStep(g_settings.dt);
	}

	glutPostRedisplay();
}

void keyboard( unsigned char key, int x, int y) 
{
	if (key == 27)
		exit(1);
    if(key == ' ')
    	g_simState.togglePause();
    if(key == '+')
    	position.z++;
    if(key == '-')
    	position.z--;

    if(key == 'g')
    	g_simState.toggleDrawGrid();
    if(key == 'v')
    	g_simState.toggleDrawVectorField();

    if(key == 'c')
    	g_simState.toggleConcentration();

    if(key == 'd')
    	g_simState.toggleDemoMode();

    if(key == 'n')
    	g_simState.toggleCollisions();

    if(key == ']')
    	g_simState.incSlice(g_settings.grid_dim.z);
    if(key == '[')
    	g_simState.decSlice();
}

int main(int argc, char **argv)
{
	try
	{

	std::cout<<"Particles reordering is: ";
#ifdef	REORDER
	std::cout<<"ON\n";
#else
	std::cout<<"OFF\n";
#endif

	int devid = 0;

	g_settings.readFromCmdLine(argc, argv);

	float	cell_volume = g_settings.cell_size.x*g_settings.cell_size.y*g_settings.cell_size.z;
	uint	expNumOfPart = g_settings.density*cell_volume/g_settings.fnum;

	std::cout<<"Expected number of particles per cell: "<<expNumOfPart<<" density is "<<g_settings.density<<"\n";

	g_settings.partPerCell = expNumOfPart;
	g_settings.numParticles = expNumOfPart*g_settings.getCellsCount();

	g_settings.print();

	cudaError_t err = cudaSuccess;

	if(!g_settings.benchmark)
		err = cudaGLSetGLDevice(devid);
	else
		err = cudaSetDevice( devid );

	if (err != cudaSuccess)
	{
		std::cout<<"Device initialization falied\n";
		return 0;
	}

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, devid);

	printf("Using device: %s\n", deviceProp.name);

	g_width = 1024;
	g_height = 768;

	if(!g_settings.benchmark)
	{
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowPosition(100,100);
		glutInitWindowSize(g_width,g_height);
		glutCreateWindow("DSMC");
	
		glutDisplayFunc(display);
		glutIdleFunc(update);
		glutReshapeFunc(reshape);
		glutKeyboardFunc(keyboard);
		glutMouseFunc(mouse);
		glutMotionFunc(mouse_motion);

		glewInit();

		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_NORMALIZE);

		glEnable(GL_TEXTURE_2D);

		glShadeModel (GL_SMOOTH);

		glFrontFace(GL_CCW);

		glCullFace(GL_BACK);

		solver = dsmc::CreateSolver(dsmc::SolverGPUVBO);
	}else
		solver = dsmc::CreateSolver(dsmc::SolverGPU);


	std::cout<<"Initializing solver...";

	solver->Init(g_settings);

	std::cout<<" done\n";

	if(g_settings.geometry)
	{
		printf("Building the model...\n");

		if(!model.LoadRaw("data/ares.raw", !g_settings.benchmark))
			throw	std::runtime_error("Unable to load model");
		else
			printf("Model loaded\n");

		model.CreateGrid(vec3i_t(5,5,20));

		model.MoveModel(vec3_t(0,0,-3));

		solver->SetModel(&model);
	}

	solver->RegisterSimState(&g_simState);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	std::cout<<"Initializing the system state...\n";

	cudaEventRecord( start, 0 );

	g_settings.numParticles = solver->InitSystemState(make_float3(g_settings.streamVelocity.x,
																	g_settings.streamVelocity.y,
																		g_settings.streamVelocity.z));

	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	
	cudaEventElapsedTime( &gpuTime, start, stop );

	std::cout<<g_settings.numParticles<<" samples were generated in "<<gpuTime<<" ms\n";

	if(!g_settings.benchmark)
	{
		std::cout<<"Starting simulation loop\n";
		glutMainLoop();
	}
	else
	{
		uint	nIter = g_settings.maxRunTime / g_settings.dt + 1;
		std::cout<<"Starting benchmarking. Will run "<<nIter<<" iterations\n";

		cudaEventRecord( start, 0 );

		solver->DoBenchmark(nIter, g_settings.dt);

		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );

		cudaEventElapsedTime( &gpuTime, start, stop );

		std::cout<<"Benchmark is done in "<<gpuTime<<" ms\n";
	}
		
	std::ofstream	fout("statistics.txt");
	solver->DumpStatistics(fout);
	}
	catch(std::exception& e)
	{
		std::cerr<<"Exception: "<<e.what()<<"\n";
	}


	if(!g_settings.benchmark)
	{
		solver->SaveColorMap("colormap.png");
		solver->SaveVelocityField("velfield.png");
	}

	delete solver;

	return 0;
}
