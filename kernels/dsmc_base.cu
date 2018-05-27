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

#include "gpu_dsmc.h"
#include "image_helper.h"
#include <fstream>
#include <iostream>
#include <math.h>

namespace	dsmc
{

DSMCSolver::DSMCSolver()
{

}

DSMCSolver::~DSMCSolver()
{

}

void	DSMCSolver::Init(const settings_t& s)
{
	m_settings = s;
	InitData();
}

void	DSMCSolver::SetModel(Model* model)
{
	const float*	vertices = &model->GetSurfaces()[0].vertices[0].x;
	const size_t	vertSize = model->GetSurfaces()[0].vertices.size();

	InitRegularGrid(model->GetGrid(), vertices, vertSize);
}

void	DSMCSolver::RunBenchmarkSimulationStep(float dt)
{
	std::cout<<"Running DSMCSolver::RunBenchmarkSimulationStep (base)\n";

	RunSimulationStep(dt);
}

void	DSMCSolver::DoBenchmark(uint itn, float dt)
{
	uint	perc = 0, lastPerc = 0;

	std::cout<<"Done so far:\n0% "<<std::flush;

	for(uint i = 0; i < itn; i++)
	{
		RunBenchmarkSimulationStep(dt);

		perc = (float)i/itn*100;

		if(perc - lastPerc >= 10)
		{
			lastPerc = perc;
			std::cout<<perc<<"% "<<std::flush;
		}
	}

	std::cout<<"100%\n";
}

void	DSMCSolver::RegisterSimState(simulation_state_t* state)
{
	m_state = state;
}

void	DSMCSolver::RenderVectorField()
{
	glEnableClientState(GL_VERTEX_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, vectorFieldBuffer);

	glVertexPointer(3, GL_FLOAT, 0, 0);

	glDrawArrays(GL_LINES, 0, m_settings.grid_dim.x*m_settings.grid_dim.y*2);

	glDisableClientState(GL_VERTEX_ARRAY);
}

void	DSMCSolver::ProcessAndOutputStatistics(dsmc::statistics_t*	stat, float4 collStat, std::ostream& out)
{
	double	tot = 0;
	float3	vel;

	out<<"Selections: "<<collStat.x<<"\nCollisions: "<<collStat.y<<"\nMean separation: "<<collStat.z/collStat.y<<"\nCollision ratio: "<<collStat.y/collStat.x<<"\n";

	out<<"Simulation time: "<<m_state->simTime<<" Running time: "<<m_state->realTime<<"\n";

	out<<"Cell N  ||  SampNum  ||  Dens  ||          Velocity          ||  Temperature\n";

	double	denn = 0;

	double	avg_denn = 0;

	uint	cell_count = m_settings.getCellsCount();

	float	cell_volume = m_settings.cell_size.x*m_settings.cell_size.y*m_settings.cell_size.z;

	for(uint i = 0; i < cell_count; i++)
	{
		dsmc::statistics_t&	s = stat[i];
                           //cell width
		float a = m_settings.fnum/(cell_volume*m_state->samplesCount);
		denn = a*s.numInSample;

		avg_denn += denn;

		vel = make_float3(s.xyzSum.x/s.numInSample, s.xyzSum.y/s.numInSample, s.xyzSum.z/s.numInSample);

		float	uu = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;

		float	tt = m_settings.gas.molmass*(s.velMag/s.numInSample - uu)/(3.0f*DSMC_BOLZ);

		tot += tt;

		out<<i<<"  ||  "<<uint(s.numInSample)<<"  ||  "<<denn<<"  ||  "<<vel.x<<","<<vel.y<<","<<vel.z<<"  ||  "<<tt<<std::endl;
	}

	tot /= cell_count; // avg temp
	avg_denn /= cell_count;

	const float		grid_volume = m_settings.boundaries.get_xsize()*m_settings.boundaries.get_ysize()*m_settings.boundaries.get_zsize();
//	const double	dens2 = m_settings.density*m_settings.density;
//	const double	diam2 = m_settings.gas.diameter*m_settings.gas.diameter;
	const double	dens2 = m_settings.density*m_settings.gas.diameter;
	const double	diam2 = m_settings.gas.diameter*m_settings.density;
	const double	molmas = m_settings.gas.molmass;
	const float		Tref = m_settings.gas.Tref;

	float tColl = 2.0f*m_state->simTime*dens2*grid_volume*DSMC_PI*diam2* \
			pow(tot/Tref,0.25)*sqrtf(DSMC_BOLZ*Tref/(DSMC_PI*molmas))/m_settings.fnum;

	out<<"\n\nExpected collisions: "<<tColl<<" Ratio: "<<tColl/collStat.y<<"\n";

	std::cout<<"Expected collisions: "<<tColl<<" Actual collisions: "<<collStat.y<<"\n";
	std::cout<<"Collision ratio: "<<tColl/collStat.y<<"\n";

	out<<"Average temperature: "<<tot<<" Reference temperature is: "<<m_settings.temp<<"\nAverage density: "<<denn<<" Reference density: "<<m_settings.density<<"\n\n"<<std::endl;
}

void	DSMCSolver::DumpStatistics(std::ostream& out)
{
	float4	collStat;

	statistics_t*	stat = SampleStatistics(collStat);

	ProcessAndOutputStatistics(stat, collStat, out);

	delete [] stat;
}

void	DSMCSolver::SaveColorMap(const	std::string& filename)
{
	uint	w, h, bpp;

	unsigned char*	data = GetColorMap(w, h, bpp);

	save_png(filename.c_str(), data, w, h, bpp);

	delete	[] data;
}

const	uint	FieldWidth  = 128;
const	uint	FieldHeight = 128;

void	DSMCSolver::SaveVelocityField(const std::string& filename)
{
	float3*	velmap = GetVelocityField();

	unsigned char*	data = new unsigned char[FieldWidth*FieldHeight];

	memset(data,255, FieldWidth*FieldHeight);

	uint	ccount = m_settings.grid_dim.x*m_settings.grid_dim.y;

//	std::ofstream	out("vellog.txt");

	for(uint i = 0; i < ccount; i++)
	{
		float3	p1 = velmap[i*2+0];
		float3	p2 = velmap[i*2+1];

//		out<<p1.x<<" "<<p1.y<<" "<<p2.x<<" "<<p2.y<<"\n";

		draw_line(data, FieldWidth, p1.x*FieldWidth,p1.y*FieldHeight,p2.x*FieldWidth,p2.y*FieldHeight, 0);
	}

//	out.close();

	save_png(filename.c_str(), data, FieldWidth, FieldHeight, 1);

	delete	[] data;

	delete [] velmap;
}

DSMCSolver*	CreateSolver(SolverType type)
{
	switch(type)
	{
	case	SolverCPU:
		return	0;
	case	SolverGPUVBO:
		return	new	GPUDSMCSolver<particles_state_vbo_t>();
	case	SolverGPU:
		return	new	GPUDSMCSolver<particles_state_t>();
	}
	return 0;
}

}
