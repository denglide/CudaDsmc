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

#include "dsmc_kernel.cuh"
#include "dsmc.h"
#include "MersenneTwister.h"
#include <stdexcept>
#include <cuda_gl_interop.h>

namespace	dsmc
{

template	<class	DataStrategy>
GPUDSMCSolver<DataStrategy>::GPUDSMCSolver(): gridHashes(0), cellGridIndices(0), cellStartEnd(0), geomHashes(0)
{

}

template	<class	DataStrategy>
GPUDSMCSolver<DataStrategy>::~GPUDSMCSolver()
{
	if(gridHashes)
		cudaFree(gridHashes);

	if(cellGridIndices)
		cudaFree(cellGridIndices);

	if(cellStartEnd)
		cudaFree(cellStartEnd);

	if(geomHashes)
		cudaFree(geomHashes);

	if(concentration)
		cudaFree(concentration);

	if(!m_settings.benchmark)
	{
		cudaGLUnregisterBufferObject(vectorFieldBuffer);
		glDeleteBuffers(1, &vectorFieldBuffer);
		colorMap.Deinit();
	}

	DeinitRegularGrid();
	DeleteBirdData();
}

template	<class	DataStrategy>
uint	GPUDSMCSolver<DataStrategy>::InitSystemState(const float3& streamVel)
{
	float4*	pos = particles.lock_positions();
	float3*	cols = particles.lock_colors();

	if(m_settings.geometry)
	{
		m_settings.numParticles = GenerateInitialVelocities(pos, particles.velocities, cols, streamVel, m_settings.partPerCell,
																m_settings.numParticles,geomHashes);
		checkCUDAError("GenerateInitialVelocities");
		IntegrateAndProcessCollisions(pos, particles.velocities, cols, 0, geomHashes, m_settings.numParticles, !m_settings.benchmark,m_settings.periodicCondition);
		checkCUDAError("IntegrateAndProcessCollisions");
	}
	else
	{
		GenInitialVelsUniform(pos, particles.velocities, cols, streamVel, m_settings.partPerCell, m_settings.numParticles);
		Integrate(pos, particles.velocities, cols, 0, m_settings.numParticles, !m_settings.benchmark,m_settings.periodicCondition);
	}

	particles.unlock_positions();
	particles.unlock_colors();

	return	m_settings.numParticles;
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::UpdateSortedGrid(float4* pos)
{
	CalculateCellIDs(pos, gridHashes, cellGridIndices, m_settings.numParticles);

	if(!checkCUDAError("CalculateCellIDs"))
		throw	std::runtime_error("CUDA error");

	SortParticlesIndices(gridHashes, cellGridIndices, m_settings.numParticles);

	cudaMemset(cellStartEnd, 0, m_settings.getCellsCount()*sizeof(uint2));
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::ReorderParticles(float4* pos, float4* vel, float3* col,
						float4*	old_pos, float4* old_vel, float3* old_col)
{
    FindCellBoundaries(
				cellStartEnd,      // output: cell start index
				pos,
				vel,
				col,

				gridHashes,		// input: sorted grid hashes
				cellGridIndices,	// input: sorted particle indices
				old_pos,       // input: sorted position array
				old_vel,
				old_col,
				m_settings.numParticles, !m_settings.benchmark);

	if(!checkCUDAError("ReorderDataAndFindCellStart"))
		throw	std::runtime_error("CUDA error");
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::RunSimulationStep(float dt)
{
	m_state->advanceSimStep();

#ifdef	REORDER
	particles.swap_positions();
	particles.swap_colors();
	particles.swap_velocities();

	float4* pos 		= particles.lock_positions();
	float4* old_pos 	= particles.lock_old_positions();
	float3*	colors 		= particles.lock_colors();
	float3*	old_colors  = particles.lock_old_colors();

	UpdateSortedGrid(old_pos);

	ReorderParticles(pos, particles.velocities, colors, old_pos, particles.old_velocities, old_colors);

#else
	float4* pos 		= particles.lock_positions();
	float3*	colors 		= particles.lock_colors();

	UpdateSortedGrid(pos);

	ReorderParticles(0, 0, 0, 0, 0, 0);

#endif

	if(!checkCUDAError("UpdateGrid"))
		throw	std::runtime_error("CUDA error");

	if(m_state->drawVectorField)
	{
		float3*	vf = 0;
		cudaGLMapBufferObject((void**)&vf, vectorFieldBuffer);

		CreateVectorField(vf, particles.velocities, cellStartEnd,cellGridIndices, m_state->slice);

		cudaGLUnmapBufferObject(vectorFieldBuffer);
	}

	if(m_state->processCollisions)
	{
		ComputeBird(particles.velocities,cellStartEnd,cellGridIndices,m_settings.partPerCell, m_settings.dt, m_settings.mt);

		if(!checkCUDAError("ComputeBird"))
			throw	std::runtime_error("CUDA error");
	}

	if(m_settings.geometry)
		IntegrateAndProcessCollisions(pos, particles.velocities, colors, dt, geomHashes, m_settings.numParticles, true,m_settings.periodicCondition);
	else
		Integrate(pos, particles.velocities, colors, dt, m_settings.numParticles, true,m_settings.periodicCondition);

	if(!checkCUDAError("IntegrateColors3D"))
		throw	std::runtime_error("CUDA error");

	//TODO: move 100 somewhere
	if(!(m_state->stepCount%100))
	{
		SampleCellStatistics(particles.velocities, cellStartEnd,cellGridIndices);
		m_state->advanceSamplesCount();

		if(!checkCUDAError("SampleCellStatistics"))
			throw	std::runtime_error("CUDA error");
	}

	particles.unlock_positions();
	particles.unlock_colors();

#ifdef	REORDER
	particles.unlock_old_positions();
	particles.unlock_old_colors();
#endif

	m_state->advanceSimTime(dt);

	if(m_state->sampleConcentration)
	{
		SampleConcentrationSliced(concentration, m_state->slice, cellStartEnd);

		uchar4* cols11;
		cudaGLMapBufferObject((void**)&cols11, colorMap.GetDataProvider().colorBuffer);

		BuildColorField(concentration,cols11,32,32,1,1);

		cudaGLUnmapBufferObject(colorMap.GetDataProvider().colorBuffer);
	}
}

//TODO: Get rid of colors
template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::RunBenchmarkSimulationStep(float dt)
{
	m_state->advanceSimStep();

#ifdef	REORDER
	particles.swap_positions();
	particles.swap_colors();
	particles.swap_velocities();

	float4* pos 		= particles.lock_positions();
	float4* old_pos 	= particles.lock_old_positions();
	float3*	colors 		= particles.lock_colors();
	float3*	old_colors  = particles.lock_old_colors();

	UpdateSortedGrid(old_pos);

	ReorderParticles(pos, particles.velocities, colors, old_pos, particles.old_velocities, old_colors);

#else
	float4* pos 		= particles.lock_positions();
	float3*	colors 		= particles.lock_colors();

	UpdateSortedGrid(pos);

	ReorderParticles(0, 0, 0, 0, 0, 0);

#endif

	if(m_state->processCollisions)
		ComputeBird(particles.velocities,cellStartEnd,cellGridIndices,m_settings.partPerCell, m_settings.dt, m_settings.mt);

	if(m_settings.geometry)
		IntegrateAndProcessCollisions(pos, particles.velocities, colors, dt, geomHashes, m_settings.numParticles, false,m_settings.periodicCondition);
	else
		Integrate(pos, particles.velocities, colors, dt, m_settings.numParticles, false, m_settings.periodicCondition);

	//TODO: move 100 somewhere
	if(!(m_state->stepCount%100))
	{
		SampleCellStatistics(particles.velocities, cellStartEnd,cellGridIndices);
		m_state->advanceSamplesCount();
	}

	particles.unlock_positions();
	particles.unlock_colors();

#ifdef	REORDER
	particles.unlock_old_positions();
	particles.unlock_old_colors();
#endif


	m_state->advanceSimTime(dt);
}

template	<class	DataStrategy>
statistics_t*	GPUDSMCSolver<DataStrategy>::SampleStatistics(float4& collStat)
{
	collStat = GetStatistics();

	return	SampleCellStatistics(particles.velocities, cellStartEnd,cellGridIndices,true);
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::InitRegularGrid(reg_grid_t& grid, const float* v, uint size)
{
	geomHashes = ::InitRegularGrid(grid, v, size);
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::RenderConcMap(uint	x, uint	y, uint	w, uint	h)
{
	colorMap.Render(x,y,w,h);
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::Render()
{
	particles.render(m_settings.numParticles);
}

template	<class	DataStrategy>
unsigned char*	GPUDSMCSolver<DataStrategy>::GetColorMap(uint& w, uint& h, uint&	bpp)
{
	uchar4* cols11;
	cudaGLMapBufferObject((void**)&cols11, colorMap.GetDataProvider().colorBuffer);

	unsigned	char*	data = new	unsigned char[4*32*32];

	cudaMemcpy(data, cols11, 4*32*32*sizeof(unsigned char), cudaMemcpyDeviceToHost);

	cudaGLUnmapBufferObject(colorMap.GetDataProvider().colorBuffer);

	w = 32;
	h = 32;
	bpp = 4;

	return data;
}

template	<class	DataStrategy>
float3*		GPUDSMCSolver<DataStrategy>::GetVelocityField()
{
	float3*	vf = 0;
	cudaGLMapBufferObject((void**)&vf, vectorFieldBuffer);

	PreProcessVectorField(vf,m_settings.boundaries.get_xsize(),m_settings.boundaries.get_ysize(),
								m_settings.boundaries.min.x, m_settings.boundaries.min.y);

	uint	size = m_settings.grid_dim.x*m_settings.grid_dim.x*2;

	float3*	data = new float3[size];

	cudaMemcpy(data, vf, size*sizeof(float3), cudaMemcpyDeviceToHost);

	cudaGLUnmapBufferObject(vectorFieldBuffer);

	return data;
}

template	<class	DataStrategy>
void	GPUDSMCSolver<DataStrategy>::InitData()
{
	const char *dat_path = "data/MersenneTwister.dat";

	InitGPUTwisters(dat_path, SEED);

	checkCUDAError("InitGPUTwisters");

	particles.init(m_settings.numParticles);

	cudaMalloc((void **)&gridHashes, m_settings.numParticles * sizeof(uint));
	cudaMalloc((void **)&cellGridIndices, m_settings.numParticles * sizeof(uint));

	const	uint	cellsCount = m_settings.getCellsCount();

	cudaMalloc((void **)&concentration, m_settings.grid_dim.x*m_settings.grid_dim.y * sizeof(float));

	cudaMalloc((void **)&cellStartEnd, cellsCount * sizeof(uint2));

	InitSimulationProperties(m_settings);
	checkCUDAError("InitSimulationProperties");

	InitBirdData(cellsCount, m_settings.partPerCell);
	checkCUDAError("InitBirdData");

	InitMersenneTwisters();
	checkCUDAError("InitMersenneTwisters");

	if(!m_settings.benchmark)
	{
		uint vbsz = m_settings.grid_dim.x*m_settings.grid_dim.y*2*sizeof(float3);

		vectorFieldBuffer = CreateCUDABufferObject(vbsz);

		colorMap.Init(32,32);
	}
}
}
