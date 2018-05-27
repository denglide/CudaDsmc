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

#include "config.h"
#include "reggrid.h"
#include "statistics.h"
#include "sim_settings.h"

__host__	void	InitSimulationProperties(const	settings_t& sett);
  
__host__	void	SampleConcentration(float* concentration, const uint2*	cellStartEnd);
__host__	void	SampleConcentrationSliced(float* concentration, const uint slice, const uint2*	cellStartEnd);
__host__	void	SampleVelocitiesSliced(float* velocities, const uint slice, const uint2*	cellStartEnd, const uint* indices);
__host__	void	BuildColorField(float* data, uchar4* colors, uint width, uint height, uint pix_x, uint pix_y);
__host__	void	ComputeBird(float4 *velocities, const uint2*	cellStartEnd, const uint* indices, const uint avgConc, const float dt, const bool mt);
__host__ 	uint*	InitRegularGrid(reg_grid_t& grid, const float*	vertices, uint	vertCount);
__host__ 	void	DeinitRegularGrid();
__host__	void	InitGPUTwisters(const char* fname, unsigned int seed);
__host__	void	CreateVectorField(float3*	field, const float4* velocities, const uint2*	cellStartEnd, const uint* indices, const uint slice);
__host__	size_t	GenerateInitialVelocities(float4 *pos, float4* vel, float3* colors, const float3 streamVector, uint partPerCell, uint numParticles, uint* geom_hashes);
__host__	void	GenInitialVelsUniform(float4 *pos, float4* vel, float3* colors, const float3 streamVector, uint partPerCell, uint numParticles);
__host__	void	IntegrateAndProcessCollisions(float4 *positions, float4 *velocities, float3* colors, float dt, const uint* geom_hashes, const uint numParticles, const bool col, char bc='n');
__host__	void	Integrate(float4 *positions, float4 *velocities, float3* colors, float dt, const uint numParticles, const bool col, char bc='n');
__host__	uint	CreateRandomSeeds(uint gDim, uint bDim, uint** seeds, uint&);
__host__	uint	CreateRandomSeedsUint4(uint4** seeds, uint size);

__host__	void	InitMersenneTwisters();

__host__	void	InitBirdData(uint	cellsCount, uint partPerCells);
__host__	void	DeleteBirdData();

__host__	float4	GetStatistics();

__host__	dsmc::statistics_t*	SampleCellStatistics(const float4* velocities, const uint2* cellStartEnd,const uint* indices, bool readBack = false);

__host__	void	PreProcessVectorField(float3*	field, float width, float height, float minx, float miny);

__host__ void DumpBirdData(const uint2* cellStartEnd, const float avgConc, const float dt, const uint cellsCount);

__host__	void	SortParticlesIndices(uint*	hashes, uint*	indices, uint size);

__host__	void	CalculateCellIDs(float4* positions, uint* d_hashes, uint* d_indices, uint numParticles);

__host__	void	FindCellBoundaries(uint2*   cellStartEnd,
													 float4* sorted_pos,        // output: sorted positions
													 float4* sorted_vel,
													 float3* sorted_cols,

													 uint *  gridParticleHash, // input: sorted grid hashes
													 uint *  gridParticleIndex,// input: sorted particle indices
													 float4* old_pos,            // input: sorted position array
													 float4* old_vel,
													 float3* old_cols,
													 uint numParticles, const bool col);

//TODO: move to cuda helper
bool checkCUDAError(const char *msg);
void checkCUDAErrorAndThrow(const char *msg);

//return size of seeds


