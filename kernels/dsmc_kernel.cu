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
#include "dsmc.h"

#include <stdio.h>
#include <vector_types.h>
#include "dsmc_kernel.cuh"
#include "WarpStandard.cuh"

#ifndef _MSC_VER

#include <fcntl.h>
#include <unistd.h>

#else

#include <windows.h>

#endif

#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/remove.h>
#include <thrust/device_vector.h>

#include <vector>

#include <fstream>

using namespace std;

#define BOUNCE_DELTA 0.001

//#define	NVIDIA_RADIX_SORT

#ifdef	NVIDIA_RADIX_SORT

#include "radixsort.h"

nvRadixSort::RadixSort*	radixSort = 0;

uint	CalculateSortingBits(uint gridCells)
{
	int	i = -1;
	while(gridCells)
	{
		gridCells>>=1;
		i++;
	}

	return i;
}

uint	g_sortBits = 0;

#endif

struct	gpu_sim_settings_t
{
	uint	numParticles;
	uint	partPerCell;

	//molmass,diam, temp
	float3	gas_props;
	//density, fnum, temp
	float3	phys_state;

	float3	grid_min;
	float3	grid_max;

	float3	cell_size;
	uint3	grid_dim;

	float	cell_volume;
	uint	cells_count;

	float	crossSection;

	float	maxRelVel;
};

struct	sim_grid_t
{
	uint3	dim;
	uint	cells_count;
};

struct	collision_grid_t
{
	float3	mmin;
	float3	mmax;
	float3	cellSize;
	uint3	dim;
};

sim_grid_t	g_hostSimGrid;

uint3*	devIndices = 0;
uint2*	devIndicesStartEnd = 0;
float3*	devVertices = 0;
float3*	devFaceNormals = 0;

uint4*	devSeeds = 0;

float4*	devStatistics = 0;

dsmc::statistics_t*	devCellStatistics = 0;

uint*	counts = 0;
float2*	maximums = 0;
float*	devConcentration = 0;
uchar2*	pairs = 0;

__device__ __constant__	collision_grid_t	g_collGrid;

__device__ __constant__ gpu_sim_settings_t	g_simSettings;

//TODO: Fix this (add headers)
#include "mtwister.cu"
#include "tri_intersect.cu"
#include "vectors.cu"

uint	GetNumberOfBlocks(uint threads, uint numParticles)
{
	uint	numBlocks = numParticles/threads;

	if(numParticles%threads)
		numBlocks++;

	return numBlocks;
}

__device__	uint	GetCollisionCellID(float3 pos)
{
	if(pos.x < g_collGrid.mmin.x ||
		pos.y < g_collGrid.mmin.y ||
		pos.z < g_collGrid.mmin.z ||

		pos.x > g_collGrid.mmax.x ||
		pos.y > g_collGrid.mmax.y ||
		pos.z > g_collGrid.mmax.z
	)

	return 0xffffffff;

	uint x = (pos.x - g_collGrid.mmin.x) / g_collGrid.cellSize.x;
	uint y = (pos.y - g_collGrid.mmin.y) / g_collGrid.cellSize.y;
	uint z = (pos.z - g_collGrid.mmin.z) / g_collGrid.cellSize.z;

	x = min(x, g_collGrid.dim.x-1);
	y = min(y, g_collGrid.dim.y-1);
	z = min(z, g_collGrid.dim.z-1);

	return z*g_collGrid.dim.x*g_collGrid.dim.y +
			y*g_collGrid.dim.x +
			x;
}

__device__	uint	GetCollisionCellID(float4 pos)
{
	if(pos.x < g_collGrid.mmin.x ||
		pos.y < g_collGrid.mmin.y ||
		pos.z < g_collGrid.mmin.z ||

		pos.x > g_collGrid.mmax.x ||
		pos.y > g_collGrid.mmax.y ||
		pos.z > g_collGrid.mmax.z
	)

	return 0xffffffff;

	uint x = (pos.x - g_collGrid.mmin.x) / g_collGrid.cellSize.x;
	uint y = (pos.y - g_collGrid.mmin.y) / g_collGrid.cellSize.y;
	uint z = (pos.z - g_collGrid.mmin.z) / g_collGrid.cellSize.z;

	x = min(x, g_collGrid.dim.x-1);
	y = min(y, g_collGrid.dim.y-1);
	z = min(z, g_collGrid.dim.z-1);

	return z*g_collGrid.dim.x*g_collGrid.dim.y +
			y*g_collGrid.dim.x +
			x;
}

__host__	uint*	UploadGeometryHashes(reg_grid_t& grid)
{
	uint*	host_has_geom = (uint*)malloc(grid.cellCount()*sizeof(uint));

	for(uint i = 0; i < grid.cellCount(); i++)
		host_has_geom[i] = grid[i].second ? 1 : 0;

	uint*	geom_hashes = 0;

	cudaError_t err = cudaMalloc(&geom_hashes, grid.cellCount()*sizeof(uint));
	err = cudaMemcpy(geom_hashes,host_has_geom,grid.cellCount()*sizeof(uint), cudaMemcpyHostToDevice);

	free(host_has_geom);

	return geom_hashes;
}

void	UploadIndices(reg_grid_t& grid)
{
	uint	idxCountTotal = 0;
	for(uint i = 0; i < grid.cellCount(); i++)
		idxCountTotal += grid[i].first.size();

	uint3*	host_indices = (uint3*)malloc(idxCountTotal*sizeof(uint));
	uint2*	host_cellStEnd = (uint2*)malloc(grid.cellCount()*sizeof(uint2));

	float3*	host_normals = (float3*)malloc(idxCountTotal/3*sizeof(float3));

	memset(host_cellStEnd,0xff,grid.cellCount()*sizeof(uint2));

	uint offset = 0;
	for(uint i = 0; i < grid.cellCount(); i++)
	{
		uint	idxCount = grid[i].first.size();

		if(idxCount)
		{
			reg_grid_t::idx_vector_t&	node = grid[i].first;
			host_cellStEnd[i].x = offset;
			for(uint k = 0; k < idxCount; k += 3)
			{
				host_indices[offset + k/3].x = node[k+0];
				host_indices[offset + k/3].y = node[k+1];
				host_indices[offset + k/3].z = node[k+2];

				vec3_t& norm = grid.normals[i][k/3];

				host_normals[offset + k/3].x = norm.x;
				host_normals[offset + k/3].y = norm.y;
				host_normals[offset + k/3].z = norm.z;
			}
			offset += idxCount/3;
			host_cellStEnd[i].y = offset;
		}
	}

	cudaMalloc(&devIndices, idxCountTotal*sizeof(uint));
	cudaMemcpy(devIndices,host_indices,idxCountTotal*sizeof(uint), cudaMemcpyHostToDevice);

	cudaMalloc(&devIndicesStartEnd, grid.cellCount()*sizeof(uint2));
	cudaMemcpy(devIndicesStartEnd,host_cellStEnd,grid.cellCount()*sizeof(uint2), cudaMemcpyHostToDevice);

	cudaMalloc(&devFaceNormals, idxCountTotal/3*sizeof(float3));
	cudaMemcpy(devFaceNormals, host_normals, idxCountTotal/3*sizeof(float3), cudaMemcpyHostToDevice);

	free(host_indices);
	free(host_cellStEnd);
	free(host_normals);
}

__global__	void	DumpConstantsDev(float3*	dump)
{
	dump[0] = g_simSettings.gas_props;
	dump[1] = g_simSettings.phys_state;
	dump[2] = g_simSettings.grid_min;
	dump[3] = g_simSettings.grid_max;
	dump[4] = g_simSettings.cell_size;
	dump[5] = make_float3(g_simSettings.grid_dim.x,g_simSettings.grid_dim.y,g_simSettings.grid_dim.z);
	dump[6] = make_float3(g_simSettings.cell_volume, g_simSettings.cells_count, g_simSettings.maxRelVel);
}

void	DumpConstants()
{
	printf("Dumping constants\n");

	float3*	dump = 0;
	float3*	cpu_dump = new float3[7];

	cudaMalloc(&dump, 7*sizeof(float3));
	DumpConstantsDev<<<1,1>>>(dump);

	cudaMemcpy(cpu_dump,dump,7*sizeof(float3),cudaMemcpyDeviceToHost);

	for(uint i = 0; i < 7; i++)
	{
		printf("%e %e %e\n",cpu_dump[i].x,cpu_dump[i].y,cpu_dump[i].z);
	}

	cudaFree(dump);

	delete [] cpu_dump;
}

__host__	void	InitSimulationProperties(const	settings_t& sett)
{
	gpu_sim_settings_t	gpu_sett;

	gpu_sett.numParticles = sett.numParticles;
	gpu_sett.partPerCell = sett.partPerCell;
	gpu_sett.gas_props = make_float3(sett.gas.molmass, sett.gas.diameter, sett.gas.Tref);
	gpu_sett.phys_state = make_float3(sett.density, sett.fnum, sett.temp);
	gpu_sett.grid_min = make_float3(sett.boundaries.min.x,sett.boundaries.min.y, sett.boundaries.min.z);
	gpu_sett.grid_max = make_float3(sett.boundaries.max.x,sett.boundaries.max.y, sett.boundaries.max.z);
	gpu_sett.cell_size = make_float3(sett.cell_size.x,sett.cell_size.y,sett.cell_size.z);
	gpu_sett.grid_dim = make_uint3(sett.grid_dim.x, sett.grid_dim.y, sett.grid_dim.z);
	gpu_sett.cells_count = sett.getCellsCount();

	double	crossSect = DSMC_PI*sett.gas.diameter*sett.gas.diameter;

	gpu_sett.crossSection = crossSect;

	gpu_sett.cell_volume = sett.cell_size.x*sett.cell_size.y*sett.cell_size.z;

	double	rv = sett.gas.diameter*sett.gas.diameter*DSMC_PI*sqrtf(sett.temp/300.0)*300;

	gpu_sett.maxRelVel = rv;

	if(cudaSuccess != cudaMemcpyToSymbol(g_simSettings, &gpu_sett, sizeof(gpu_sett)))
	{
		printf("\nCan't init simulation settings\n");
	}

	g_hostSimGrid.cells_count = sett.getCellsCount();
	g_hostSimGrid.dim 		  = make_uint3(sett.grid_dim.x,sett.grid_dim.y,sett.grid_dim.z);

#ifdef	NVIDIA_RADIX_SORT

	g_sortBits = CalculateSortingBits(sett.getCellsCount());

#endif

//	DumpConstants();
}

__host__ uint*	InitRegularGrid(reg_grid_t& grid, const float*	vertices, uint	vertCount)
{
	collision_grid_t* g;

	g = (collision_grid_t*)malloc(sizeof(collision_grid_t));

	g->mmin  = make_float3(grid.gridSize.min.x,grid.gridSize.min.y,grid.gridSize.min.z);
	g->mmax  = make_float3(grid.gridSize.max.x,grid.gridSize.max.y,grid.gridSize.max.z);
	g->dim  = make_uint3(grid.dim.x,grid.dim.y,grid.dim.z);
	g->cellSize = make_float3(grid.cellSize.x,grid.cellSize.y, grid.cellSize.z);

	if(cudaSuccess != cudaMemcpyToSymbol(g_collGrid, g, sizeof(collision_grid_t)))
	{
		printf("\nCan't init collision grid\n");
	}

	uint* gh = UploadGeometryHashes(grid);

	UploadIndices(grid);

	free(g);

	cudaMalloc(&devVertices, vertCount*sizeof(float3));
	cudaMemcpy(devVertices,vertices,vertCount*sizeof(float3), cudaMemcpyHostToDevice);

	return gh;
}

__host__ void	DeinitRegularGrid()
{
	if(devIndices)
		cudaFree(devIndices);

	if(devIndicesStartEnd)
		cudaFree(devIndicesStartEnd);

	if(devVertices)
		cudaFree(devVertices);

	if(devFaceNormals)
		cudaFree(devFaceNormals);
}

//gives VMP
__device__	float	SampleDensityFunction(unsigned	int	id, mt_state* mts, mt_struct_stripped* config, float VMP)
{
//according to Bird's code
	return sqrtf(2*DSMC_BOLZ*g_simSettings.phys_state.z/g_simSettings.gas_props.x*(-log(GetRandomFloatFast(id, mts, config))));//*GetRandomFloatFast(id, mts, config);
}

__device__ float4	GetRandomVelocity(int tid, mt_state* mts, mt_struct_stripped* config, float3 s)
{
	float	vv = SampleDensityFunction(tid, mts, config, g_simSettings.phys_state.z);

	const float	tetha = 2.0f*DSMC_PI;

	float4	v;

	v.x = vv*sinf(GetRandomFloatFast(tid, mts, config)*tetha) + s.x;
	v.y = vv*sinf(GetRandomFloatFast(tid, mts, config)*tetha) + s.y;
	v.z = vv*sinf(GetRandomFloatFast(tid, mts, config)*tetha) + s.z;

	v.w = 0;

	return	v;
}

__device__ float4	GetRandomPosition(int tid, mt_state* mts, mt_struct_stripped* config, float3 cell)
{
	float4	v;

	v.x = cell.x + GetRandomFloatFast(tid, mts, config)*g_simSettings.cell_size.x;
	v.y = cell.y + GetRandomFloatFast(tid, mts, config)*g_simSettings.cell_size.y;
	v.z = cell.z + GetRandomFloatFast(tid, mts, config)*g_simSettings.cell_size.z;
	v.w = 1;

	return	v;
}
//add slice
__global__	void	CreateVectorFieldDev(float3*	field, const float4* velocities, const uint2*	cellStartEnd,const uint* indices,const uint slice)
{
	const	uint	cellid = gridDim.x*blockDim.x*slice + threadIdx.x*gridDim.x + blockIdx.x;

	const float3	cell = make_float3(g_simSettings.grid_min.x + blockIdx.x*g_simSettings.cell_size.x,
									g_simSettings.grid_min.y + threadIdx.x*g_simSettings.cell_size.y,
									g_simSettings.grid_min.z + slice*g_simSettings.cell_size.z);

	uint2	cellStEnd = cellStartEnd[cellid];

	const uint	part_in_cell = cellStEnd.y - cellStEnd.x;

	uint	index = 0;

	float4	vel;

	const uint	cstart = cellStEnd.x;

	float3	p1 = make_float3(cell.x + g_simSettings.cell_size.x/2.0f, cell.y + g_simSettings.cell_size.y/2.0f, cell.z + g_simSettings.cell_size.z/2.0f);

	float3	p2 = make_float3(0,0,0);

	for(uint i = 0; i < part_in_cell; i++)
	{
		index = cstart + i;

#ifndef	REORDER
		index = indices[index];
#endif

		vel = velocities[index];

		p2.x += vel.x;
		p2.y += vel.y;
		p2.z += vel.z;
	}

	p2 = normalize(p2);

	uint	idx = threadIdx.x*gridDim.x + blockIdx.x;

	field[idx*2+0] = p1;
	field[idx*2+1] = make_float3(p1.x + p2.x*g_simSettings.cell_size.x/2.0f, p1.y + p2.y*g_simSettings.cell_size.y/2.0f, p1.z + p2.z*g_simSettings.cell_size.z/2.0f);
}

__host__	void	CreateVectorField(float3*	field, const float4* velocities, const uint2*	cellStartEnd,const uint* indices,const uint slice)
{
	CreateVectorFieldDev<<<g_hostSimGrid.dim.x, g_hostSimGrid.dim.y>>>(field, velocities, cellStartEnd,indices,slice);
}

__device__ bool IsGeometryInCell(float3 cmin, float3 csize, uint* geom_hashes)
{
	float3	p = cmin;
	uint cid = 0xffffffff;
	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p.x += csize.x;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p.y += csize.y;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p.z += csize.z;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p = cmin;

	p.y += csize.y;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p.z += csize.z;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}

	p.y = cmin.y;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}
	p.x += csize.x;

	if((cid = GetCollisionCellID(p)) != 0xffffffff)
	{
		if(geom_hashes[cid] != 0)
			return true;
	}
	return false;
}

__global__	void	InitMersenneTwistersDev(uint	count)
{
	uint idx = blockDim.x*blockIdx.x + threadIdx.x;

	if(idx >= count)
		return;

	mt_state	mts =	ds_mtState[idx];

	mt_struct_stripped config = ds_MT[idx];

	InitTwisterFast(idx, &mts, &config);

	ds_mtState[idx] = mts;
}

__host__	void	InitMersenneTwisters()
{
	uint	numThreads = 128;
	uint	numBlocks = GetNumberOfBlocks(numThreads, MT_RNG_COUNT);

	InitMersenneTwistersDev<<<numBlocks, numThreads>>>(MT_RNG_COUNT);
}

#define INITIAL_COLOR	0,1,0

__global__	void	GenerateUniformVelocitiesDevBlocks(float4 *pos, float4* vel, float3* colors, const float3 streamVector, uint partPerCell, uint z_dim)
{
	const	uint cellX = blockIdx.x;
	const	uint cellY = blockIdx.y;

	//unique tid for MersenneTwister
	const	uint	tid = gridDim.x*gridDim.y*threadIdx.x + blockIdx.y*gridDim.x + blockIdx.x;

	mt_state	mts =	ds_mtState[tid];

	mt_struct_stripped config = ds_MT[tid];

	float3	cell = make_float3(g_simSettings.grid_min.x + cellX*g_simSettings.cell_size.x,
									g_simSettings.grid_min.y + cellY*g_simSettings.cell_size.y,
									0);

	float3	col = make_float3(0,0,0);

	float4	v;

	for(uint cellZ = 0; cellZ < z_dim; cellZ++)
	{
		uint idx = (gridDim.x*gridDim.y*cellZ + blockIdx.y*gridDim.x + blockIdx.x)*partPerCell;

		cell.z = g_simSettings.grid_min.z + cellZ*g_simSettings.cell_size.z;

		for(uint i = threadIdx.x; i < partPerCell; i += blockDim.x)
		{
			v = GetRandomVelocity(tid, &mts, &config, streamVector);

			vel[idx+i] = v;
			pos[idx+i] = GetRandomPosition(tid, &mts, &config, cell);

			float	l = length(make_float3(v.x,v.y,v.z));

			col.x = l;
			col.y = 1-l;

			colors[idx+i] = col;
		}
	}

	ds_mtState[tid] = mts;
}

__host__	void	GenInitialVelsUniform(float4 *pos, float4* vel, float3* colors, const float3 streamVector, uint partPerCell, uint numParticles)
{
	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.x;
	dimGrid.y = g_hostSimGrid.dim.y;
	dimGrid.z = 1;

	GenerateUniformVelocitiesDevBlocks<<<dimGrid, partPerCell>>>(pos, vel, colors, streamVector, partPerCell, g_hostSimGrid.dim.z);
}

__global__	void	GenerateInitialVelocitiesDevBlocks(float4 *pos, float4* vel, float3* colors, const float3 streamVector, uint partPerCell, uint* geom_hashes, uint z_dim)
{
	const	uint cellX = blockIdx.x;
	const	uint cellY = blockIdx.y;

	//unique tid for MersenneTwister
	const	uint	tid = gridDim.x*gridDim.y*threadIdx.x + blockIdx.y*gridDim.x + blockIdx.x;

	mt_state	mts =	ds_mtState[tid];

	mt_struct_stripped config = ds_MT[tid];

	float3	cell = make_float3(g_simSettings.grid_min.x + cellX*g_simSettings.cell_size.x,
									g_simSettings.grid_min.y + cellY*g_simSettings.cell_size.y,
									0);

	float4	v;

	for(uint cellZ = 0; cellZ < z_dim; cellZ++)
	{
		uint idx = (gridDim.x*gridDim.y*cellZ + blockIdx.y*gridDim.x + blockIdx.x)*partPerCell;

		cell.z = g_simSettings.grid_min.z + cellZ*g_simSettings.cell_size.z;

		bool k = IsGeometryInCell(cell, make_float3(g_simSettings.cell_size.x, g_simSettings.cell_size.y, g_simSettings.cell_size.z), geom_hashes);

		if(!k)
		{
			for(uint i = threadIdx.x; i < partPerCell; i += blockDim.x)
			{
				v = GetRandomVelocity(tid, &mts, &config, streamVector);

				vel[idx+i] = v;
				pos[idx+i] = GetRandomPosition(tid, &mts, &config, cell);

				colors[idx+i] = make_float3(INITIAL_COLOR);
			}
		}
		else
		{
			for(uint i = threadIdx.x; i < partPerCell; i += blockDim.x)
			{
				vel[idx+i] = make_float4(0,0,0,0);//7777777,7777777,7777777,7777777);
				pos[idx+i] = make_float4(0,0,0,0);//7777777,7777777,7777777,7777777);

				colors[idx+i] = make_float3(INITIAL_COLOR);
			}
		}
	}

	ds_mtState[tid] = mts;
}

__global__	void	SampleConcentrationDev(float* concentration, const uint2*	cellStartEnd)
{
	const	uint	cellid = gridDim.x*gridDim.y*threadIdx.x + blockIdx.y*gridDim.x + blockIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];

	concentration[cellid] = cellStEnd.y - cellStEnd.x;
}

__global__	void	SampleConcentrationSlicedDev(float* concentration, uint slice,  const uint2*	cellStartEnd)
{
	const	uint	cellid = gridDim.x*blockDim.x*slice + threadIdx.x*gridDim.x + blockIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];

	concentration[threadIdx.x*gridDim.x + blockIdx.x] = cellStEnd.y - cellStEnd.x;
}

__host__	void	SampleConcentrationSliced(float* concentration, const uint slice, const uint2*	cellStartEnd)
{
	SampleConcentrationSlicedDev<<<g_hostSimGrid.dim.x, g_hostSimGrid.dim.y>>>(concentration,slice,cellStartEnd);
}

__host__	void	SampleConcentration(float* concentration, const uint2*	cellStartEnd)
{
	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.x;
	dimGrid.y = g_hostSimGrid.dim.y;
	dimGrid.z = 1;

	SampleConcentrationDev<<<dimGrid, g_hostSimGrid.dim.z>>>(concentration,cellStartEnd);
}

__global__	void	SampleVelocitiesSlicedDev(float* velocities, const uint slice, const float4* vels_data, const uint2*	cellStartEnd,const uint* indices)
{
	const	uint	cellid = gridDim.x*blockDim.x*slice + threadIdx.x*gridDim.x + blockIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];

	const uint	part_in_cell = cellStEnd.y - cellStEnd.x;

	if(part_in_cell <= 0)
	{
		velocities[threadIdx.x*gridDim.x + blockIdx.x] = 0;
		return;
	}

	float4	vel,p = make_float4(0,0,0,0);

	for(uint	index = cellStEnd.x; index < cellStEnd.y; index++)
	{
#ifndef	REORDER
		uint	idx = indices[index];
		vel = vels_data[idx];
#else
		vel = vels_data[index];
#endif

		p.x += vel.x;
		p.y += vel.y;
		p.z += vel.z;
	}

	velocities[threadIdx.x*gridDim.x + blockIdx.x] = length4(p) / part_in_cell;
}

__host__	void	SampleVelocitiesSliced(float* velocities, const uint slice, const float4* vels_data, const uint2*	cellStartEnd,const uint* indices)
{
	SampleVelocitiesSlicedDev<<<g_hostSimGrid.dim.x, g_hostSimGrid.dim.y>>>(velocities,slice,vels_data,cellStartEnd, indices);
}

__global__	void	BuildColorFieldDev(float* data, uchar4* colors, float* minmax, uint xx, uint yy)
{
	float	mn = minmax[0];
	float	mx = minmax[1];

	float	median = (mx - mn)/2.0f;

	const uint idx = threadIdx.x*gridDim.x/yy/yy + blockIdx.x/xx;

	float	val = data[idx];

	uchar4	col;

#if	1

	if(val < median)
	{
		float alpha = (val - mn)/(median - mn);

		col.x = 0;
		col.y = 255*(1-alpha);
		col.z = 255*alpha;

	}else
	{
		float alpha = (val - median)/(mx - median);

		col.x = 255*alpha;
		col.y = 0;
		col.z = 255*(1-alpha);
	}
#else

	float	alpha = 1;

	if(!(val < 0.1 || mn == mx || mx < 0.1))
		alpha = val/(mx-mn);

	col.x = 255*(1-alpha);
	col.y = 255*(1-alpha);
	col.z = 255*(1-alpha);

#endif
	col.w = 255;

	const	uint	col_idx = threadIdx.x*gridDim.x + blockIdx.x;

	colors[col_idx] = col;

}
//change to thurst
__host__	void	BuildColorField(float* data, uchar4* colors, uint width, uint height, uint pix_x, uint pix_y)
{
	uint size = width*height;

	thrust::device_ptr<float> dev_data(data);
	thrust::device_ptr<float> dev_data_end(data + size);

	thrust::device_vector<float>	output(size);
	thrust::device_vector<float>	minmax(2);

	thrust::inclusive_scan(dev_data, dev_data_end, output.begin(), thrust::minimum<float>());

	minmax[0] = *(output.end()-1);

	//i'm not really sure about this part. Probably we don't need it
	thrust::fill(output.begin(), output.end(), 0.0f);

	thrust::inclusive_scan(dev_data, dev_data_end, output.begin(), thrust::maximum<float>());

	minmax[1] = *(output.end()-1);

	float*	ptr = thrust::raw_pointer_cast(&minmax[0]);
    BuildColorFieldDev<<<width*pix_x, height*pix_y>>>(data, colors, ptr, pix_x, pix_y);
}

namespace	dsmc
{
	struct	float4plus
	{
		__device__	__host__	float4	operator()(float4 a, float4 b)
		{
			return make_float4(a.x+b.x, a.y+b.y, a.z + b.z, a.w + b.w);
		}
	};
}

__host__	float4	GetStatistics()
{
	float4	stat;

	thrust::device_ptr<float4> dev_ptr(devStatistics);
	thrust::device_ptr<float4> dev_ptr_end(devStatistics + g_hostSimGrid.cells_count);

	dsmc::float4plus pf;

	inclusive_scan(dev_ptr, dev_ptr_end, dev_ptr, pf);

	cudaMemcpy(&stat,devStatistics+g_hostSimGrid.cells_count-1, sizeof(float4), cudaMemcpyDeviceToHost);

	return stat;
}

__global__	void	SampleCellStatisticsDev(dsmc::statistics_t*	stats, const float4*	velocities, const uint2* cellStartEnd, const uint* indices)
{
	const	uint	cellid = blockDim.x*gridDim.x*blockIdx.y + blockIdx.x*blockDim.x + threadIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];

	const uint	part_in_cell = cellStEnd.y - cellStEnd.x;

	if(part_in_cell <= 0)
	{
		return;
	}

	dsmc::statistics_t	s =	stats[cellid];

	for(uint	index = cellStEnd.x; index < cellStEnd.y; index++)
	{

#ifndef	REORDER
		uint	idx = indices[index];
		float4	v = velocities[idx];
#else
		float4	v = velocities[index];
#endif
		s.numInSample++;

		s.xyzSum.x += v.x;
		s.xyzSum.y += v.y;
		s.xyzSum.z += v.z;

		s.velMag += v.x*v.x + v.y*v.y + v.z*v.z;
	}

	stats[cellid] = s;
}

__host__	dsmc::statistics_t*	SampleCellStatistics(const float4* velocities, const uint2* cellStartEnd,const uint* indices, bool readBack)
{
	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.y;
	dimGrid.y = g_hostSimGrid.dim.z;
	dimGrid.z = 1;

	SampleCellStatisticsDev<<<dimGrid, g_hostSimGrid.dim.x>>>(devCellStatistics, velocities, cellStartEnd, indices);

	dsmc::statistics_t*	stat = 0;

	if(readBack)
	{
		stat = new dsmc::statistics_t[g_hostSimGrid.cells_count];
		cudaMemcpy(stat, devCellStatistics, g_hostSimGrid.cells_count*sizeof(dsmc::statistics_t), cudaMemcpyDeviceToHost);
	}

	return stat;
}

__global__	void	ComputeBirdCollisionsNoMTDev(float4* velocities, float2* maximums, float4* statistics, uint4* seeds, const uint2* cellStartEnd, const uint* indices, const float avgConc, const float dt)
{
	const	uint	cellid = blockDim.x*gridDim.x*blockIdx.y + blockIdx.x*blockDim.x + threadIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];
	int		partCount = cellStEnd.y - cellStEnd.x;

	if(partCount < 2)
		return;

	uint4	seed = seeds[cellid];

	float4	stat = statistics[cellid];

	float4	v1,v2, v11, v22;

	float2	vmax = maximums[cellid];

	float3	angles;

	uint2	pair;

	float	fpc = (0.5*partCount*avgConc*g_simSettings.phys_state.y*vmax.x*dt)/g_simSettings.cell_volume + vmax.y;

	uint	pairsCount = fpc;

	stat.x += pairsCount;

	vmax.y = fpc - pairsCount;

	for(uint i = 0; i < pairsCount; i++)
	{

		pair.x = cellStEnd.x+HybridTausRngInt(&seed)%partCount;

#ifndef	REORDER
		pair.x = indices[pair.x];
#endif
		v1 = velocities[pair.x];

		pair.y = cellStEnd.x+HybridTausRngInt(&seed)%partCount;

#ifndef	REORDER
		pair.y = indices[pair.y];
#endif

		v2 = velocities[pair.y];

		float	mag = length(make_float3(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z));

		float v = mag*g_simSettings.crossSection;

		if(v/vmax.x > HybridTausRng(&seed))
		{
			stat.y++;
			stat.z += mag;

			float	b = 2.0f*HybridTausRng(&seed) - 1.0f;
			float	a = sqrtf(1-b*b);

			float	c = 2*DSMC_PI*HybridTausRng(&seed);//GetRandomFloatFast(cellid, &mts, &config);

			angles = make_float3(b*mag/2.0f, a*mag*__cosf(c)/2.0f, a*mag*__sinf(c)/2.0f);

			v11.x = (v1.x+v2.x)/2.0f + angles.x;
			v11.y = (v1.y+v2.y)/2.0f + angles.y;
			v11.z = (v1.z+v2.z)/2.0f + angles.z;

			velocities[pair.x] = v11;

			v22.x = (v1.x+v2.x)/2.0f - angles.x;
			v22.y = (v1.y+v2.y)/2.0f - angles.y;
			v22.z = (v1.z+v2.z)/2.0f - angles.z;

			velocities[pair.y] = v22;
		}

		if(v > vmax.x)
			vmax.x = v;
	}

	maximums[cellid] = vmax;

	statistics[cellid] = stat;

	seeds[cellid] = seed;

//	ds_mtState[cellid] = mts;
}
//TODO
__global__	void	ComputeBirdCollisionsDev(float4* velocities, float2* maximums, float4* statistics, const uint2* cellStartEnd, const uint* indices, const float avgConc, const float dt)
{
	const	uint	cellid = blockDim.x*gridDim.x*blockIdx.y + blockIdx.x*blockDim.x + threadIdx.x;

	uint2	cellStEnd = cellStartEnd[cellid];
	int		partCount = cellStEnd.y - cellStEnd.x;

	if(partCount < 2)
		return;

	float4	stat = statistics[cellid];

	float4	v1,v2, v11, v22;

	float2	vmax = maximums[cellid];

	float3	angles;

	uint2	pair;

	float	fpc = (0.5*partCount*avgConc*g_simSettings.phys_state.y*vmax.x*dt)/g_simSettings.cell_volume + vmax.y;

	uint	pairsCount = fpc;

	stat.x += pairsCount;

	vmax.y = fpc - pairsCount;

	mt_state	mts =	ds_mtState[cellid];
	mt_struct_stripped config = ds_MT[cellid];

	for(uint i = 0; i < pairsCount; i++)
	{
		pair.x = cellStEnd.x+GetRandomIntegerFast(cellid, &mts, &config)%partCount;

#ifndef	REORDER
		pair.x = indices[pair.x];
#endif

		v1 = velocities[pair.x];

		pair.y = cellStEnd.x+GetRandomIntegerFast(cellid, &mts, &config)%partCount;

#ifndef	REORDER
		pair.y = indices[pair.y];
#endif

		v2 = velocities[pair.y];

		float	mag = length(make_float3(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z));

		float v = mag*g_simSettings.crossSection;//DSMC_PI*DSMC_DIAM*DSMC_DIAM;

		if(v/vmax.x > GetRandomFloatFast(cellid, &mts, &config))
		{
			stat.y++;
			stat.z += mag;

			float	b = 2.0f*GetRandomFloatFast(cellid, &mts, &config) - 1.0f;
			float	a = sqrtf(1-b*b);

			float	c = 2*DSMC_PI*GetRandomFloatFast(cellid, &mts, &config);

			angles = make_float3(b*mag/2.0f, a*mag*__cosf(c)/2.0f, a*mag*__sinf(c)/2.0f);

			v11.x = (v1.x+v2.x)/2.0f + angles.x;
			v11.y = (v1.y+v2.y)/2.0f + angles.y;
			v11.z = (v1.z+v2.z)/2.0f + angles.z;

			velocities[pair.x] = v11;

			v22.x = (v1.x+v2.x)/2.0f - angles.x;
			v22.y = (v1.y+v2.y)/2.0f - angles.y;
			v22.z = (v1.z+v2.z)/2.0f - angles.z;

			velocities[pair.y] = v22;
		}

		if(v > vmax.x)
			vmax.x = v;
	}

	maximums[cellid] = vmax;

	statistics[cellid] = stat;

	ds_mtState[cellid] = mts;
}

void	ComputeBird(float4 *velocities, const uint2*	cellStartEnd, const uint* indices, const uint avgConc, const float dt, const bool mt)
{
	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.y;
	dimGrid.y = g_hostSimGrid.dim.z;
	dimGrid.z = 1;

	if(mt)
		ComputeBirdCollisionsDev<<<dimGrid, g_hostSimGrid.dim.x>>>(velocities, maximums, devStatistics, cellStartEnd, indices, avgConc, dt);
	else
		ComputeBirdCollisionsNoMTDev<<<dimGrid, g_hostSimGrid.dim.x>>>(velocities, maximums, devStatistics, devSeeds, cellStartEnd,indices,avgConc, dt);
}

__global__	void	InitMaximums(float2* maximum)
{
	const	uint	cellid = gridDim.x*gridDim.y*threadIdx.x + blockIdx.y*gridDim.x + blockIdx.x;

	float2	v;

	mt_state	mts =	ds_mtState[cellid];

	mt_struct_stripped config = ds_MT[cellid];

	v.x = g_simSettings.maxRelVel;//DSMC_PI*DSMC_DIAM*DSMC_DIAM*sqrtf(DSMC_T/300.0)*300;
	v.y = GetRandomFloatFast(cellid, &mts, &config);

	ds_mtState[cellid] = mts;

	maximum[cellid] = v;
}

void	InitBirdData(uint	cellsCount, uint partPerCells)
{
	cudaMalloc(&maximums, cellsCount*sizeof(float2));
	cudaMalloc(&devConcentration, cellsCount*sizeof(float));

	cudaMalloc(&devCellStatistics, cellsCount*sizeof(dsmc::statistics_t));

	cudaMalloc(&devStatistics, cellsCount*sizeof(float4));

	thrust::device_ptr<float4>	dev(devStatistics);
	thrust::device_ptr<float4>	end_dev(devStatistics + cellsCount);

	thrust::fill(dev, end_dev, make_float4(0.0f, 0.0f, 0.0f, 0.0f));

	thrust::device_ptr<dsmc::statistics_t>	stat_dev(devCellStatistics);
	thrust::device_ptr<dsmc::statistics_t>	stat_dev_end(devCellStatistics + cellsCount);

	dsmc::statistics_t	stat_null;

	stat_null.numInSample = 0;
	stat_null.velMag = 0;
	stat_null.xyzSum = make_float3(0.0f, 0.0f, 0.0f);

	thrust::fill(stat_dev, stat_dev_end, stat_null);

	CreateRandomSeedsUint4(&devSeeds, g_hostSimGrid.dim.x*g_hostSimGrid.dim.y*g_hostSimGrid.dim.z);

	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.x;
	dimGrid.y = g_hostSimGrid.dim.y;
	dimGrid.z = 1;

	InitMaximums<<<dimGrid, g_hostSimGrid.dim.z>>>(maximums);
}

void	DeleteBirdData()
{
	if(counts)
		cudaFree(counts);

	if(maximums)
		cudaFree(maximums);

	if(devConcentration)
		cudaFree(devConcentration);

	if(devStatistics)
		cudaFree(devStatistics);

	if(devCellStatistics)
		cudaFree(devCellStatistics);

	if(devSeeds)
		cudaFree(devSeeds);

	if(pairs)
		cudaFree(pairs);

#ifdef	NVIDIA_RADIX_SORT
	if(radixSort)
		delete radixSort;
#endif
}

struct	DevAverage
{
	__device__	float	operator()(float v)
	{
		return v/g_simSettings.cells_count;
	}
};

//1041120 generated
//1041120

namespace	dsmc
{
	struct test_mask
	{
	    __host__ __device__
	    bool operator()(const uint4 x)
	    {
	      return x.x == 1 && x.y ==1 && x.z == 1 && x.w == 1;
	    }
	};
}

__device__	__host__	bool	operator==(const float4 a, const float4 b)
{
	return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

__host__	size_t	GenerateInitialVelocities(float4 *pos, float4* vel, float3* col, const float3 streamVector, uint partPerCell, uint numParticles, uint* geom_hashes)
{
	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.x;
	dimGrid.y = g_hostSimGrid.dim.y;
	dimGrid.z = 1;

	GenerateInitialVelocitiesDevBlocks<<<dimGrid, partPerCell>>>(pos, vel, col, streamVector, partPerCell, geom_hashes,g_hostSimGrid.dim.z);

	cudaDeviceSynchronize();
	checkCUDAError("GenerateInitialVelocitiesDevBlocks");

	return numParticles;

	thrust::device_ptr<float4> dev_pos(pos);
	thrust::device_ptr<float4> dev_pos_end(pos+numParticles);

	thrust::device_ptr<float4>	res = thrust::remove(dev_pos, dev_pos_end, make_float4(7777777,7777777,7777777,7777777));

	int rem = dev_pos_end - res;

	uint	count = numParticles - rem;

	thrust::device_ptr<float4> dev_vel(vel);
	thrust::device_ptr<float4> dev_vel_end(vel+numParticles);

	thrust::remove(dev_vel, dev_vel_end, make_float4(7777777,7777777,7777777,7777777));

	return count;

}
/*
 * test if there is a geometry in particles's cube
 * test if there is a geometry in next step particle's geometry cube
 * if line's end in other cube - track intermediate cubes
 */
#define COLLISION_LIMIT	10
#define EPS 0.0001f

#define	MODEL_WIDTH 0.0001f

//TODO: Pre compute plane coefficients for each tri
__device__	bool	CollideWithGeometry(float4* pos, float4* vel, float4* old_pos, const uint3* indices, const uint2* indicesStartEnd, const float3*	vertices, const float3* normals, uint cellid)
{
	uint2 	idx_stend = indicesStartEnd[cellid];

	float3	v1,v2,v3, p, norm;

	float3	dir = normalize(make_float3(vel->x, vel->y, vel->z));
	float	t,u,v;

	float	dist = distance4(*old_pos, *pos);

	if(dist < EPS)
		return false;

	float3	width;

	for(uint i = idx_stend.x; i < idx_stend.y; i++)
	{
		uint3	idx = indices[i];

		v1 = vertices[idx.x];
		v2 = vertices[idx.y];
		v3 = vertices[idx.z];

		if(intersect_triangle(*old_pos, dir, v1, v2, v3, &t, &u, &v))
		{
			if(t >= 0 && t <= dist)//intersection
			{
				float uv = 1-u-v;

				norm = normals[i];
				width = make_float3(norm.x*MODEL_WIDTH, norm.y*MODEL_WIDTH, norm.z*MODEL_WIDTH);

				p.x = v1.x*uv + v2.x*u + v3.x*v;
				p.y = v1.y*uv + v2.y*u + v3.y*v;
				p.z = v1.z*uv + v2.z*u + v3.z*v;

				*vel = reflect4(*vel, norm);

				float3	p1p = normalize(make_float3(vel->x, vel->y, vel->z));

				float	od = dist - t;

				pos->x = p.x + p1p.x*od;
				pos->y = p.y + p1p.y*od;
				pos->z = p.z + p1p.z*od;

				old_pos->x = p.x + width.x;
				old_pos->y = p.y + width.y;
				old_pos->z = p.z + width.z;

				return true;
			}
		}
	}

	return false;

}

struct	PeriodicZBoundaryCondition
{
	__device__	bool	operator()(float4*	pos, float4* vel, float4	old_pos) const
	{
		bool	b = false;
		//right side - respawn at left one
		if(pos->z > g_simSettings.grid_max.z)
		{
			pos->z = g_simSettings.grid_min.z+BOUNCE_DELTA + pos->z - g_simSettings.grid_max.z;
			b = true;
		}

		//top side
		if(pos->y > g_simSettings.grid_max.y)
		{
			pos->y = g_simSettings.grid_max.y - BOUNCE_DELTA;// - (pos->y-gridBoundaries_w);

			vel->y *= -1;

			b = true;
		}

		//left side - respawn at the right one
		if(pos->z < g_simSettings.grid_min.z)
		{
			pos->z = g_simSettings.grid_max.z-BOUNCE_DELTA + pos->z - g_simSettings.grid_min.z;
			b = true;
		}
		//bottom side
		if(pos->y < g_simSettings.grid_min.y)
		{
			pos->y = g_simSettings.grid_min.y + BOUNCE_DELTA;// - (pos->y-gridBoundaries_y);
			vel->y *= -1;

			b = true;
		}

		//far plane
		if(pos->x < g_simSettings.grid_min.x)
		{
			pos->x = g_simSettings.grid_min.x + BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMin);
			vel->x *= -1;

			b = true;
		}

		//near plane
		if(pos->x > g_simSettings.grid_max.x)
		{
			pos->x = g_simSettings.grid_max.x - BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMax);
			vel->x *= -1;

			b = true;
		}

		return b;
	}
};

struct	PeriodicXBoundaryCondition
{
	__device__	bool	operator()(float4*	pos, float4* vel, float4	old_pos) const
	{
		bool	b = false;
		//right side - respawn at left one
		if(pos->x > g_simSettings.grid_max.x)
		{
			pos->x = g_simSettings.grid_min.x+BOUNCE_DELTA + pos->x - g_simSettings.grid_max.x;
			b = true;
		}

		//top side
		if(pos->y > g_simSettings.grid_max.y)
		{
			pos->y = g_simSettings.grid_max.y - BOUNCE_DELTA;// - (pos->y-gridBoundaries_w);

			vel->y *= -1;

			b = true;
		}

		//left side - respawn at the right one
		if(pos->x < g_simSettings.grid_min.x)
		{
			pos->x = g_simSettings.grid_max.x-BOUNCE_DELTA + pos->x - g_simSettings.grid_min.x;
			b = true;
		}
		//bottom side
		if(pos->y < g_simSettings.grid_min.y)
		{
			pos->y = g_simSettings.grid_min.y + BOUNCE_DELTA;// - (pos->y-gridBoundaries_y);
			vel->y *= -1;

			b = true;
		}

		//far plane
		if(pos->z < g_simSettings.grid_min.z)
		{
			pos->z = g_simSettings.grid_min.z + BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMin);
			vel->z *= -1;

			b = true;
		}

		//near plane
		if(pos->z > g_simSettings.grid_max.z)
		{
			pos->z = g_simSettings.grid_max.z - BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMax);
			vel->z *= -1;

			b = true;
		}

		return b;
	}
};

struct	PeriodicYBoundaryCondition
{
	__device__	bool	operator()(float4*	pos, float4* vel, float4	old_pos) const
	{
		bool	b = false;
		//right side - respawn at left one
		if(pos->y > g_simSettings.grid_max.y)
		{
			pos->y = g_simSettings.grid_min.y+BOUNCE_DELTA + pos->y - g_simSettings.grid_max.y;
			b = true;
		}

		//top side
		if(pos->x > g_simSettings.grid_max.x)
		{
			pos->x = g_simSettings.grid_max.x - BOUNCE_DELTA;// - (pos->y-gridBoundaries_w);

			vel->x *= -1;

			b = true;
		}

		//left side - respawn at the right one
		if(pos->y < g_simSettings.grid_min.y)
		{
			pos->y = g_simSettings.grid_max.y-BOUNCE_DELTA + pos->y - g_simSettings.grid_min.y;
			b = true;
		}
		//bottom side
		if(pos->x < g_simSettings.grid_min.x)
		{
			pos->x = g_simSettings.grid_min.x + BOUNCE_DELTA;// - (pos->y-gridBoundaries_y);
			vel->x *= -1;

			b = true;
		}

		//far plane
		if(pos->z < g_simSettings.grid_min.z)
		{
			pos->z = g_simSettings.grid_min.z + BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMin);
			vel->z *= -1;

			b = true;
		}

		//near plane
		if(pos->z > g_simSettings.grid_max.z)
		{
			pos->z = g_simSettings.grid_max.z - BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMax);
			vel->z *= -1;

			b = true;
		}

		return b;
	}
};

struct	NoPeriodicBoundaryCondition
{
	__device__	bool	operator()(float4*	pos, float4* vel, float4	old_pos) const
	{
		bool	b = false;
		//right side - respawn at left one
		if(pos->y > g_simSettings.grid_max.y)
		{
			pos->y = g_simSettings.grid_max.y - BOUNCE_DELTA;// - (pos->y-gridBoundaries_w);

			vel->y *= -1;

			b = true;
		}

		//top side
		if(pos->x > g_simSettings.grid_max.x)
		{
			pos->x = g_simSettings.grid_max.x - BOUNCE_DELTA;// - (pos->y-gridBoundaries_w);

			vel->x *= -1;

			b = true;
		}

		//left side - respawn at the right one
		if(pos->y < g_simSettings.grid_min.y)
		{
			pos->y = g_simSettings.grid_min.y + BOUNCE_DELTA;// - (pos->y-gridBoundaries_y);
			vel->y *= -1;

			b = true;
		}
		//bottom side
		if(pos->x < g_simSettings.grid_min.x)
		{
			pos->x = g_simSettings.grid_min.x + BOUNCE_DELTA;// - (pos->y-gridBoundaries_y);
			vel->x *= -1;

			b = true;
		}

		//far plane
		if(pos->z < g_simSettings.grid_min.z)
		{
			pos->z = g_simSettings.grid_min.z + BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMin);
			vel->z *= -1;

			b = true;
		}

		//near plane
		if(pos->z > g_simSettings.grid_max.z)
		{
			pos->z = g_simSettings.grid_max.z - BOUNCE_DELTA;// - (pos->z-gridBoundaries_zMax);
			vel->z *= -1;

			b = true;
		}

		return b;
	}
};

template <class BoundaryConditions>
__global__	void IntegrateDev(float4* positions, float4* velocities, float3* colors, float dt, const uint numParticles, const bool useCols, BoundaryConditions BC)
{
	uint idx = blockDim.x*blockIdx.x + threadIdx.x;

	if(idx >= numParticles)
		return;

	float4	pos, old_pos;
	float4	vel;
	float3	col = make_float3(0,0,0);

	bool b = false;

	old_pos = positions[idx];
	vel = velocities[idx];

	pos.x = old_pos.x + vel.x*dt;
	pos.y = old_pos.y + vel.y*dt;
	pos.z = old_pos.z + vel.z*dt;

	b = BC(&pos, &vel, old_pos);

	positions[idx] = pos;

	if(b)
		velocities[idx] = vel;

	if(useCols)
	{
		float l = length(make_float3(vel.x, vel.y, vel.z));

		col.x = l;
		col.y = 1-l;

		colors[idx] = col;
	}
}

template <class BoundaryConditions>
__global__	void	ProcessCollisionsDev(float4 *positions, float4 *velocities, float3* colors, float dt, const uint* geom_hashes, const uint3* indices, const uint2* indicesStartEnd, const float3*	vertices, const float3* normals, const uint numParticles, const bool useCols, BoundaryConditions BC)
{
	const	uint idx = blockDim.x*blockIdx.x + threadIdx.x;

	if(idx >= numParticles)
		return;

	float4	pos, old_pos;
	float4	vel;
	float3	col = make_float3(0,0,0);

	uint st_cubid = 0xffffffff, end_cubid = 0xffffffff;

	bool	collProcessed = false;

	collProcessed = false;

	old_pos = positions[idx];
	vel = velocities[idx];

	pos.x = old_pos.x + vel.x*dt;
	pos.y = old_pos.y + vel.y*dt;
	pos.z = old_pos.z + vel.z*dt;

	uint i = 0;

	do
	{
		st_cubid  = GetCollisionCellID(old_pos);
		end_cubid = GetCollisionCellID(pos);

		if(st_cubid != 0xffffffff)
		{
			if(geom_hashes[st_cubid] != 0)
				collProcessed = CollideWithGeometry(&pos, &vel, &old_pos, indices, indicesStartEnd, vertices, normals, st_cubid);
		}

		if(!collProcessed && st_cubid != end_cubid && end_cubid != 0xffffffff)
		{
			if(geom_hashes[end_cubid] != 0)
				collProcessed = CollideWithGeometry(&pos, &vel, &old_pos, indices, indicesStartEnd, vertices, normals, end_cubid);
		}
	}while(collProcessed && i++ < COLLISION_LIMIT);

	BC(&pos, &vel, old_pos);

	positions[idx] = pos;
	velocities[idx] = vel;

	if(useCols && i > 0)
	{
		col.x = 1;
		col.y = 0;
		col.z = 0;

		if(i == COLLISION_LIMIT)
			col.y = 1;

		colors[idx] = col;
	}
}

__host__	void	Integrate(float4 *positions, float4 *velocities, float3* colors, float dt, const uint numParticles, const bool col, char bc)
{
	uint	numThreads = 128;
	uint	numBlocks = GetNumberOfBlocks(numThreads, numParticles);

	switch(bc)
	{
		case	'x':
			IntegrateDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, numParticles, col, PeriodicXBoundaryCondition());
			break;
		case	'y':
			IntegrateDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, numParticles, col, PeriodicYBoundaryCondition());
			break;
		case	'z':
			IntegrateDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, numParticles, col, PeriodicZBoundaryCondition());
			break;
		default:
			IntegrateDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, numParticles, col, NoPeriodicBoundaryCondition());
	}
}

__host__	void	IntegrateAndProcessCollisions(float4 *positions, float4 *velocities, float3* colors, float dt, const uint* geom_hashes, const uint numParticles, const bool col, char bc)
{
	uint	numThreads = 128;
	uint	numBlocks = GetNumberOfBlocks(numThreads, numParticles);

	switch(bc)
	{
		case	'x':
			ProcessCollisionsDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, geom_hashes,devIndices, devIndicesStartEnd, devVertices, devFaceNormals, numParticles, col, PeriodicXBoundaryCondition());
			break;
		case	'y':
			ProcessCollisionsDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, geom_hashes,devIndices, devIndicesStartEnd, devVertices, devFaceNormals, numParticles, col, PeriodicYBoundaryCondition());
			break;
		case	'z':
			ProcessCollisionsDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, geom_hashes,devIndices, devIndicesStartEnd, devVertices, devFaceNormals, numParticles, col, PeriodicZBoundaryCondition());
			break;
		default:
			ProcessCollisionsDev<<<numBlocks, numThreads>>>(positions, velocities, colors, dt, geom_hashes,devIndices, devIndicesStartEnd, devVertices, devFaceNormals, numParticles, col, NoPeriodicBoundaryCondition());
	}
}

__device__	unsigned int	GetCellID3D(float4 pos)
{
	uint x = (pos.x - g_simSettings.grid_min.x) / g_simSettings.cell_size.x;
	uint y = (pos.y - g_simSettings.grid_min.y) / g_simSettings.cell_size.y;
	uint z = (pos.z - g_simSettings.grid_min.z) / g_simSettings.cell_size.z;

	return min(z,g_simSettings.grid_dim.z-1)*g_simSettings.grid_dim.x*g_simSettings.grid_dim.y +
			min(y,g_simSettings.grid_dim.y-1)*g_simSettings.grid_dim.x +
				min(x,g_simSettings.grid_dim.x-1);
}

__global__	void	CalculateCellIDsDev(float4* positions, uint* d_hashes, uint* d_indices, uint numParticles)
{
	const uint idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx >= numParticles)
		return;

	float4	pos = positions[idx];

	d_hashes[idx] = GetCellID3D(pos);

	d_indices[idx] = idx;
}

__host__	void	CalculateCellIDs(float4* positions, uint* d_hashes, uint* d_indices, uint numParticles)
{
	uint	numThreads = 128;
	uint	numBlocks = GetNumberOfBlocks(numThreads, numParticles);

	CalculateCellIDsDev<<<numBlocks,numThreads>>>(positions, d_hashes, d_indices, numParticles);
}

__global__	void	FindCellBoundariesDev(uint2*   cellStartEnd,
													 float4* sorted_pos,        // output: sorted positions
													 float4* sorted_vel,
													 float3* sorted_cols,

													 uint *  gridParticleHash, // input: sorted grid hashes
													 uint *  gridParticleIndex,// input: sorted particle indices
													 float4* old_pos,            // input: sorted position array
													 float4* old_vel,
													 float3* old_cols,
													 uint numParticles, const bool col)
{
	extern __shared__ uint sharedHash[];

	uint idx = blockIdx.x*blockDim.x + threadIdx.x;

	if(idx >= numParticles)
		return;

	uint hash = gridParticleHash[idx];

	sharedHash[threadIdx.x+1] = hash;

	if (idx > 0 && threadIdx.x == 0)
	{
		sharedHash[0] = gridParticleHash[idx-1];
	}

	__syncthreads();
                                 //previous one
	if(idx == 0 || hash != sharedHash[threadIdx.x])
	{
		cellStartEnd[hash].x = idx;
		if (idx > 0)
			cellStartEnd[sharedHash[threadIdx.x]].y = idx;
	}

	if (idx == numParticles - 1)
	{
		cellStartEnd[hash].y = idx + 1;
	}

#ifdef	REORDER

	// Now use the sorted index to reorder the pos and vel data
	uint sortedIndex = gridParticleIndex[idx];

	sorted_pos[idx] = old_pos[sortedIndex];
	sorted_vel[idx] = old_vel[sortedIndex];

	if(col)
		sorted_cols[idx] = old_cols[sortedIndex];
#endif
}

__host__	void	FindCellBoundaries(uint2*   cellStartEnd,
													 float4* sorted_pos,        // output: sorted positions
													 float4* sorted_vel,
													 float3* sorted_cols,

													 uint *  gridParticleHash, // input: sorted grid hashes
													 uint *  gridParticleIndex,// input: sorted particle indices
													 float4* old_pos,            // input: sorted position array
													 float4* old_vel,
													 float3* old_cols,
													 uint numParticles, const bool col)
{
	const uint	numThreads = 256;
    const uint	smemSize   = sizeof(uint)*(numThreads+1);
    const uint	numBlocks  = GetNumberOfBlocks(numThreads, numParticles);

    FindCellBoundariesDev<<<numBlocks, numThreads, smemSize>>>(cellStartEnd, sorted_pos, sorted_vel, sorted_cols,
																gridParticleHash, gridParticleIndex, old_pos, old_vel, old_cols,
																numParticles, col);
}

#ifndef _MSC_VER

__host__	uint	CreateRandomSeeds(uint gDim, uint bDim, uint** seeds, uint& sharedMemBytesPerBlock)
{
	uint totalThreads = gDim*bDim;
	uint totalRngs    = totalThreads/WarpStandard_K;

	uint seedBytes=totalRngs*4*WarpStandard_STATE_WORDS;

	sharedMemBytesPerBlock = (bDim/WarpStandard_K)*WarpStandard_K*4;

	int ff = open("/dev/urandom", O_RDONLY);

	if(!ff)
	{
		printf("Can't open /dev/urandom\n");
		return 0;
	}

	std::vector<uint>	cpuSeeds(seedBytes/4);

	if(seedBytes != read(ff, &cpuSeeds[0], seedBytes))
	{
		printf("Can't read random seeds\n");
		return 0;
	}

	close(ff);

	cudaMalloc(seeds, seedBytes);

	cudaMemcpy(*seeds, &cpuSeeds[0], seedBytes, cudaMemcpyHostToDevice);

	return seedBytes/4;
}

__host__	uint	CreateRandomSeedsUint4(uint4** seeds, uint size)
{
	int ff = open("/dev/urandom", O_RDONLY);

	if(!ff)
	{
		printf("Can't open /dev/urandom\n");
		return 0;
	}

	std::vector<uint4>	cpuSeeds(size);

	if(size*16 != read(ff, &cpuSeeds[0], size*16))
	{
		printf("Can't read random seeds\n");
		return 0;
	}

	close(ff);

	for(uint i = 0; i < size; i++)
	{
		cpuSeeds[i].x %= (4294967296-129);
		cpuSeeds[i].x += 128;

		cpuSeeds[i].y %= (4294967296-129);
		cpuSeeds[i].y += 128;

		cpuSeeds[i].z %= (4294967296-129);
		cpuSeeds[i].z += 128;

		cpuSeeds[i].w %= (4294967296-129);
		cpuSeeds[i].w += 128;
	}

	cudaMalloc(seeds, size*sizeof(uint4));

	cudaMemcpy(*seeds, &cpuSeeds[0], size*sizeof(uint4),cudaMemcpyHostToDevice);

	return size;
}

#else

__host__	uint	CreateRandomSeedsUint4(uint4** seeds, uint size)
{
	std::vector<uint4> cpuSeeds(size);


	//Use RtlGenRandom function (SystemFunction036)
	HMODULE dll=LoadLibrary("Advapi32.dll");
	if(dll==NULL)
	{
		printf("Can't read random seeds: can't load Advapi32.dll\n");
		return 0;
	}


	typedef BOOLEAN (WINAPI *RtlFunc_ptr)(void *, ULONG);


	RtlFunc_ptr RtlFunc=(RtlFunc_ptr)GetProcAddress(dll, "SystemFunction036");

	if(RtlFunc==NULL)
	{
		printf("Can't read random seeds: can't find RtlGenRandom address\n");
		return 0;
	}


	if(!(RtlFunc)(&cpuSeeds[0], size*sizeof(uint4)))
	{
		FreeLibrary(dll);
		printf("Can't read random seeds: can't invoke RtlGenRandom function\n");
		return 0;
	}

	FreeLibrary(dll);

	cudaMalloc(seeds, size*sizeof(uint4));

	cudaMemcpy(*seeds, &cpuSeeds[0], size*sizeof(uint4),cudaMemcpyHostToDevice);

	return size;
}

#endif

__global__	void	GenerateRandomSamplesFloat4Warp(uint*	seeds, float4*	samples, uint numSamples)
{
	extern __shared__ unsigned rngShmem[];
	unsigned rngRegs[WarpStandard_REG_COUNT];

	const uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	const uint numThreads = blockDim.x*gridDim.x;

	WarpStandard_LoadState(seeds, rngRegs, rngShmem);

	float4	sample = make_float4(0.0f,0.0f,0.0f,0.0f);

	for(uint idx = tid; idx < numSamples; idx += numThreads)
	{
		sample.x = (WarpStandard_Generate(rngRegs, rngShmem) + 1.0f) / 4294967296.0f;
		sample.y = (WarpStandard_Generate(rngRegs, rngShmem) + 1.0f) / 4294967296.0f;
		sample.z = 2.0f*DSMC_PI*((WarpStandard_Generate(rngRegs, rngShmem) + 1.0f) / 4294967296.0f);
		sample.w = acosf(2.0f*((WarpStandard_Generate(rngRegs, rngShmem) + 1.0f) / 4294967296.0f)-1.0f);

		samples[idx] = sample;

		__syncthreads();
	}

	WarpStandard_SaveState(rngRegs, rngShmem, seeds);
}

__global__	void	PreProcessVectorFieldDev(float3*	field, float width, float height, float minx, float miny)
{
	uint	idx = threadIdx.x*gridDim.x + blockIdx.x;

	float3	p1 = field[idx*2+0];
	float3	p2 = field[idx*2+1];

	p1.x -= minx;
	p1.y -= miny;
	p1.z = 0;

	p1.x /= width;
	p1.y /= height;

	p2.x -= minx;
	p2.y -= miny;
	p2.z = 0;

	p2.x /= width;
	p2.y /= height;

	field[idx*2+0] = p1;
	field[idx*2+1] = p2;
}

__global__ void DumpBirdDataDev(float4* data, const float2* maximums, const uint2* cellStartEnd, const float avgConc, const float dt)
{
	const uint cellid = blockDim.x*gridDim.x*blockIdx.y + blockIdx.x*blockDim.x + threadIdx.x;
	uint2 cellStEnd = cellStartEnd[cellid];

	uint partCount = cellStEnd.y - cellStEnd.x;

	float2 vmax = maximums[cellid];

	float fpc = (0.5*partCount*avgConc*g_simSettings.phys_state.y*vmax.x*dt)/g_simSettings.cell_volume + vmax.y;
	data[cellid] = make_float4(cellStEnd.x, cellStEnd.y, partCount, fpc);
}

__host__ void DumpBirdData(const uint2* cellStartEnd, const float avgConc, const float dt, const uint cellsCount)
{

	float4* bdata = 0, *hdata = 0;

	if(bdata == 0)
		cudaMalloc(&bdata, cellsCount*sizeof(uint4));

	dim3 dimGrid;
	dimGrid.x = g_hostSimGrid.dim.x;
	dimGrid.y = g_hostSimGrid.dim.y;
	dimGrid.z = 1;

	DumpBirdDataDev<<<dimGrid, g_hostSimGrid.dim.z>>>(bdata, maximums, cellStartEnd, avgConc, dt);

	if(hdata == 0)
		hdata =	new float4[cellsCount];

	cudaMemcpy(hdata, bdata, cellsCount*sizeof(uint4), cudaMemcpyDeviceToHost);

	static uint step = 0;

	std::ofstream out("BirdData.txt", std::ios::app);
	out<<"Step: "<<step++<<"***************************************************************************************************************************\n";

	for(uint i = 0; i < cellsCount; i++)
		out<<hdata[i].x<<" "<<hdata[i].y<<" "<<hdata[i].z<<" "<<hdata[i].w<<"\n";
}

void	PreProcessVectorField(float3*	field, float width, float height, float minx, float miny)
{
	PreProcessVectorFieldDev<<<g_hostSimGrid.dim.x, g_hostSimGrid.dim.y>>>(field, width, height, minx, miny);
}

__host__	void	SortParticlesIndices(uint*	hashes, uint*	indices, uint size)
{
#ifdef	NVIDIA_RADIX_SORT

	if(radixSort == 0)
		radixSort = new nvRadixSort::RadixSort(size);

	radixSort->sort(hashes, indices, size, g_sortBits);

#else

	thrust::device_ptr<uint> keys_ptr(hashes);
	thrust::device_ptr<uint> keys_ptr_end(hashes + size);
	thrust::device_ptr<uint> values_ptr(indices);

	thrust::sort_by_key(keys_ptr, keys_ptr_end, values_ptr);
#endif
	checkCUDAErrorAndThrow("CUDA error performing RadixSort");
}

bool checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();

    if( cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
		return false;
	}

	return true;
}

void checkCUDAErrorAndThrow(const char *msg)
{
	cudaError_t err = cudaGetLastError();

	if( cudaSuccess != err)
		throw	std::runtime_error(msg);
}

