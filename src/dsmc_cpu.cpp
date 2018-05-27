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

#include "dsmc_cpu.h"
#include "mtwister.h"

#include <math.h>
#include <algorithm>
#include <string.h>
#include <vector>
#include <stdio.h>

#include "dsmc.h"

using namespace std;

typedef	std::vector<unsigned int>	particle_indices_t;

particle_indices_t	g_hashTable[N_GRID_CELLS];

static const	float	cellSize_x = (gridBoundaries_z - gridBoundaries_x)/GRID_DIM_X;
static const	float	cellSize_y = (gridBoundaries_w - gridBoundaries_y)/GRID_DIM_Y;

void	GetRandomSpeed(unsigned int id, float*);
void	ComputePostCollisionVel1(unsigned int id, float*, float*);
void	ComputePostCollisionVel2(unsigned int id, float*, float*);

struct float2 
{
	float	x;
	float	y;

	float2(float _x = 0, float _y = 0): x(_x), y(_y) {}
};

float2	make_float2(float x, float y)
{
	return	float2(x,y);
}

//										test line			model line
bool	IsLinesIntersect(float2 p1, float2 p2, float2 p3, float2 p4)
{
	float2	r_min =	make_float2(min(p3.x,p4.x), min(p3.y,p4.y));
	float2	r_max =	make_float2(max(p3.x,p4.x), max(p3.y,p4.y));
	
	float ux = ((p4.x-p3.x)*(p1.y-p3.y) - (p4.y-p3.y)*(p1.x-p3.x))/
	((p4.y-p3.y)*(p2.x-p1.x) - (p4.x-p3.x)*(p2.y-p1.y));
	
	float uy = ((p2.x-p1.x)*(p1.y-p3.y)-(p2.y-p1.y)*(p1.x-p3.x))/
	((p4.y-p3.y)*(p2.x-p1.x)-(p4.x-p3.x)*(p2.y-p1.y));
	
	float2 pN = make_float2(p3.x + ux*(p4.x-p3.x),p3.y + uy*(p4.y-p3.y));
	
	if(pN.x >= r_min.x && pN.x <= r_max.x &&
	   pN.y >= r_min.y && pN.y <= r_max.y)
		return true;
	
	return false;
}

bool	IsInsideModel(float2 pos, float* model, int vert_count)
{
	float2	p1 = pos;
	float2	p2 = make_float2(pos.x+10,pos.y+10);
	
	int c = 0;
	for(int i = 0; i < vert_count-1; i++)
	{
		float2	p3 = make_float2(model[i*2 + 0], model[i*2 + 1]);
		float2	p4 = make_float2(model[(i+1)*2 + 0], model[(i+1)*2 + 1]);
		
		if(IsLinesIntersect(p1,p2,p3,p4))
			c++;
	}
	
	return c%2;
}

bool	Collide(float* vels, float* pos, float* model, int vert_count)
{
	if(IsInsideModel(make_float2(pos[0],pos[1]), model, vert_count))
	{
		vels[0] *= -1;
		vels[1] *= -1;
		pos[0] = 0;
		pos[1] = 0;
		
		return true;
	}
	
	return false;
}

unsigned int	GetCellID(float2 pos)
{
	int x = (pos.x - gridBoundaries_x) / cellSize_x;
	int y = (pos.y - gridBoundaries_y) / cellSize_y;
	
	return min(y,GRID_DIM_Y-1)*GRID_DIM_X + min(x,GRID_DIM_X-1);
}

void		ResetPositionsCPU(float* d_pos, unsigned int count)
{
	for(unsigned int i = 0; i < count*2; i++)
		d_pos[i] = 0;
}

#define BOUNCE_DELTA 0.001

void	IntegrateCPU(float*	d_vels, float*	d_pos, float* d_cols, float dt, unsigned int count)
{
	for(unsigned int i = 0; i < count; i++)
	{
		float2 vels = make_float2(d_vels[i*2 + 0],d_vels[i*2 + 1]);
		
		float old_pos = d_pos[i*2 + 0];

		d_pos[i*2 + 0] += vels.x*dt;
		d_pos[i*2 + 1] += vels.y*dt;
		
		float2 npos = make_float2(d_pos[i*2 + 0],d_pos[i*2+ 1]);
		
		if(npos.x > gridBoundaries_z)
		{
			d_pos[i*2 + 0] = gridBoundaries_x+BOUNCE_DELTA + npos.x - gridBoundaries_z;
		}
		
		if(npos.y > gridBoundaries_w)
		{
			d_pos[i*2 + 0] += npos.x - old_pos;
			d_pos[i*2 + 1] = gridBoundaries_w - BOUNCE_DELTA - (npos.y-gridBoundaries_w);
			d_vels[i*2 + 1] *= -1;
		}
		
		if(npos.x < gridBoundaries_x)
		{
			d_pos[i*2 + 0] = gridBoundaries_z+BOUNCE_DELTA + npos.x - gridBoundaries_x;
			d_vels[i*2 + 0] *= -1;						
		} 
		if(npos.y < gridBoundaries_y)
		{
			d_pos[i*2 + 0] += npos.x - old_pos;
			d_pos[i*2 + 1] = gridBoundaries_y + BOUNCE_DELTA - (npos.y-gridBoundaries_y);
			d_vels[i*2 + 1] *= -1;								
		}
		
		float	alpha = min(sqrtf(vels.x*vels.x + vels.y*vels.y),1.0f);
		
		d_cols[i*3 + 0] = alpha;								
		d_cols[i*3 + 1] = 1-alpha;								
		d_cols[i*3 + 2] = 0;
	}	
}

void	UpdateGridCPU(float *d_pos, unsigned int numParticles)
{
	for(int j = 0; j < N_GRID_CELLS; j++)
		g_hashTable[j].clear();
	
	for (int i = 0; i < numParticles; i++) 
	{
		float2 pos = make_float2(d_pos[i*2+0],d_pos[i*2+1]);
		unsigned int hash = GetCellID(pos);
		g_hashTable[hash].push_back(i);
	}
}

void	ComputeNanbuBabovskyCPU(float *vels, unsigned int bin_size)
{

}

void	ComputeNanbuCPU(float *d_speeds, unsigned int bin_size)
{
	const float p = DSMC_DT*DSMC_M/DSMC_E;
	
	for (int i = 0; i < GRID_DIM_Y; i++) 
		for (int j = 0; j < GRID_DIM_X; j++) 
		{
			unsigned int	cellid = i * GRID_DIM_X + j;
			
			unsigned int	part_in_cell = g_hashTable[cellid].size();

			for(unsigned int ii = 0; ii < part_in_cell; ii++)
			{
				unsigned int i1 = g_hashTable[cellid][ii]; 
				
				unsigned int mt_id = cellid;
				
				if(GetRandomFloat(mt_id) >= p)
				{
					unsigned int hash_j = GetRandomInteger(mt_id) % part_in_cell;

					unsigned int j = g_hashTable[cellid][hash_j];
					
					if(i1 != j)
						ComputePostCollisionVel1(mt_id,&d_speeds[i1*2],&d_speeds[j*2]);
				}
			}		
		}	
}

float	MaxwellDistr(float v)
{
	return pow(DSMC_mm/(2*DSMC_PI*DSMC_K*DSMC_T),3/2.0)*4*DSMC_PI*v*v*exp(-DSMC_mm*v*v/(2*DSMC_K*DSMC_T));
}

void	GetRandomPosition(unsigned int id, float*	pos)
{
	pos[0] = gridBoundaries_x + GetRandomFloat(id)*(gridBoundaries_z - gridBoundaries_x);
	pos[1] = gridBoundaries_y + GetRandomFloat(id)*(gridBoundaries_w - gridBoundaries_y);
}

void	GetRandomPositionCell(unsigned int id, float*	pos, float x, float y, float with, float height)
{
	pos[0] = x + GetRandomFloat(id)*with;
	pos[1] = y + GetRandomFloat(id)*height;
}

void	GenerateInitialSpeedsPerCellCPU(float* vels, float* pos, float sx, float sy, unsigned int partPerCell, unsigned int bin_size)
{
	printf("Particles per cell: %d\n", partPerCell);

	for (int i = 0; i < GRID_DIM_Y; i++)
		for (int j = 0; j < GRID_DIM_X; j++)
		{
			float	cell_x = gridBoundaries_x + j*cellSize_x;
			float	cell_y = gridBoundaries_y + i*cellSize_y;

			for(int k = 0; k < partPerCell; k++)
			{
				unsigned int idx = (i*GRID_DIM_X + j)*partPerCell + k;

				GetRandomSpeed(idx, &vels[idx*2]);

				vels[idx*2+0] += sx;
				vels[idx*2+1] += sy;

				GetRandomPositionCell(idx, &pos[idx*2], cell_x, cell_y, cellSize_x, cellSize_y);
			}
		}
}
void	GenerateInitialSpeedsCPU(float* vels, float* pos, float sx, float sy, unsigned int count)
{
	for(unsigned int k = 0; k < count; k++)
	{
		GetRandomSpeed(k,&vels[k*2]);

		vels[k*2 + 0] += sx;
		vels[k*2 + 1] += sy;

		GetRandomPosition(k,&pos[k*2]);
	}
}

float	SampleDensityFunction(unsigned int id)
{
	const float	speed_min = 0;
	const float speed_max = 2;
	
	for(;;)
	{
		float x1 = GetRandomFloat(id);
		float x2 = GetRandomFloat(id); //n2
		
		float n1 = speed_min + x1*(speed_max - speed_min);
		
		if(MaxwellDistr(n1) <= x2)
			return n1;
	}
	
}

void	GetRandomSpeed(unsigned int id, float*	v)
{
	float	vv = SampleDensityFunction(id);
	
	v[0] = vv*cos(GetRandomFloat(id)*2*DSMC_PI)*sin(GetRandomFloat(id)*DSMC_PI);
	v[1] = vv*sin(GetRandomFloat(id)*2*DSMC_PI)*sin(GetRandomFloat(id)*DSMC_PI);
}
//bug here
/*
void	ComputePostCollisionVel1(unsigned int id, float* v1, float* v2)
{
	float angle = 2*PI*GetRandomFloat(id);

	float nv1x = (v1[0]+v2[0])/2 + fabs(v1[0]-v2[0])/2*cos(angle);
	float nv1y = (v1[1]+v2[1])/2 + fabs(v1[1]-v2[1])/2*sin(angle);

	v1[0] = nv1x;
	v1[1] = nv1y;
}
*/

void	ComputePostCollisionVel1(unsigned int id, float* v1, float* v2)
{
	float angle = 2*DSMC_PI*GetRandomFloat(id);

	float mag = sqrtf((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]))/2;
	
	float nv1x = (v1[0]+v2[0])/2 + mag*cos(angle);
	float nv1y = (v1[1]+v2[1])/2 + mag*sin(angle);

	v1[0] = nv1x;
	v1[1] = nv1y;
}

void	ComputePostCollisionVel2(unsigned int id, float* v1, float* v2)
{
	float angle = 2*DSMC_PI*GetRandomFloat(id);
	
	float nv1x = (v1[0]+v2[0])/2 + fabs(v1[0]-v2[0])/2*cos(angle);
	float nv1y = (v1[1]+v2[1])/2 + fabs(v1[1]-v2[1])/2*sin(angle);
	
	float nv2x = (v1[0]+v2[0])/2 - fabs(v1[0]-v2[0])/2*cos(angle);
	float nv2y = (v1[1]+v2[1])/2 - fabs(v1[1]-v2[1])/2*sin(angle);
	
	v1[0] = nv1x;
	v1[1] = nv1y;

	v2[0] = nv2x;
	v2[1] = nv2y;
}
