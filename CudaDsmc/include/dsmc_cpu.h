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

#ifndef	DSMC_CPU_H
#define DSMC_CPU_H

void	InitGridCPU(unsigned int numCells, unsigned int numParticles);
void	ComputeNanbuCPU(float *d_speeds, unsigned int bin_size);
void	GenerateInitialSpeedsCPU(float* vels, float* pos, float sx, float sy, unsigned int count);
void	GenerateInitialSpeedsPerCellCPU(float*, float*, float sx, float sy, unsigned int partPerCell, unsigned int bin_size);
void	ComputeNanbuBabovskyCPU(float *vels, unsigned int bin_size);
void	ComputeBirdCPU(float *vels, float* d_pos, float* d_model, int vert_count, unsigned int bin_size, unsigned int count);

void	ResetPositionsCPU(float* d_pos, unsigned int count);
void	ProcessModelCPU(float* d_model, int vert_count);
void	IntegrateCPU(float*	d_vels, float*	d_pos, float* d_cols, float dt, unsigned int count);
void	UpdateGridCPU(float *d_pos, unsigned int count);


#endif
