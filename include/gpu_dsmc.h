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

#ifndef GPU_DSMC_H_
#define GPU_DSMC_H_

#include "dsmc_base.h"
#include "particles.h"

#include "color_map.h"

namespace	dsmc
{

	template	<class	DataStrategy>
	class	GPUDSMCSolver: public	DSMCSolver
	{
	public:
		GPUDSMCSolver();
		~GPUDSMCSolver();

	protected:

		void	RunSimulationStep(float dt);
		void	RunBenchmarkSimulationStep(float dt);

		void	Render();
		void	InitData();

		statistics_t*	SampleStatistics(float4& collStat);

		uint	InitSystemState(const float3& streamVel);

		void	ReorderParticles(float4* pos, float4* vel, float3* col,
								float4*	old_pos, float4* old_vel, float3* old_col);

		void	UpdateSortedGrid(float4* pos);

		void	InitRegularGrid(reg_grid_t& grid, const float* v, uint size);

		void	RenderConcMap(uint	x, uint	y, uint	w, uint	h);

		unsigned char*	GetColorMap(uint& w, uint& h, uint&	bpp);

		float3*	GetVelocityField();

		DataStrategy	particles;

		ColorMap<gpu_data_provider_t>	colorMap;

		float*	concentration;
		uint*	gridHashes;
		uint*	cellGridIndices;
		uint2*	cellStartEnd;
		uint*	geomHashes;
	};
}

#include "gpu_dsmc.inl"

#endif /* GPU_DSMC_H_ */
