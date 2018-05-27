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

#ifndef DSMC_BASE_H_
#define DSMC_BASE_H_

#include "sim_settings.h"
#include <string>

#include "draw_helper.h"

#include "model.h"

#include <ostream>

#include "statistics.h"

namespace	dsmc
{
	class	DSMCSolver
	{
	public:
		DSMCSolver();
		virtual ~DSMCSolver();

		void	Init(const settings_t& s);
		void	DoBenchmark(uint itn, float dt);

		virtual void	RunSimulationStep(float dt) = 0;
		virtual void	Render() = 0;
		virtual	uint	InitSystemState(const float3&	streamVel) = 0;

		void	RegisterSimState(simulation_state_t* state);
		void	DumpStatistics(std::ostream& out);
		void	SetModel(Model* model);

		virtual	void	RenderConcMap(uint	x, uint	y, uint	w, uint	h) = 0;

		void	RenderVectorField();

		void	SaveColorMap(const	std::string& filename);
		void	SaveVelocityField(const	std::string& filename);

	protected:

		virtual	void	RunBenchmarkSimulationStep(float dt);

		void	ProcessAndOutputStatistics(statistics_t* stat, float4 collStat, std::ostream& out);

		virtual	statistics_t*	SampleStatistics(float4& collStat) = 0;
		virtual	void	InitData() = 0;
		virtual void	InitRegularGrid(reg_grid_t& grid, const float* v, uint size) = 0;

		virtual	unsigned char*	GetColorMap(uint& w, uint& h, uint&	bpp) = 0;
		virtual	float3*	GetVelocityField() = 0;

		settings_t		m_settings;
		simulation_state_t*	m_state;
		uint	vectorFieldBuffer;
	};

	enum	SolverType
	{
		SolverGPU,
		SolverGPUVBO,
		SolverCPU
	};

	DSMCSolver*	CreateSolver(SolverType);

}

#endif /* DSMC_BASE_H_ */
