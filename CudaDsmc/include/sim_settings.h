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

#ifndef SIM_SETTINGS_H
#define SIM_SETTINGS_H

#include "config.h"
#include "gas_props.h"
#include "aabb.h"
#include <string>

struct	settings_t
{
	bool	benchmark;
	bool	cpu;
	bool	geometry;
	bool	mt;

	uint	numParticles;
	uint	partPerCell;

	float	maxRunTime;

	float	dt;

	dsmc::gas_props_t	gas;

	std::string	gas_name;

	float	density;
	float	fnum;
	float	temp;

	vec3_t	streamVelocity;

	char	periodicCondition;

	bbox_t	boundaries;
	vec3_t	cell_size;
	vec3i_t	grid_dim;

	uint	getCellsCount() const { return grid_dim.x*grid_dim.y*grid_dim.z; }

	void	readFromFile(const char* filename);
	void	readFromCmdLine(int argc, char** argv);
	void	print();

	settings_t();
};

#endif
