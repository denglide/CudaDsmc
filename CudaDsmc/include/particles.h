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

#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <algorithm>

#include "config.h"

namespace	dsmc
{

struct	base_particles_t
{
	float4*	velocities;
	float4*	old_velocities;

	uint	numParticles;

	base_particles_t(): velocities(0), old_velocities(0), numParticles(0) {}
	~base_particles_t();

	void	init(uint samples);
	void	swap_velocities()  { std::swap(velocities, old_velocities); }
};

struct	particles_state_t: public base_particles_t
{
	float4*	positions;
	float4*	old_positions;

	float3*	colors;
	float3*	old_colors;

	particles_state_t(): positions(0), old_positions(0), colors(0), old_colors(0) {}
	~particles_state_t();

	float4*	lock_positions()       { return positions; }
	float4*	lock_old_positions()   { return old_positions; }
	void	unlock_positions()     {}
	void	unlock_old_positions() {}
	void	swap_positions()       { std::swap(positions, old_positions); }

	float3*	lock_colors()       { return colors; }
	float3*	lock_old_colors()   { return old_colors; }
	void	unlock_colors()     {}
	void	unlock_old_colors() {}
	void	swap_colors()       { std::swap(colors, old_colors); }

	void	init(uint samples);

	void	render(uint numParticles);
};

struct	particles_state_vbo_t: public base_particles_t
{
	uint	positionsVBO, oldPositionsVBO;
	uint	colorsVBO, oldColorsVBO;

	particles_state_vbo_t(): positionsVBO(0), oldPositionsVBO(0), colorsVBO(0), oldColorsVBO(0) {}
	~particles_state_vbo_t();

	float4*	lock_positions();
	float4*	lock_old_positions();
	void	unlock_positions();
	void	unlock_old_positions();
	void	swap_positions();

	float3*	lock_colors();
	float3*	lock_old_colors();
	void	unlock_colors();
	void	unlock_old_colors();
	void	swap_colors();

	void	init(uint samples);

	void	render(uint numParticles);
};

}

#endif/* PARTICLES_H_ */
