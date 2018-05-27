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

#include "particles.h"
#include "draw_helper.h"

#include <stdexcept>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <cuda_gl_interop.h>

namespace dsmc
{

	base_particles_t::~base_particles_t()
	{
		if(velocities)
			cudaFree(velocities);

		if(old_velocities)
			cudaFree(old_velocities);
	}

	void	base_particles_t::init(uint samples)
	{
		numParticles = samples;

		cudaMalloc((void**)&velocities, samples*sizeof(float4));
		cudaMalloc((void**)&old_velocities, samples*sizeof(float4));
	}

	particles_state_t::~particles_state_t()
	{
		if(positions)
			cudaFree(positions);

		if(old_positions)
			cudaFree(old_positions);

		if(colors)
			cudaFree(colors);

		if(old_colors)
			cudaFree(old_colors);
	}

	void	particles_state_t::init(uint samples)
	{
		base_particles_t::init(samples);

		cudaMalloc((void**)&positions, samples*sizeof(float4));
		cudaMalloc((void**)&old_positions, samples*sizeof(float4));

		cudaMalloc((void**)&colors, samples*sizeof(float3));
		cudaMalloc((void**)&old_colors, samples*sizeof(float3));
	}

	void	particles_state_t::render(uint	np)
	{
		glVertexPointer(3, GL_FLOAT, 16, positions);
		glEnableClientState(GL_VERTEX_ARRAY);

		glColorPointer(3, GL_FLOAT, 0, colors);
		glEnableClientState(GL_COLOR_ARRAY);

		glDrawArrays(GL_POINTS, 0, np);

		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	}

//
//
//


	particles_state_vbo_t::~particles_state_vbo_t()
	{
		cudaGLUnregisterBufferObject(positionsVBO);
		glDeleteBuffers(1, &positionsVBO);

		cudaGLUnregisterBufferObject(colorsVBO);
		glDeleteBuffers(1, &colorsVBO);

		cudaGLUnregisterBufferObject(oldColorsVBO);
		glDeleteBuffers(1, &oldColorsVBO);

		cudaGLUnregisterBufferObject(oldPositionsVBO);
		glDeleteBuffers(1, &oldPositionsVBO);
	}

	float4*	particles_state_vbo_t::lock_positions()
	{
		float4*	positions = 0;

		cudaGLMapBufferObject((void**)&positions, positionsVBO);

		return positions;
	}

	float4*	particles_state_vbo_t::lock_old_positions()
	{
		float4*	old_positions = 0;

		cudaGLMapBufferObject((void**)&old_positions, oldPositionsVBO);

		return old_positions;
	}

	void	particles_state_vbo_t::unlock_positions()
	{
		cudaGLUnmapBufferObject(positionsVBO);
	}

	void	particles_state_vbo_t::unlock_old_positions()
	{
		cudaGLUnmapBufferObject(oldPositionsVBO);
	}

	void	particles_state_vbo_t::swap_positions()
	{
		std::swap(positionsVBO, oldPositionsVBO);
	}

	float3*	particles_state_vbo_t::lock_colors()
	{
		float3*	colors = 0;

		cudaGLMapBufferObject((void**)&colors, colorsVBO);

		return colors;
	}

	float3*	particles_state_vbo_t::lock_old_colors()
	{
		float3*	colors = 0;

		cudaGLMapBufferObject((void**)&colors, oldColorsVBO);

		return colors;
	}

	void	particles_state_vbo_t::unlock_colors()
	{
		cudaGLUnmapBufferObject(colorsVBO);
	}

	void	particles_state_vbo_t::unlock_old_colors()
	{
		cudaGLUnmapBufferObject(oldColorsVBO);
	}

	void	particles_state_vbo_t::swap_colors()
	{
		std::swap(colorsVBO, oldColorsVBO);
	}

	void	particles_state_vbo_t::init(uint samples)
	{
		base_particles_t::init(samples);

		uint size = samples*sizeof(float)*4;

		positionsVBO	= CreateCUDABufferObject(size);
		oldPositionsVBO = CreateCUDABufferObject(size);

		size = samples*sizeof(float)*3;

		colorsVBO = CreateCUDABufferObject(size);
		oldColorsVBO = CreateCUDABufferObject(size);
	}

	void	particles_state_vbo_t::render(uint np)
	{
		glBindBuffer(GL_ARRAY_BUFFER, positionsVBO);

		glVertexPointer(3, GL_FLOAT, 16, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

		glBindBuffer(GL_ARRAY_BUFFER, colorsVBO);
		glColorPointer(3, GL_FLOAT, 0, 0);
		glEnableClientState(GL_COLOR_ARRAY);

		glDrawArrays(GL_POINTS, 0, np);

		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	}
}
