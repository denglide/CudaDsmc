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

#ifndef REGGRID_H_
#define REGGRID_H_

#include "vector.h"
#include "aabb.h"

#include <vector>

struct	reg_grid_t
{
	bbox_t	gridSize;
	vec3_t	cellSize;
	vec3i_t	dim;

	typedef	std::vector<uint>				idx_vector_t;
	typedef	std::vector<vec3_t>				normals_vector_t;
	typedef	std::vector<normals_vector_t>	normals_t;
	typedef	std::pair<idx_vector_t, bool>	node_t;
	typedef	std::vector<node_t>				nodes_t;

	nodes_t		nodes;
	normals_t	normals;

	void	reset(const	bbox_t& size, const vec3i_t&	dm);
	AABB_t	getAABB(const vec3i_t&	idx);
	uint	getCellIdx(const vec3i_t& idx);
	uint	cellCount()	{ return dim.x*dim.y*dim.z; }

	void	encloseGrid();

	node_t&	operator[](const vec3i_t& idx) { return nodes[getCellIdx(idx)]; }
	node_t&	operator[](uint i) { return nodes[i]; }
};

#endif /* REGGRID_H_ */
