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

#include "reggrid.h"
#include "stdio.h"

void	reg_grid_t::reset(const	bbox_t& size, const vec3i_t&	dm)
{
	printf("Grid size min: %f %f %f\n", size.min.x, size.min.y, size.min.z);
	printf("Grid size max: %f %f %f\n", size.max.x, size.max.y, size.max.z);

	gridSize = size;
	dim		 = dm;

	vec3_t gridWidth = gridSize.max - gridSize.min;

	cellSize.x = gridWidth.x / dim.x;
	cellSize.y = gridWidth.y / dim.y;
	cellSize.z = gridWidth.z / dim.z;

	nodes.resize(cellCount());
	normals.resize(cellCount());

	for(int i = 0; i < nodes.size(); i++)
		nodes[i].second = false;

}

void	reg_grid_t::encloseGrid()
{
	vec3i_t	gridIdx;

	for(gridIdx.x = 0; gridIdx.x < dim.x; gridIdx.x++)
		for(gridIdx.y = 0; gridIdx.y < dim.y; gridIdx.y++)
		{
			vec3i_t	lastIdx = gridIdx;
			bool	lastFull = false;
			for(gridIdx.z = 0; gridIdx.z < dim.z; gridIdx.z++)
			{
				if((*this)[gridIdx].second == true)
				{
					if(lastFull == false)
					{
						lastIdx = gridIdx;
						lastFull = true;
					}else
					{
						while(lastIdx.z != gridIdx.z)
						{
							(*this)[lastIdx].second = true;
							lastIdx.z++;
						}

						lastFull = false;
					}
				}
			}
		}
}

AABB_t	reg_grid_t::getAABB(const vec3i_t&	idx)
{
	bbox_t	bbox;
	bbox.min = gridSize.min + vec3_t(idx.x*cellSize.x, idx.y*cellSize.y, idx.z*cellSize.z);
	bbox.max = bbox.min + cellSize;

	return AABB_t::from_bbox(bbox);
}

uint	reg_grid_t::getCellIdx(const vec3i_t& idx)
{
	return	dim.x*dim.y*idx.z + dim.x*idx.y + idx.x;
}
