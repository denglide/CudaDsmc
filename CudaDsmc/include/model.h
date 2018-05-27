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

#ifndef	MODEL_H
#define MODEL_H

#include <vector>

#include <string>

#include "config.h"

#include "vector.h"

#include "aabb.h"

#include "reggrid.h"

/*
 * store vertices, texcoords, normals etc in separate arrays for all subsets flat
 * store 2 sets of indices:
 * 1. subsurface indices
 * 2. grid structure
 */

class Model
{
public:
	Model();
	~Model();
	
	struct	surface_t;

	bool	Load(const char*, bool calcNormals);
//TODO: use collada???
	bool	LoadRaw(const char*, bool createBuffers = true);
	bool	SaveRaw(const char*);

	const	std::vector<Model::surface_t>&	GetSurfaces()
	{
		return	surfaces;
	}

	bbox_t	GetBBox()
	{
		return bbox;
	}

	void	Render();
	void	Reset();

	void	RenderTestGrid();
	void	RenderTestGrid(uint);
	void	RenderTestGrid(const vec3i_t& idx);

	void	TestAABB(const AABB_t& aabb);

	reg_grid_t&	GetGrid() { return grid; }

	void	MoveModel(const vec3_t&	translate);

	void	CreateGrid(const vec3i_t&	dim);

	struct	surface_t
	{
		unsigned int				vertBuffer;
		unsigned int				indexBuffer;
		unsigned int				texCoordBuffer;
		unsigned int				normalsBuffer;
		unsigned int				colorsBuffer;

		bbox_t						bbox;

		std::vector<uint> 			indices;
		std::vector<vec3_t>			vertices;
		std::vector<vec3_t>			colors;
		std::vector<vec3_t>			normals;
		std::vector<vec3_t>			facesNorm;
		std::vector<vec2_t>			texCoords;

		std::string					name;

		surface_t();
		~surface_t();

		void	Render();
		void	Reset();
		void	CalcNormals();
		void	CreateBuffers();
		void	CalcBBox();

		void	CalcFacesNormals();

		static	Model::surface_t	merge(const	std::vector<surface_t>& surfaces, const std::string& name);
		static	bbox_t				merge_bbox(const	std::vector<surface_t>& surfaces);
	};

private:

	reg_grid_t						grid;
	std::vector<surface_t>			surfaces;
	bbox_t							bbox;
};

#endif

