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

#include "model.h"

#include <vector>
#include <GL/glew.h>
#include <sys/stat.h>
#include <math.h>
#include <algorithm>
#include <fstream>

using namespace std;

Model::Model()
{
	
}

Model::~Model()
{
	Reset();
}

Model::surface_t::surface_t(): vertBuffer(0), indexBuffer(0), texCoordBuffer(0), normalsBuffer(0), colorsBuffer(0)
{

}

Model::surface_t::~surface_t()
{
	Reset();
}

long filelength(int f)
{
    struct stat buf;

    fstat(f, &buf);

    return(buf.st_size);
}

//TODO: Optimize
void	Model::surface_t::CalcNormals()
{
	normals.resize(vertices.size());

	for(int i = 0; i < vertices.size(); i++)
	{
		for(int j = 0; j < indices.size(); j+=3)
		{
			if(indices[j] == i)
				normals[i] += vec3_t::cross(vertices[indices[j+1]]-vertices[indices[j]],vertices[indices[j+2]]-vertices[indices[j]]);
			else if(indices[j+1] == i)
				normals[i] += vec3_t::cross(vertices[indices[j]]-vertices[indices[j+1]],vertices[indices[j+2]]-vertices[indices[j+1]]);
			else if(indices[j+2] == i)
				normals[i] += vec3_t::cross(vertices[indices[j]]-vertices[indices[j+2]],vertices[indices[j+1]]-vertices[indices[j+2]]);
		}
		normals[i].normalize();
	}
}

void	Model::surface_t::CalcFacesNormals()
{
	facesNorm.resize(indices.size()/3);
	for(int j = 0; j < indices.size(); j+=3)
	{
		facesNorm[j/3] = calc_normal(vertices[indices[j]],vertices[indices[j+1]],vertices[indices[j+2]]);
	}
}

struct	adj_indices
{
	adj_indices(uint	v): val(v) {}

	void	operator()(uint& idx) {idx += val;}

	uint	val;
};

Model::surface_t	Model::surface_t::merge(const	std::vector<Model::surface_t>& surfaces, const std::string& name)
{
	Model::surface_t	surf;

	if(!surfaces.empty())
	{
		uint	total_idx  = 0;
		uint	total_vert = 0;

		for(uint i = 0; i < surfaces.size(); i++)
		{
			total_idx  += surfaces[i].indices.size();
			total_vert += surfaces[i].vertices.size();
		}

		surf.indices.resize(total_idx);
		surf.facesNorm.resize(total_idx/3);

		surf.vertices.resize(total_vert);
		surf.texCoords.resize(total_vert);
		surf.normals.resize(total_vert);
		surf.colors.resize(total_vert);

		surf.name = name;

		uint offset = 0;
		uint vert_offset = 0;
		uint	max_idx = 0;

		for(int i = 0; i < surfaces.size(); i++)
		{
			surf.vertices.insert(surf.vertices.begin()+vert_offset, surfaces[i].vertices.begin(), surfaces[i].vertices.end());
			surf.texCoords.insert(surf.texCoords.begin()+vert_offset, surfaces[i].texCoords.begin(),surfaces[i].texCoords.end());
			surf.normals.insert(surf.normals.begin()+vert_offset, surfaces[i].normals.begin(),surfaces[i].normals.end());
			surf.colors.insert(surf.colors.begin()+vert_offset, surfaces[i].colors.begin(),surfaces[i].colors.end());

			vert_offset += surfaces[i].vertices.size();

			surf.indices.insert(surf.indices.begin()+offset, surfaces[i].indices.begin(),surfaces[i].indices.end());

			surf.facesNorm.insert(surf.facesNorm.begin()+offset/3, surfaces[i].facesNorm.begin(),surfaces[i].facesNorm.end());

			std::for_each(surf.indices.begin()+offset, surf.indices.begin()+offset+surfaces[i].indices.size(),adj_indices(max_idx));

			max_idx += surfaces[i].vertices.size();

			offset += surfaces[i].indices.size();
		}

		surf.bbox = surface_t::merge_bbox(surfaces);
	}

	return surf;
}

void	Model::surface_t::CalcBBox()
{
	bbox.min.x = bbox.min.y = bbox.min.z = 100000;
	bbox.max.x = bbox.max.y = bbox.max.z = -100000;

	for(int i = 0; i < vertices.size(); i++)
	{
		if(vertices[i].x < bbox.min.x)
			bbox.min.x = vertices[i].x;

		if(vertices[i].y < bbox.min.y)
			bbox.min.y = vertices[i].y;

		if(vertices[i].z < bbox.min.z)
			bbox.min.z = vertices[i].z;

		if(vertices[i].x > bbox.max.x)
			bbox.max.x = vertices[i].x;

		if(vertices[i].y > bbox.max.y)
			bbox.max.y = vertices[i].y;

		if(vertices[i].z > bbox.max.z)
			bbox.max.z = vertices[i].z;
	}

}

/*
 * 	unsigned int				vertBuffer;
	unsigned int				indexBuffer;
	unsigned int				texCoordBuffer;
 */

template	<class	T>
GLuint	CreateBuffer(const std::vector<T>&	data, bool elementArray, bool stat = true)
{
	GLuint	buffId = 0;

	GLenum	buffType = (elementArray) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;

	glGenBuffers(1, &buffId);
	glBindBuffer(buffType, buffId);

	GLuint	size = data.size()*sizeof(T);

	GLenum	drawType = (stat)? GL_STATIC_DRAW: GL_DYNAMIC_DRAW;

	glBufferData(buffType, size, &data[0], drawType);

	return buffId;
}

template	<class	T>
void	ModifyBuffer(uint buffId, const std::vector<T>&	data)
{
	glBindBuffer(GL_ARRAY_BUFFER, buffId);
	glBufferData(GL_ARRAY_BUFFER, data.size()*sizeof(T), &data[0], GL_DYNAMIC_DRAW);
}

void	Model::surface_t::CreateBuffers()
{
	vertBuffer 		= CreateBuffer(vertices, false);
	indexBuffer 	= CreateBuffer(indices, true);
	texCoordBuffer  = CreateBuffer(texCoords, false);
	normalsBuffer	= CreateBuffer(normals, false);
	colorsBuffer	= CreateBuffer(colors, false, false);
}

void	Model::TestAABB(const AABB_t& aabb)
{
	std::vector<uint>&	ind	= surfaces[0].indices;
	std::vector<vec3_t>&	vert = surfaces[0].vertices;
	std::vector<vec3_t>&	colors = surfaces[0].colors;

	for(int i = 0; i < colors.size(); i++)
		colors[i].reset(1,1,1);

	for(int i = 0; i < ind.size(); i+=3)
	{
		if(aabb.tri_overlaps(vert[ind[i]], vert[ind[i+1]], vert[ind[i+2]]))
		{
			colors[ind[i]].reset(1,0,0);
			colors[ind[i+1]].reset(1,0,0);
			colors[ind[i+2]].reset(1,0,0);
		}
	}

	ModifyBuffer(surfaces[0].colorsBuffer, colors);
}

enum	FileVer {FileVerMinor = 1, FileVerMajor = 0};

const	char	magic[] = "RW";

template	<class	T>
void	stream_write(ofstream& out, const T& val)
{
	out.write((char*)&val, sizeof(T));
}

template	<class	T>
void	stream_write(ofstream& out, T* v, int count)
{
	out.write((char*)v, count*sizeof(T));
}

template	<class	T>
void	stream_write(ofstream& out, const vector<T>& v)
{
	stream_write(out, v.size());
	stream_write(out, &v[0], v.size());
}

template	<class	T>
T	stream_read(ifstream& in, T& val)
{
	in.read((char*)&val, sizeof(T));
	return val;
}

template	<class	T>
T*	stream_read(ifstream& in, T* v, int count)
{
	in.read((char*)v, count*sizeof(T));
	return v;
}

template	<class	T>
void	stream_read(ifstream& in, vector<T>& v)
{
	uint32_t size = 0;

	stream_read(in, size);

	if(size)
	{
		v.resize(size);
		stream_read(in, &v[0], v.size());
	}
}

bool	Model::SaveRaw(const char* fname)
{
	ofstream	out(fname, ios_base::binary);

	if(!out.is_open())
		return false;


	stream_write(out, magic, 2);
	stream_write(out, FileVerMajor);
	stream_write(out, FileVerMinor);

	stream_write(out, bbox);
	stream_write(out, surfaces[0].vertices);
	stream_write(out, surfaces[0].texCoords);
	stream_write(out, surfaces[0].normals);
	stream_write(out, surfaces[0].facesNorm);
	stream_write(out, surfaces[0].colors);

	stream_write(out, surfaces[0].indices);

	stream_write(out, grid.gridSize);
	stream_write(out, grid.cellSize);
	stream_write(out, grid.dim);

	stream_write(out, grid.nodes.size());

	for(uint i = 0; i < grid.nodes.size(); i++)
	{
		stream_write(out, grid.nodes[i].second);
		stream_write(out, grid.nodes[i].first);
		stream_write(out, grid.normals[i]);
	}

	return true;
}
/*
 * 	bbox_t	gridSize;
	vec3_t	cellSize;
	vec3i_t	dim;
	nodes_t	nodes;

	typedef	std::vector<uint>				idx_vector_t;
	typedef	std::pair<idx_vector_t, bool>	node_t;
	typedef	std::vector<node_t>				nodes_t;

 */

bool	Model::LoadRaw(const char* fname, bool createBuffers)
{
	ifstream	in(fname, ios_base::binary);

	if(!in.is_open())
		return false;

	char	mg[2];
	stream_read(in, mg, 2);

	if(mg[0] != 'R' || mg[1] != 'W')
		return false;

	FileVer	ver1, ver2;

	stream_read(in, ver1);
	stream_read(in, ver2);

	surfaces.resize(1);

	stream_read(in, bbox);

	surfaces[0].bbox = bbox;

	stream_read(in, surfaces[0].vertices);
	stream_read(in, surfaces[0].texCoords);
	stream_read(in, surfaces[0].normals);
	stream_read(in, surfaces[0].facesNorm);
	stream_read(in, surfaces[0].colors);

	stream_read(in, surfaces[0].indices);

	stream_read(in, grid.gridSize);
	stream_read(in, grid.cellSize);
	stream_read(in, grid.dim);

	size_t	size = 0;

	stream_read(in, size);
	grid.nodes.resize(size);
	grid.normals.resize(size);

	for(uint i = 0; i < grid.nodes.size(); i++)
	{
		stream_read(in, grid.nodes[i].second);
		stream_read(in, grid.nodes[i].first);
		stream_read(in, grid.normals[i]);
	}

	if(createBuffers)
		surfaces[0].CreateBuffers();

	return true;
}

bool	Model::Load(const char* filename, bool calcNormals)
{
	int i; //Index variable

	FILE *l_file; //File pointer

	unsigned short l_chunk_id; //Chunk identifier
	unsigned int l_chunk_lenght; //Chunk lenght

	unsigned char l_char; //Char variable
	unsigned short l_qty; //Number of elements in each chunk

	unsigned short l_face_flags; //Flag that stores some face information

	vec3_t	vert;
	vec2_t	texc;

	unsigned	short	tri_data[4];

	if ((l_file=fopen (filename, "rb"))== NULL)
		return false;

	int	curr_surface = -1;

	std::string	name = "";

	while (ftell (l_file) < filelength (fileno (l_file))) //Loop to scan the whole file
	{
		fread (&l_chunk_id, 2, 1, l_file); //Read the chunk header
//		printf("ChunkID: %x\n",l_chunk_id);
		fread (&l_chunk_lenght, 4, 1, l_file); //Read the lenght of the chunk
//		printf("ChunkLenght: %x\n",l_chunk_lenght);

		switch (l_chunk_id)
        {
			//----------------- MAIN3DS -----------------
			// Description: Main chunk, contains all the other chunks
			// Chunk ID: 4d4d
			// Chunk Lenght: 0 + sub chunks
			//-------------------------------------------
			case 0x4d4d:
			break;

			//----------------- EDIT3DS -----------------
			// Description: 3D Editor chunk, objects layout info
			// Chunk ID: 3d3d (hex)
			// Chunk Lenght: 0 + sub chunks
			//-------------------------------------------
			case 0x3d3d:
			break;

			//--------------- EDIT_OBJECT ---------------
			// Description: Object block, info for each object
			// Chunk ID: 4000 (hex)
			// Chunk Lenght: len(object name) + sub chunks
			//-------------------------------------------
			case 0x4000:
				name = "";
				do
				{
					fread (&l_char, 1, 1, l_file);
                    name += l_char;
                }while(l_char != '\0');

//				printf("Surface found: %s\n", name.c_str());

				curr_surface++;
				surfaces.resize(curr_surface+1);
				surfaces[curr_surface].name = name;
			break;

			//--------------- OBJ_TRIMESH ---------------
			// Description: Triangular mesh, contains chunks for 3d mesh info
			// Chunk ID: 4100 (hex)
			// Chunk Lenght: 0 + sub chunks
			//-------------------------------------------
			case 0x4100:
			break;

			//--------------- TRI_VERTEXL ---------------
			// Description: Vertices list
			// Chunk ID: 4110 (hex)
			// Chunk Lenght: 1 x unsigned short (number of vertices)
			//             + 3 x float (vertex coordinates) x (number of vertices)
			//             + sub chunks
			//-------------------------------------------
			case 0x4110:
				fread (&l_qty, sizeof (unsigned short), 1, l_file);
				surfaces[curr_surface].vertices.resize(l_qty);
//                printf("Number of vertices: %d\n",l_qty);
                for (i=0; i<l_qty; i++)
                {
					fread (&vert.x, 3*sizeof(float), 1, l_file);
					surfaces[curr_surface].vertices[i] = vert;
				}
				break;

			//--------------- TRI_FACEL1 ----------------
			// Description: Polygons (faces) list
			// Chunk ID: 4120 (hex)
			// Chunk Lenght: 1 x unsigned short (number of polygons)
			//             + 3 x unsigned short (polygon points) x (number of polygons)
			//             + sub chunks
			//-------------------------------------------
			case 0x4120:
				fread (&l_qty, sizeof (unsigned short), 1, l_file);
				surfaces[curr_surface].indices.resize(l_qty*3);
//                printf("Number of polygons: %d\n",l_qty);
                for (i=0; i<l_qty; i++)
                {
                	fread (tri_data, 4*sizeof (unsigned short), 1, l_file);
                	surfaces[curr_surface].indices[3*i+0] = tri_data[0];
                	surfaces[curr_surface].indices[3*i+1] = tri_data[1];
                	surfaces[curr_surface].indices[3*i+2] = tri_data[2];
                }
                break;

			//------------- TRI_MAPPINGCOORS ------------
			// Description: Vertices list
			// Chunk ID: 4140 (hex)
			// Chunk Lenght: 1 x unsigned short (number of mapping points)
			//             + 2 x float (mapping coordinates) x (number of mapping points)
			//             + sub chunks
			//-------------------------------------------
			case 0x4140:
				fread (&l_qty, sizeof (unsigned short), 1, l_file);
				surfaces[curr_surface].texCoords.resize(l_qty);
				for (i=0; i<l_qty; i++)
				{
					fread (&texc.x, 2*sizeof (float), 1, l_file);
					surfaces[curr_surface].texCoords[i] = texc;
				}
                break;

			//----------- Skip unknow chunks ------------
			//We need to skip all the chunks that currently we don't use
			//We use the chunk lenght information to set the file pointer
			//to the same level next chunk
			//-------------------------------------------
			default:
				 fseek(l_file, l_chunk_lenght-6, SEEK_CUR);
        }
	}
	fclose (l_file);

	printf("%d surfaces loaded\n", surfaces.size());

	for(int i = 0; i < surfaces.size(); i++)
	{
		if(calcNormals)
		{
			surfaces[i].CalcNormals();
			surfaces[i].CalcFacesNormals();
		}

		surfaces[i].colors.resize(surfaces[i].vertices.size());

		for(int j = 0; j < surfaces[i].colors.size(); j++)
			surfaces[i].colors[j].reset(1,1,1);

		surfaces[i].CalcBBox();
	}

	surface_t	surf = surface_t::merge(surfaces, "MergedSurface");

	surfaces.clear();

	surfaces.push_back(surf);

	surfaces[0].CreateBuffers();

	printf("Model loaded. Tris: %d, vertices %d\n",surfaces[0].indices.size()/3, surfaces[0].vertices.size());

	bbox = surface_t::merge_bbox(surfaces);

	return true;
}

bbox_t	Model::surface_t::merge_bbox(const	std::vector<surface_t>& surfaces)
{
	bbox_t	bbox;

	if(!surfaces.empty())
	{
		bbox = surfaces[0].bbox;

		for(int i = 1; i < surfaces.size(); i++)
		{
			bbox.min.x = std::min(bbox.min.x,surfaces[i].bbox.min.x);
			bbox.min.y = std::min(bbox.min.y,surfaces[i].bbox.min.y);
			bbox.min.z = std::min(bbox.min.z,surfaces[i].bbox.min.z);

			bbox.max.x = std::max(bbox.max.x,surfaces[i].bbox.max.x);
			bbox.max.y = std::max(bbox.max.y,surfaces[i].bbox.max.y);
			bbox.max.z = std::max(bbox.max.z,surfaces[i].bbox.max.z);
		}
	}

	return bbox;
}

void	Model::surface_t::Render()
{
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, colorsBuffer);
	glColorPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, normalsBuffer);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, texCoordBuffer);
	glTexCoordPointer(2, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void	Model::Render()
{
	for(int i = 0; i < surfaces.size(); i++)
		surfaces[i].Render();
}

void	Model::RenderTestGrid()
{
	std::vector<vec3_t>&	vert = surfaces[0].vertices;

	vec3i_t	gridDim = grid.dim;

	vec3i_t	gridIdx;

	for(gridIdx.x = 0; gridIdx.x < gridDim.x; gridIdx.x++)
		for(gridIdx.y = 0; gridIdx.y < gridDim.y; gridIdx.y++)
			for(gridIdx.z = 0; gridIdx.z < gridDim.z; gridIdx.z++)
			{
				std::vector<uint>&		ind	= grid[grid.getCellIdx(gridIdx)].first;

				glColor3f(gridIdx.x/5.0,gridIdx.y/5.0,gridIdx.z/5.0);

				for(int i = 0; i < ind.size(); i += 3)
				{
					glBegin(GL_TRIANGLES);

						glVertex3fv(&vert[ind[i]].x);
						glVertex3fv(&vert[ind[i+1]].x);
						glVertex3fv(&vert[ind[i+2]].x);

					glEnd();
				}
			}
}

void	Model::RenderTestGrid(const vec3i_t& idx)
{
	std::vector<vec3_t>&	vert = surfaces[0].vertices;
	std::vector<uint>&		ind	= grid[grid.getCellIdx(idx)].first;

	for(int i = 0; i < ind.size(); i += 3)
	{
		glBegin(GL_TRIANGLES);

			glVertex3fv(&vert[ind[i]].x);
			glVertex3fv(&vert[ind[i+1]].x);
			glVertex3fv(&vert[ind[i+2]].x);

		glEnd();
	}
}

void	Model::RenderTestGrid(uint idx)
{
	std::vector<vec3_t>&	vert = surfaces[0].vertices;

	std::vector<uint>&		ind	= grid[idx].first;

	for(int i = 0; i < ind.size(); i += 3)
	{
		glBegin(GL_TRIANGLES);

			glVertex3fv(&vert[ind[i]].x);
			glVertex3fv(&vert[ind[i+1]].x);
			glVertex3fv(&vert[ind[i+2]].x);

		glEnd();
	}
}

void	Model::Reset()
{
	for(int i = 0; i < surfaces.size(); i++)
		surfaces[i].Reset();
}

void	Model::MoveModel(const vec3_t&	translate)
{
	for(int i = 0; i < surfaces.size(); i++)
	{
		std::vector<vec3_t>&	vert = surfaces[i].vertices;

		for(int k = 0; k < vert.size(); k++)
			vert[k] += translate;

		if(surfaces[i].vertBuffer)
		{
			glDeleteBuffers(1, &surfaces[i].vertBuffer);
			surfaces[i].vertBuffer	= CreateBuffer(vert, false);
		}

		surfaces[i].CalcBBox();
	}

	grid.gridSize.min += translate;
	grid.gridSize.max += translate;

	bbox = surface_t::merge_bbox(surfaces);
}

void	Model::CreateGrid(const vec3i_t&	dim)
{
	std::vector<uint>&		ind	= surfaces[0].indices;
	std::vector<vec3_t>&	vert = surfaces[0].vertices;

	std::vector<vec3_t>&	faceNorm = surfaces[0].facesNorm;

	vec3i_t	gridDim = dim;

	grid.reset(surfaces[0].bbox, gridDim);

	vec3i_t	gridIdx;

	for(gridIdx.x = 0; gridIdx.x < gridDim.x; gridIdx.x++)
		for(gridIdx.y = 0; gridIdx.y < gridDim.y; gridIdx.y++)
			for(gridIdx.z = 0; gridIdx.z < gridDim.z; gridIdx.z++)
			{
				AABB_t	aabb = grid.getAABB(gridIdx);
				uint	idx = grid.getCellIdx(gridIdx);

				for(uint i = 0; i < ind.size(); i += 3)
				{
					if(aabb.tri_overlaps(vert[ind[i]], vert[ind[i+1]], vert[ind[i+2]]))
					{
						grid[idx].first.push_back(ind[i]);
						grid[idx].first.push_back(ind[i+1]);
						grid[idx].first.push_back(ind[i+2]);

						grid.normals[idx].push_back(faceNorm[i/3]);

						grid[idx].second = true;
					}
				}
			}

	grid.encloseGrid();
}

void	Model::surface_t::Reset()
{
	if(vertBuffer)
		glDeleteBuffers(1, &vertBuffer);

	if(indexBuffer)
		glDeleteBuffers(1, &indexBuffer);

	if(texCoordBuffer)
		glDeleteBuffers(1, &texCoordBuffer);

	if(normalsBuffer)
		glDeleteBuffers(1, &normalsBuffer);

	if(colorsBuffer)
		glDeleteBuffers(1, &colorsBuffer);

	normals.clear();
	indices.clear();
	vertices.clear();
	texCoords.clear();
	colors.clear();
}
