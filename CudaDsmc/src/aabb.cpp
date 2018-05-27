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

#include "aabb.h"
#include "math.h"
#include "tritest.h"

bbox_t	bbox_t::from_aabb(const	AABB_t& aabb)
{
	bbox_t	bbox;

	bbox.min = aabb.center - aabb.half_size;
	bbox.max = aabb.center + aabb.half_size;

	return bbox;
}

bool	AABB_t::tri_overlaps(const vec3_t& a, const vec3_t& b, const vec3_t& c) const
{
	float	bc[] = {center.x,center.y, center.z};
	float	hsz[] = {half_size.x, half_size.y, half_size.z};

	float	tris[][3] = {
			{a.x,a.y,a.z},
			{b.x,b.y,b.z},
			{c.x,c.y,c.z}
	};

	return triBoxOverlap(bc,hsz,tris) == 1;
}
