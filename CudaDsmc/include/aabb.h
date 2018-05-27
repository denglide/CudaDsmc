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

#ifndef AABB_H_
#define AABB_H_

#include "config.h"
#include "vector.h"

struct	AABB_t;

struct	bbox_t
{
	vec3_t	min;
	vec3_t	max;

	bbox_t() {}

	bbox_t(const vec3_t& a, const vec3_t& b)
	{
		reset(a,b);
	}

	static	bbox_t	from_aabb(const	AABB_t& aabb);

	float	get_xsize() { return max.x - min.x; }
	float	get_ysize() { return max.y - min.y; }
	float	get_zsize() { return max.z - min.z; }

	void	reset(const vec3_t& a, const vec3_t& b)
	{
		min = a;
		max = b;
	}
};

struct	AABB_t
{
	vec3_t	center;
	vec3_t	half_size;

	AABB_t() {}

	AABB_t(const vec3_t& c, const vec3_t& size)
	{
		reset(c,size);
	}

	static	AABB_t	from_bbox(const	bbox_t& bbox)
	{
		AABB_t	a;

		a.half_size = (bbox.max - bbox.min)/2;
		a.center = bbox.min + a.half_size;

		return a;
	}

	void	reset(const vec3_t& c, const vec3_t& size)
	{
		center = c;
		half_size = size;
	}

	bool	tri_overlaps(const vec3_t& a, const vec3_t& b, const vec3_t& c) const;
};

#endif /* AABB_H_ */
