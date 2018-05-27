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

#include "matrix.h"
#include "memory.h"
#include <math.h>

matrix4x4&	matrix4x4::operator*(const matrix4x4& m)
{
	matrix4x4	mm;
	return mm;
}

vec3_t&		matrix4x4::operator*(const vec3_t& v)
{
	vec3_t	vv;
	return	vv;
}

matrix4x4&	matrix4x4::inverse()
{
	matrix4x4	m;
	return m;
}

matrix4x4&	matrix4x4::transpose()
{
	matrix4x4	m;
	return m;
}

matrix4x4&	matrix4x4::invtrans()
{
	matrix4x4	m;
	return m;
}

matrix4x4	matrix4x4::inverse(const matrix4x4& m)
{
	matrix4x4	mm = m;
	return	mm.inverse();
}

matrix4x4	matrix4x4::transpose(const matrix4x4& m)
{
	matrix4x4	mm = m;
	return	mm.transpose();
}

matrix4x4	matrix4x4::invtrans(const matrix4x4& m)
{
	matrix4x4	mm = m;
	return	mm.invtrans();
}

matrix4x4	matrix4x4::identity()
{
	matrix4x4	m;
	memset(&m.m[0][0],0,16*sizeof(float));

	m.m[0][0] = m.m[1][1] = m.m[2][2] = m.m[3][3] = 1.0f;

	return m;
}

matrix4x4	matrix4x4::translation(const vec3_t& v)
{
	matrix4x4	m = identity();

	m.m[0][3] = v.x;
	m.m[1][3] = v.y;
	m.m[3][3] = v.z;

	return m;
}

matrix4x4	matrix4x4::rotation(const vec3_t& v)
{
	matrix4x4	m = identity();

	float	cx = cosf(v.x);
	float	sx = sinf(v.x);

	float	cy = cosf(v.y);
	float	sy = sinf(v.y);

	float	cz = cosf(v.z);
	float	sz = sinf(v.z);

	m.m[0][0] = cz;
	m.m[1][0] = sx*sy*cz+cx*sz;
	m.m[2][0] = -cx*sy*cz+sx*sz;

	m.m[0][1] = -cy*sz;
	m.m[1][1] = -sx*sy*sz+cx*cz;
	m.m[2][1] = cx*sy*sz+sx*cz;

	m.m[0][2] = sy;
	m.m[1][2] = -sx*cy;
	m.m[2][2] = cx*cy;

	return m;
}

matrix4x4	matrix4x4::scale(const vec3_t& v)
{
	matrix4x4	m = identity();

	m.m[0][0] = v.x;
	m.m[1][1] = v.y;
	m.m[2][2] = v.z;

	return m;
}
