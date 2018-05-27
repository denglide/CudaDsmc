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

#ifndef VECTOR_H_
#define VECTOR_H_

#include "config.h"

struct	vec3_t
{
	float	x,y,z;

	vec3_t(float _x = 0, float _y = 0, float _z = 0)
	{
		reset(_x,_y,_z);
	}

	void	reset(float _x = 0, float _y = 0, float _z = 0)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	vec3_t&	operator+=(const vec3_t& c)
	{
		x += c.x;
		y += c.y;
		z += c.z;

		return *this;
	}

	vec3_t&	operator-=(const vec3_t& c)
	{
		x -= c.x;
		y -= c.y;
		z -= c.z;

		return *this;
	}

	vec3_t&	operator*=(float v)
	{
		x *= v;
		y *= v;
		z *= v;

		return *this;
	}

	vec3_t&	operator/=(float v)
	{
		x /= v;
		y /= v;
		z /= v;

		return *this;
	}

	vec3_t&	normalize();
	float	length() const;
	vec3_t&	negate();

	static	vec3_t	neg(const vec3_t& a);
	static	float	dot(const vec3_t& a, const vec3_t& b);
	static	vec3_t	norm(const vec3_t& a);
	static	vec3_t	cross(const vec3_t& a, const vec3_t& b);
};

vec3_t	calc_normal(const vec3_t& v1,const vec3_t& v2,const vec3_t& v3);

vec3_t	operator+(const vec3_t& a, const vec3_t& b);
vec3_t	operator-(const vec3_t& a, const vec3_t& b);
vec3_t	operator/(const vec3_t& a, float v);
vec3_t	operator*(const vec3_t& a, float v);

struct	vec3i_t
{
	uint	x,y,z;

	vec3i_t(uint _x = 0, uint _y = 0, uint _z = 0)
	{
		reset(_x,_y,_z);
	}

	void	reset(uint _x = 0, uint _y = 0, uint _z = 0)
	{
		x = _x;
		y = _y;
		z = _z;
	}
};

struct	vec2_t
{
	float x,y;

	vec2_t()
	{
		reset();
	}

	void	reset(float _x = 0, float _y = 0)
	{
		x = _x;
		y = _y;
	}

};

#endif /* VECTOR_H_ */
