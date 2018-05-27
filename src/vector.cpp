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

#include "vector.h"
#include "math.h"

vec3_t	vec3_t::cross(const vec3_t& a, const vec3_t& b)
{
	vec3_t	c;

	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;

	return c;
}

float	vec3_t::dot(const vec3_t& a, const vec3_t& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

vec3_t&	vec3_t::negate()
{
	x *= -1;
	y *= -1;
	z *= -1;

	return *this;
}

vec3_t	vec3_t::neg(const vec3_t& a)
{
	vec3_t	b = a;
	return b.negate();
}

vec3_t&	vec3_t::normalize()
{
	float	l = length();

	if(l != 0)
	{
		x /= l;
		y /= l;
		z /= l;
	}

	return *this;
}

vec3_t	calc_normal(const vec3_t& v1,const vec3_t& v2,const vec3_t& v3)
{
	return vec3_t::norm(vec3_t::cross(v2-v1, v3-v1));
}

vec3_t	operator+(const vec3_t& a, const vec3_t& b)
{
	vec3_t	c = a;

	return c += b;
}

vec3_t	operator-(const vec3_t& a, const vec3_t& b)
{
	vec3_t	c = a;

	return c -= b;
}

vec3_t	operator/(const vec3_t& a, float v)
{
	vec3_t t = a;
	return t /= v;
}

vec3_t	operator*(const vec3_t& a, float v)
{
	vec3_t t = a;
	return t *= v;
}

float	vec3_t::length() const
{
	return sqrtf(x*x+y*y+z*z);
}

vec3_t	vec3_t::norm(const vec3_t& a)
{
	vec3_t	b = a;
	b.normalize();
	return b;
}
