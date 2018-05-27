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

#ifndef MATRIX_H_
#define MATRIX_H_

#include "vector.h"

struct	matrix4x4
{
	float	m[4][4];

	float&	operator()(uint i, uint j) {return m[i][j];}

	matrix4x4&	operator*(const matrix4x4& m);
	vec3_t&		operator*(const vec3_t& v);

	matrix4x4&	inverse();
	matrix4x4&	transpose();
	matrix4x4&	invtrans();

	static	matrix4x4	inverse(const matrix4x4& m);
	static	matrix4x4	transpose(const matrix4x4& m);
	static	matrix4x4	invtrans(const matrix4x4& m);
	static	matrix4x4	identity();

	static	matrix4x4	translation(const vec3_t& v);
	static	matrix4x4	rotation(const vec3_t& v);
	static	matrix4x4	scale(const vec3_t& v);
};

#endif /* MATRIX_H_ */
