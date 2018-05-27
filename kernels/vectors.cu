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

__device__	float3	normalize(float3 v)
{
	float l = sqrtf(v.x*v.x+v.y*v.y+v.z*v.z);
	if(l > 0)
		return make_float3(v.x/l, v.y/l, v.z/l);
	else
		return v;
}

__device__ float length4(float4 vec)
{
	return	sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

__device__ float length(float3 vec)
{
	return	sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

__device__	float3	cross(const float3& a, const float3& b)
{
	float3	c;

	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;

	return c;
}

__device__	float	distance(float3 v1, float3 v2)
{
	return sqrtf((v2.x-v1.x)*(v2.x-v1.x) + (v2.y-v1.y)*(v2.y-v1.y) + (v2.z-v1.z)*(v2.z-v1.z));
}

__device__	float	distance4(float4 v1, float4 v2)
{
	return sqrtf((v2.x-v1.x)*(v2.x-v1.x) + (v2.y-v1.y)*(v2.y-v1.y) + (v2.z-v1.z)*(v2.z-v1.z));
}

__device__	float4	GetPlaneFromTri(float3 v1, float3 v2, float3 v3)
{
	float3	v12 = make_float3(v2.x - v1.x,v2.y - v1.y,v2.z - v1.z);
	float3	v13 = make_float3(v3.x - v1.x,v3.y - v1.y,v3.z - v1.z);

	float3	n = normalize(cross(v12, v13));

	return make_float4(n.x, n.y, n.z, -v1.x*n.x - v1.y*n.y - v1.z*n.z);
}

__device__	float	GetDistanceFromPlane(float4 plane, float4 p)
{
	return plane.x*p.x + plane.y*p.y + plane.z*p.z + plane.w;
}

__device__	float3	GetNormalFromTri(float3 v1, float3 v2, float3 v3)
{
	float3	v12 = make_float3(v2.x - v1.x,v2.y - v1.y,v2.z - v1.z);
	float3	v13 = make_float3(v3.x - v1.x,v3.y - v1.y,v3.z - v1.z);

	return normalize(cross(v12, v13));
}

/*
 * I' = I-2 * dot(N,I) * N;
I is the incident vector
N is the normal
dot() is the dot product
 * */

__device__	float3	reflect(float3 v, float3 n)
{
	float	dt = v.x*n.x + v.y*n.y + v.z*n.z;

	return make_float3(v.x - 2*dt*n.x, v.y - 2*dt*n.y, v.z - 2*dt*n.z);
}

__device__	float4	reflect4(float4 v, float3 n)
{
	float	dt = v.x*n.x + v.y*n.y + v.z*n.z;

	return make_float4(v.x-2*dt*n.x, v.y - 2*dt*n.y, v.z - 2*dt*n.z, v.w);

}
