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

#include "color_map.h"

#include <cuda_gl_interop.h>

namespace	dsmc
{

void	gpu_data_provider_t::init(uint w, uint h)
{
	colorBuffer = CreateCUDABufferObject(w*h*4*sizeof(char), true);
}

void	gpu_data_provider_t::lock_data(uint width, uint height)
{
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, colorBuffer);
	glTexSubImage2D(GL_TEXTURE_2D,0,0,0,width,height,GL_RGBA,GL_UNSIGNED_BYTE, 0);
}

void	gpu_data_provider_t::unlock_data()
{
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
}

void	gpu_data_provider_t::deinit()
{
	cudaGLUnregisterBufferObject(colorBuffer);
	glDeleteBuffers(1, &colorBuffer);
}

}
