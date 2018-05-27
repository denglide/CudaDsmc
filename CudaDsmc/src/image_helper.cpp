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

//#include <png.h>
#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>

#include "config.h"

void	save_png(const char* filename, unsigned char* data, uint w, uint h, uint bpp)
{
#ifndef _MSC_VER
	/* create file */
	FILE *fp = fopen(filename, "wb");

	if (!fp)
		throw	std::runtime_error("DSMCSolver::SaveColorMap. Can't open file");

	png_structp	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)
		throw	std::runtime_error("DSMCSolver::SaveColorMap. png_create_write_struct failed");

	png_infop	info_ptr = png_create_info_struct(png_ptr);

	if (!info_ptr)
		throw	std::runtime_error("DSMCSolver::SaveColorMap. png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		throw	std::runtime_error("DSMCSolver::SaveColorMap. Error during init_io");

	png_init_io(png_ptr, fp);

	if (setjmp(png_jmpbuf(png_ptr)))
		throw	std::runtime_error("DSMCSolver::SaveColorMap. Error during writing header");

	int	color_type;

	switch(bpp)
	{
	case 1:
		color_type = PNG_COLOR_TYPE_GRAY;
		break;
	case 4:
		color_type = PNG_COLOR_TYPE_RGB_ALPHA;
		break;
	default:
		break;
	}

	png_set_IHDR(png_ptr, info_ptr, w, h,
		     8, color_type, PNG_INTERLACE_NONE,
		     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);


	/* write bytes */
	if (setjmp(png_jmpbuf(png_ptr)))
		throw	std::runtime_error("DSMCSolver::SaveColorMap. Error during writing bytes");

	png_bytep*	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * h);

	for (uint i=0; i<h; i++)
			row_pointers[i] = (png_bytep)&data[i*w*bpp];

	png_write_image(png_ptr, row_pointers);


	/* end write */
	if (setjmp(png_jmpbuf(png_ptr)))
		throw	std::runtime_error("DSMCSolver::SaveColorMap. Error during end of write");

	png_write_end(png_ptr, NULL);

	free(row_pointers);

	fclose(fp);
#endif
}

inline	void	putPixel(unsigned char*	data, uint width, uint x, uint y, char col)
{
	data[y*width+x] = col;
}

void	draw_line(unsigned char* data, uint width, uint x0, uint y0, uint x1, uint y1, unsigned char color)
{
    int dy = y1 - y0;
    int dx = x1 - x0;
    float t = (float) 0.5;                      // offset for rounding

    putPixel(data, width, x0, y0, color);

    if (abs(dx) > abs(dy))
    {          // slope < 1
        float m = (float) dy / (float) dx;      // compute slope
        t += y0;
        dx = (dx < 0) ? -1 : 1;
        m *= dx;
        while (x0 != x1)
        {
            x0 += dx;                           // step to next x value
            t += m;                             // add slope to y value

            putPixel(data, width, x0, (uint)t, color);
        }
    } else
    {                                    // slope >= 1
        float m = (float) dx / (float) dy;      // compute slope
        t += x0;
        dy = (dy < 0) ? -1 : 1;
        m *= dy;

        while (y0 != y1)
        {
            y0 += dy;                           // step to next y value
            t += m;                             // add slope to x value

            putPixel(data, width, (uint) t, y0, color);
        }
    }
}
