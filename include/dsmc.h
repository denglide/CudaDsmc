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

#ifndef	DSMC_H
#define DSMC_H

#define DSMC_PI 3.14159265358979f
#define DSMC_DT 2e-6f
#define	DSMC_DIAM 3.5e-10f
#define DSMC_MOL_MASS 5.0e-26f
#define DSMC_mm 5.0e-26f
#define	DSMC_K 8.617343e-5f
#define DSMC_BOLZ 1.3806e-23f
#define DSMC_T0 273.0f
#define DSMC_FNUM 0.2e14f
#define DSMC_T 300.0f

#define DSMC_DENSITY 1.0e20f

#define DSMC_M 0.6
#define DSMC_E 0.3

#define SAMPLES_COUNT (4096*256)

#define GRID_DIM_X		30
#define	GRID_DIM_Y		18
#define GRID_DIM_Z		18

#define N_GRID_CELLS	((GRID_DIM_X)*(GRID_DIM_Y)*(GRID_DIM_Z))

#define gridBoundaries_x  	 0.0f
#define gridBoundaries_y     0.0f
#define gridBoundaries_zMin  0.0f
#define gridBoundaries_z  	 0.3f
#define gridBoundaries_w  	 0.18f
#define gridBoundaries_zMax  0.18f

#define GRID_WIDTH  (gridBoundaries_z - gridBoundaries_x)
#define GRID_HEIGHT (gridBoundaries_w - gridBoundaries_y)
#define GRID_DEPTH  (gridBoundaries_zMax - gridBoundaries_zMin)

#define CELL_SIZE_X (GRID_WIDTH/GRID_DIM_X)
#define CELL_SIZE_Y (GRID_HEIGHT/GRID_DIM_Y)
#define CELL_SIZE_Z (GRID_DEPTH/GRID_DIM_Z)

#define CELL_VOLUME	(CELL_SIZE_X*CELL_SIZE_Y*CELL_SIZE_Z)
#define GRID_VOLUME (GRID_WIDTH*GRID_HEIGHT*GRID_DEPTH)

#define REORDER

#endif
