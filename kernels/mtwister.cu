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

#include "MersenneTwister.h"
#include <stdio.h>
//TODO: Fix this
#include "config.h"
#include "dsmc.h"

__device__ struct mt_state
{
	int iState;
    unsigned int mti1;
    unsigned int mt[MT_NN];
};

__device__ mt_struct_stripped ds_MT[MT_RNG_COUNT];
static mt_struct_stripped h_MT[MT_RNG_COUNT];
__device__ mt_state ds_mtState[MT_RNG_COUNT];

unsigned long lcg_rand(unsigned long a)
{
    return (a * 279470273UL) % 4294967291UL;
}

void InitGPUTwisters(const char* fname, unsigned int seed)
{
	FILE *fd = fopen(fname, "rb");
    if(!fd)
        printf("initMTGPU(): failed to open %s\n", fname);

	if( !fread(h_MT, sizeof(h_MT), 1, fd) )
        printf("initMTGPU(): failed to load %s\n", fname);

	fclose(fd);

    mt_struct_stripped *MT = (mt_struct_stripped *)malloc(MT_RNG_COUNT * sizeof(mt_struct_stripped));

    for(int i = 0; i < MT_RNG_COUNT; i++)
	{
        MT[i]      = h_MT[i];
        seed = lcg_rand(seed);
        MT[i].seed = seed;
    }

	cudaMemcpyToSymbol(ds_MT, MT, sizeof(h_MT));
    free(MT);
}

__device__	unsigned int	GetRandomIntegerFast(unsigned int tid, mt_state* mts, mt_struct_stripped* config)
{
	int iState1, iStateM;
    unsigned int mti, mtiM, x;

	iState1 = mts->iState + 1;
    iStateM = mts->iState + MT_MM;

	if(iState1 >= MT_NN) 
		iState1 -= MT_NN;

    if(iStateM >= MT_NN) 
		iStateM -= MT_NN;

	mti  = mts->mti1;
    mts->mti1 = mts->mt[iState1];
    mtiM = mts->mt[iStateM];

    x    = (mti & MT_UMASK) | (mts->mti1 & MT_LMASK);
    x    =  mtiM ^ (x >> 1) ^ ((x & 1) ? config->matrix_a : 0);
    mts->mt[mts->iState] = x;
    mts->iState = iState1;

    //Tempering transformation
    x ^= (x >> MT_SHIFT0);
    x ^= (x << MT_SHIFTB) & config->mask_b;
    x ^= (x << MT_SHIFTC) & config->mask_c;
    x ^= (x >> MT_SHIFT1);

    return x;
}

__device__	float	GetRandomFloatFast(unsigned int tid, mt_state* mts, mt_struct_stripped* config)
{
    return ((float)GetRandomIntegerFast(tid,mts,config) + 1.0f) / 4294967296.0f;
}

__device__	void	InitTwisterFast(unsigned int tid, mt_state* mts, mt_struct_stripped* config)
{
	mts->mt[0] = config->seed;

	for(int i = 1; i < MT_NN; i++)
		mts->mt[i] = (1812433253U * (mts->mt[i - 1] ^ (mts->mt[i - 1] >> 30)) + i) & MT_WMASK;

    mts->iState = 0;
    mts->mti1 = mts->mt[0];
}

#define NEXT_RAND(x)         ((x) * 1664525 + 1013904223UL)

__device__ inline uint	LCGStep(unsigned int* seed)
{
	*seed = NEXT_RAND(*seed);
	return *seed;
}

__device__ inline uint	TausStep(uint* z, int S1, int S2, int S3, unsigned long M)
{
	uint	b = ((((*z) << S1) ^ (*z)) >> S2);

	*z = ((((*z) & M) << S3) ^ b);

	return  *z;
}

__device__	float	HybridTausRng(uint4*	state)
{
	return 2.3283064365387e-10 * (
			TausStep(&state->x, 13, 19, 12, 4294967294UL) ^
			TausStep(&state->y, 2, 25, 4, 4294967288UL) ^
			TausStep(&state->z, 3, 11, 17, 4294967280UL) ^
			LCGStep(&state->w)
	);
}

__device__ uint HybridTausRngInt(uint4*	state)
{
	return (
			TausStep(&state->x, 13, 19, 12, 4294967294UL) ^
			TausStep(&state->y, 2, 25, 4, 4294967288UL) ^
			TausStep(&state->z, 3, 11, 17, 4294967280UL) ^
			LCGStep(&state->w)
	);
}
/*
 * unsigned long xor128(){
static unsigned long x=123456789,y=362436069,z=521288629,w=88675123;
unsigned long t;
t=(xˆ(x<<11));x=y;y=z;z=w; return( w=(wˆ(w>>19))ˆ(tˆ(t>>8)) );
}
 *
 */
__device__	inline unsigned int	rng_xor128(uint4&	s)
{
	unsigned int t;

	t = s.x^(s.x<<11);

	s.x = s.y;
	s.y = s.z;
	s.z = s.w;

	s.w = (s.w^(s.w>>19))^(t^(t>>8));

	return s.w;
}


/*
__global__	void	GenerateRandomSamplesLCGUint4(uint4*	samples, uint numSamples)
{
	const uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	const uint numThreads = blockDim.x*gridDim.x;
	const uint bin_size = numSamples / numThreads;
	
	uint4	sample;

	uint	seed = tid;

	for(int i = 0; i < bin_size; i++)
	{
		seed = LCG_NEXT_RAND(seed);
		sample.x = seed;
		
		seed = LCG_NEXT_RAND(seed);
		sample.y = seed;

		seed = LCG_NEXT_RAND(seed);
		sample.z = seed;

		seed = LCG_NEXT_RAND(seed);
		sample.w = seed;

		samples[tid+i*numThreads] = sample;
	}
}
*/
__global__	void	GenerateRandomSamplesUint4(uint4*	samples, uint numSamples)
{
	const uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	const uint numThreads = blockDim.x*gridDim.x;
	const uint bin_size = numSamples / numThreads;

	mt_state	mts =	ds_mtState[tid];
	mt_struct_stripped config = ds_MT[tid];

	uint4	sample;

	for(int i = 0; i < bin_size; i++)
	{
		sample.x = GetRandomIntegerFast(tid,&mts, &config);
		sample.y = GetRandomIntegerFast(tid,&mts, &config);
		sample.z = GetRandomIntegerFast(tid,&mts, &config);
		sample.w = GetRandomIntegerFast(tid,&mts, &config);

		samples[tid+i*numThreads] = sample;
	}

	ds_mtState[tid] = mts;
}

__global__	void	GenerateRandomSamplesFloat2(float2*	samples, uint numSamples)
{
	const uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	const uint numThreads = blockDim.x*gridDim.x;
	const uint bin_size = numSamples / numThreads;

	mt_state	mts =	ds_mtState[tid];
	mt_struct_stripped config = ds_MT[tid];

	float2	sample;

	for(int i = 0; i < bin_size; i++)
	{
		sample.x = GetRandomFloatFast(tid,&mts, &config);
		sample.y = 2.0*DSMC_PI*GetRandomFloatFast(tid,&mts, &config);

		samples[tid+i*numThreads] = sample;
	}

	ds_mtState[tid] = mts;
}

__global__	void	GenerateRandomSamplesFloat4(float4*	samples, uint numSamples)
{
	const uint tid = blockDim.x * blockIdx.x + threadIdx.x;
	const uint numThreads = blockDim.x*gridDim.x;
	const uint bin_size = numSamples / numThreads;

	mt_state	mts =	ds_mtState[tid];
	mt_struct_stripped config = ds_MT[tid];

	float4	sample = make_float4(0.0f,0.0f,0.0f,0.0f);

//	uint	seed = tid;

	for(int i = 0; i < bin_size; i++)
	{

		sample.x = GetRandomFloatFast(tid,&mts, &config);
		sample.y = GetRandomFloatFast(tid,&mts, &config);
		sample.z = 2.0*DSMC_PI*GetRandomFloatFast(tid,&mts, &config);
/*
		seed = NEXT_RAND(seed);
		sample.x = (seed + 1.0f) / 4294967296.0f;

		seed = NEXT_RAND(seed);
		sample.y = (seed + 1.0f) / 4294967296.0f;

		seed = NEXT_RAND(seed);
		sample.z = 2.0*PI*(seed + 1.0f) / 4294967296.0f;
*/
		samples[tid+i*numThreads] = sample;
	}

	ds_mtState[tid] = mts;
}
