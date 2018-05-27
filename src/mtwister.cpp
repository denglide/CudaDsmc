#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dci.h"

#include "mtwister.h"

//#include "MersenneTwister.h"


#define      DCMT_SEED 4172
#define  MT_RNG_PERIOD 607


typedef struct{
    unsigned int matrix_a;
    unsigned int mask_b;
    unsigned int mask_c;
    unsigned int seed;
} mt_struct_stripped;


#define   MT_RNG_COUNT 4096
#define          MT_MM 9
#define          MT_NN 19
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0xFFFFFFFEU
#define       MT_LMASK 0x1U
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18


static mt_struct MT[MT_RNG_COUNT];
static uint32_t state[MT_NN];

void	InitTwisters(const char* fname, unsigned int seed)
{
    FILE *fd = fopen(fname, "rb");
    if(!fd){
        printf("initMTRef(): failed to open %s\n", fname);
        printf("TEST FAILED\n");
    }

    for (int i = 0; i < MT_RNG_COUNT; i++)
	{
        if( !fread(MT + i, 16 /* sizeof(mt_struct) */ * sizeof(int), 1, fd) )
		{
            printf("initMTRef(): failed to load %s\n", fname);
            printf("TEST FAILED\n");
        }

		MT[i].state = state;
		sgenrand_mt(seed, &MT[i]);
    }

    fclose(fd);
}

float	GetRandomFloat(unsigned int tid)
{
	return ((float)genrand_mt(&MT[tid%MT_RNG_COUNT]) + 1.0f) / 4294967296.0f;
}

unsigned int	GetRandomInteger(unsigned int tid)
{
	return genrand_mt(&MT[tid%MT_RNG_COUNT]);
}

extern "C" void RandomRef(float *h_Random, int NPerRng, unsigned int seed)
{
    int iRng, iOut;
    for(iRng = 0; iRng < MT_RNG_COUNT; iRng++)
        for(iOut = 0; iOut < NPerRng; iOut++)
           h_Random[iRng * NPerRng + iOut] = ((float)genrand_mt(&MT[iRng]) + 1.0f) / 4294967296.0f;
}
