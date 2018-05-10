///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for the "box pruning revisited" project.
 *	\file		IceBoxPruning.cpp
 *	\author		Pierre Terdiman
 *	\date		February 2017
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Precompiled Header
#include "Stdafx.h"
#include <xmmintrin.h>
#include <emmintrin.h>

using namespace Meshmerizer;

// InsertionSort has better coherence, RadixSort is better for one-shot queries.
#define PRUNING_SORTER	RadixSort
//#define PRUNING_SORTER	InsertionSort
#define SAFE_VERSION
#define USE_RYG_FUNCTION

//#ifdef USE_RYG_FUNCTION
// Munge the float bits to return produce an unsigned order-preserving
// ranking of floating-point numbers.
// (Old trick: http://stereopsis.com/radix.html FloatFlip, with a new
// spin to get rid of -0.0f)

// In /fp:precise, we can just calc "x + 0.0f" and get what we need.
// But fast math optimizes it away. Could use #pragma float_control,
// but that prohibits inlining of MungeFloat. So do this silly thing
// instead.
float g_global_this_always_zero = 0.0f;

static __forceinline udword MungeFloat(float f)
{
    union
    {
        float f;
        udword u;
        sdword s;
    } u;
    u.f = f + g_global_this_always_zero;  // NOT a nop! Canonicalizes -0.0f to +0.0f
    udword toggle = (u.s >> 31) | (1u << 31);
    return u.u ^ toggle;
}
//#endif

static __forceinline udword encodeFloat(udword ir)
{
#ifdef USE_RYG_FUNCTION
	return MungeFloat(*reinterpret_cast<const float*>(&ir));
#else
	if(ir & 0x80000000) //negative?
		return ~ir;//reverse sequence of negative numbers
	else
		return ir | 0x80000000; // flip sign
#endif
}

struct SIMD_AABB_X
{
	__forceinline	SIMD_AABB_X()	{}
	__forceinline	~SIMD_AABB_X()	{}

	void	InitFrom(const AABB& b)
	{
		const udword* Binary = reinterpret_cast<const udword*>(&b.mMin.x);
		mMinX	= encodeFloat(Binary[0]);
		mMaxX	= encodeFloat(Binary[3]);
	}

    udword mMinX;
    udword mMaxX;
};

struct SIMD_AABB_YZ
{
	__forceinline	SIMD_AABB_YZ()	{}
	__forceinline	~SIMD_AABB_YZ()	{}

	void	InitFrom(const AABB& b)
	{
#ifdef SAFE_VERSION
		mMinY	= -b.mMin.y;
		mMinZ	= -b.mMin.z;
		mMaxY	= b.mMax.y;
		mMaxZ	= b.mMax.z;
#else
		mMinY	= b.mMin.y;
		mMinZ	= b.mMin.z;
		mMaxY	= b.mMax.y;
		mMaxZ	= b.mMax.z;
#endif
	}

    float mMinY;
    float mMinZ;
    float mMaxY;
    float mMaxZ;
};

#ifdef SAFE_VERSION
	#define SIMD_OVERLAP_INIT(box)	\
		   __m128 b = _mm_shuffle_ps(_mm_load_ps(&box.mMinY), _mm_load_ps(&box.mMinY), 78);\
			const float Coeff = -1.0f;\
			b = _mm_mul_ps(b, _mm_load1_ps(&Coeff));

	#define SIMD_OVERLAP_TEST	_mm_movemask_ps(_mm_cmpngt_ps(b, _mm_load_ps(box)))==15
#else
	#define SIMD_OVERLAP_INIT(box)	\
		   const __m128 b = _mm_shuffle_ps(_mm_load_ps(&box.mMinY), _mm_load_ps(&box.mMinY), 78);

	#define SIMD_OVERLAP_TEST	_mm_movemask_ps(_mm_cmpnle_ps(_mm_load_ps(box), b))==12
#endif

static /*__forceinline*/ void /*__cdecl*/ outputPair(udword id0, udword id1, Pairs& pairs, const udword* remap)	// MODIFIED
{
	pairs.AddPair(id0, remap[id1]);		// MODIFIED
}

static /*__forceinline*/ void outputPair(udword id0, Pairs& pairs, const udword* remap, const char* const CurrentBoxListX, SIMD_AABB_X* BoxListX, udword Offset)	// MODIFIED
{
	const udword id1 = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
	pairs.AddPair(id0, remap[id1]);	// MODIFIED
}

#define BIP_VERSION4
#define BIP_VERSION2_UNROLL

template<int codepath>
static void bipartiteKernel(udword nb0, udword nb1,
							const SIMD_AABB_X* BoxListX0, const SIMD_AABB_X* BoxListX1,
							const SIMD_AABB_YZ* BoxListYZ0, const SIMD_AABB_YZ* BoxListYZ1,
							const udword* Remap0, const udword* Remap1, Pairs& pairs)	// MODIFIED
{
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const SIMD_AABB_X& Box0X = BoxListX0[Index0];

		const udword MinLimit = Box0X.mMinX;
		if(!codepath)
		{
			while(BoxListX1[RunningAddress1].mMinX<MinLimit)
				RunningAddress1++;
		}
		else
		{
			while(BoxListX1[RunningAddress1].mMinX<=MinLimit)
				RunningAddress1++;
		}

		const SIMD_AABB_YZ& Box0YZ = BoxListYZ0[Index0];
		SIMD_OVERLAP_INIT(Box0YZ)

		const udword RIndex0 = Remap0[Index0];
		const udword MaxLimit = Box0X.mMaxX;

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ1[RunningAddress1];
		const char* const CurrentBoxListX = (const char*)&BoxListX1[RunningAddress1];

#ifdef BIP_VERSION4
#define BLOCK4(x, label)	{const float* box = (const float*)(CurrentBoxListYZ + Offset*2 + x*2);	\
							if(SIMD_OVERLAP_TEST)													\
								goto label;	}
		goto StartLoop4;
		_asm	align 16
FoundOverlap3:
		Offset += 8;
		_asm	align 16
FoundOverlap2:
		Offset += 8;
		_asm	align 16
FoundOverlap1:
		Offset += 8;
		_asm	align 16
FoundOverlap0:
		Offset += 8;
		_asm	align 16
FoundOverlap:
		{
		const udword Index1 = (CurrentBoxListX + Offset - 1 - (const char*)BoxListX1)>>3;
		pairs.AddPair(RIndex0, Remap1[Index1]);	// MODIFIED
		}
		_asm	align 16
StartLoop4:
		while(*(const udword*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0)
			BLOCK4(8, FoundOverlap1)
			BLOCK4(16, FoundOverlap2)
			BLOCK4(24, FoundOverlap3)
			Offset += 40;
			BLOCK4(-8, FoundOverlap)
		}
#undef BLOCK4
#endif

#ifdef BIP_VERSION2_UNROLL
#define BLOCK	if(*(const udword*)(CurrentBoxListX + Offset)<=MaxLimit)			\
				{const float* box = (const float*)(CurrentBoxListYZ + Offset*2);\
					if(SIMD_OVERLAP_TEST)										\
						goto OverlapFound;										\
					Offset += 8;

		goto LoopStart;
		_asm	align	16
OverlapFound:
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
		outputPair(RIndex0, Index1, pairs, Remap1);
		Offset += 8;
		_asm	align	16
LoopStart:
		BLOCK
			BLOCK
				BLOCK
				}
			}
			goto LoopStart;
		}
#undef BLOCK
#endif

		Index0++;
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Bipartite box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to a different set.
 *	\param		nb0		[in] number of boxes in the first set
 *	\param		list0	[in] list of boxes for the first set
 *	\param		nb1		[in] number of boxes in the second set
 *	\param		list1	[in] list of boxes for the second set
 *	\param		pairs	[out] list of overlapping pairs
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB* list0, udword nb1, const AABB* list1, Container& pairs_)	// MODIFIED
{
	// Checkings
	if(!nb0 || !list0 || !nb1 || !list1)
		return false;

#ifdef BIP_VERSION4
	SIMD_AABB_X* BoxListX0 = new SIMD_AABB_X[nb0+1+5];
#else
	SIMD_AABB_X* BoxListX0 = new SIMD_AABB_X[nb0+1];
#endif
	SIMD_AABB_YZ* BoxListYZ0 = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*(nb0+1), 16);
#ifdef BIP_VERSION4
	SIMD_AABB_X* BoxListX1 = new SIMD_AABB_X[nb1+1+5];
#else
	SIMD_AABB_X* BoxListX1 = new SIMD_AABB_X[nb1+1];
#endif
	SIMD_AABB_YZ* BoxListYZ1 = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*(nb1+1), 16);

	udword* Remap0;
	udword* Remap1;
	{
		// Allocate some temporary data
		float* PosList0 = new float[nb0+1];
		float* PosList1 = new float[nb1+1];

		// 1) Build main lists using the primary axis
		for(udword i=0;i<nb0;i++)
			PosList0[i] = list0[i].mMin.x;
		PosList0[nb0] = FLT_MAX;
		for(udword i=0;i<nb1;i++)
			PosList1[i] = list1[i].mMin.x;
		PosList1[nb1] = FLT_MAX;

		// 2) Sort the lists
		static PRUNING_SORTER RS0, RS1;	// Static for coherence.
		Remap0 = RS0.Sort(PosList0, nb0+1).GetRanks();
		Remap1 = RS1.Sort(PosList1, nb1+1).GetRanks();

		for(udword i=0;i<nb0;i++)
		{
			const udword SortedIndex = Remap0[i];
			BoxListX0[i].InitFrom(list0[SortedIndex]);
			BoxListYZ0[i].InitFrom(list0[SortedIndex]);
		}

		const float SentinelValue = FLT_MAX;
		const udword SentinelIntegerValue = encodeFloat(*reinterpret_cast<const udword*>(&SentinelValue));

		BoxListX0[nb0].mMinX = SentinelIntegerValue;
#ifdef BIP_VERSION4
		BoxListX0[nb0+1].mMinX = SentinelIntegerValue;
		BoxListX0[nb0+2].mMinX = SentinelIntegerValue;
		BoxListX0[nb0+3].mMinX = SentinelIntegerValue;
		BoxListX0[nb0+4].mMinX = SentinelIntegerValue;
		BoxListX0[nb0+5].mMinX = SentinelIntegerValue;
#endif

		for(udword i=0;i<nb1;i++)
		{
			const udword SortedIndex = Remap1[i];
			BoxListX1[i].InitFrom(list1[SortedIndex]);
			BoxListYZ1[i].InitFrom(list1[SortedIndex]);
		}
		BoxListX1[nb1].mMinX = SentinelIntegerValue;
#ifdef BIP_VERSION4
		BoxListX1[nb1+1].mMinX = SentinelIntegerValue;
		BoxListX1[nb1+2].mMinX = SentinelIntegerValue;
		BoxListX1[nb1+3].mMinX = SentinelIntegerValue;
		BoxListX1[nb1+4].mMinX = SentinelIntegerValue;
		BoxListX1[nb1+5].mMinX = SentinelIntegerValue;
#endif

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}


	// 3) Prune the lists
	// NEW
	Pairs pairs(pairs_);

	bipartiteKernel<0>(nb0, nb1, BoxListX0, BoxListX1, BoxListYZ0, BoxListYZ1, Remap0, Remap1, pairs);
	bipartiteKernel<1>(nb1, nb0, BoxListX1, BoxListX0, BoxListYZ1, BoxListYZ0, Remap1, Remap0, pairs);

	_aligned_free(BoxListYZ1);
	DELETEARRAY(BoxListX1);
	_aligned_free(BoxListYZ0);
	DELETEARRAY(BoxListX0);
	return true;
}

/*static void test()
{
	float ZeroP = +0.0f;
	float ZeroN = -0.0f;
	if(ZeroP==ZeroN)
		printf("SAME ZERO\n");

	const float Values[] = { 1.0f, -1.0f, 2.0f, -2.0f, 0.0f, -0.0f };

	RadixSort RS;
	const udword* Sorted = RS.Sort(Values, 6).GetRanks();
	udword PrevVal0 = 0;
	udword PrevVal1 = 0;
	for(udword i=0;i<6;i++)
	{
		const udword SortedIndex = *Sorted++;
		printf("%f\n", Values[SortedIndex]);

		udword Val0 = MungeFloat(Values[SortedIndex]);
		udword Val1 = encodeFloat(*(udword*)&Values[SortedIndex]);
		if(i)
		{
			if(Val0<PrevVal0)
				printf("ERROR0!\n");
			if(Val1<PrevVal1)
				printf("ERROR1!\n");
		}
		PrevVal0 = Val0;
		PrevVal1 = Val1;
	}
	exit(0);
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Complete box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to the same set.
 *	\param		nb		[in] number of boxes
 *	\param		list	[in] list of boxes
 *	\param		pairs	[out] list of overlapping pairs
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs_)	// MODIFIED
{
//	test();

	// Checkings
	if(!nb || !list)
		return false;

	SIMD_AABB_X* BoxListX = new SIMD_AABB_X[nb+1+5];
	SIMD_AABB_YZ* BoxListYZ = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*(nb+1), 16);
//#define SORT_INTS
	udword* Remap;
//	{
		// Allocate some temporary data
#ifdef SORT_INTS
		udword* PosList = new udword[nb+1];
#else
		float* PosList = new float[nb+1];
#endif
		// 1) Build main list using the primary axis
#ifdef SORT_INTS
		for(udword i=0;i<nb;i++)
			PosList[i] = encodeFloat(*reinterpret_cast<const udword*>(&list[i].mMin.x));
		const float SentinelValue = FLT_MAX;
		const udword SentinelIntegerValue = encodeFloat(*reinterpret_cast<const udword*>(&SentinelValue));
		PosList[nb] = SentinelIntegerValue;
#else
		for(udword i=0;i<nb;i++)
			PosList[i] = list[i].mMin.x;
		PosList[nb] = FLT_MAX;
#endif

		// 2) Sort the list
		static PRUNING_SORTER RS;	// Static for coherence
#ifdef SORT_INTS
		Remap = RS.Sort(PosList, nb+1, false).GetRanks();
#else
		Remap = RS.Sort(PosList, nb+1).GetRanks();
#endif
		for(udword i=0;i<nb;i++)
		{
			const udword SortedIndex = Remap[i];
#ifdef SORT_INTS
			BoxListX[i].mMinX = PosList[SortedIndex];
			BoxListX[i].mMaxX = encodeFloat(*reinterpret_cast<const udword*>(&list[SortedIndex].mMax.x));
#else
			BoxListX[i].InitFrom(list[SortedIndex]);
#endif
			BoxListYZ[i].InitFrom(list[SortedIndex]);
		}

#ifndef SORT_INTS
		const float SentinelValue = FLT_MAX;
		const udword SentinelIntegerValue = encodeFloat(*reinterpret_cast<const udword*>(&SentinelValue));
#endif
		BoxListX[nb].mMinX = SentinelIntegerValue;
		BoxListX[nb+1].mMinX = SentinelIntegerValue;
		BoxListX[nb+2].mMinX = SentinelIntegerValue;
		BoxListX[nb+3].mMinX = SentinelIntegerValue;
		BoxListX[nb+4].mMinX = SentinelIntegerValue;
		BoxListX[nb+5].mMinX = SentinelIntegerValue;
		DELETEARRAY(PosList);
//	}

	// NEW
	Pairs pairs(pairs_);

	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const SIMD_AABB_X& Box0X = BoxListX[Index0];

		const udword MinLimit = Box0X.mMinX;
		while(BoxListX[RunningAddress++].mMinX<MinLimit);

		const SIMD_AABB_YZ& Box0YZ = BoxListYZ[Index0];
		SIMD_OVERLAP_INIT(Box0YZ)

		const udword MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ[RunningAddress];
		const char* const CurrentBoxListX = (const char*)&BoxListX[RunningAddress];

#define VERSION4c
#ifdef VERSION4c
#define VERSION3	// Enable this as our safe loop
#define BLOCK4(x, label)	{const float* box = (const float*)(CurrentBoxListYZ + Offset*2 + x*2);	\
							if(SIMD_OVERLAP_TEST)													\
								goto label;	}
		goto StartLoop4;
		_asm	align 16
FoundOverlap3:
		Offset += 8;
		_asm	align 16
FoundOverlap2:
		Offset += 8;
		_asm	align 16
FoundOverlap1:
		Offset += 8;
		_asm	align 16
FoundOverlap0:
		Offset += 8;
		_asm	align 16
FoundOverlap:
		{const udword Index = (CurrentBoxListX + Offset - 8 - (const char*)BoxListX)>>3;
		pairs.AddPair(RIndex0, Remap[Index]);	// MODIFIED
		}
		_asm	align 16
StartLoop4:
		while(*(const udword*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0)
			BLOCK4(8, FoundOverlap1)
			BLOCK4(16, FoundOverlap2)
			BLOCK4(24, FoundOverlap3)
			Offset += 40;
			BLOCK4(-8, FoundOverlap)
		}
#endif

// The following code from Version 14a becomes our safe loop
//#define VERSION3
#ifdef VERSION3

#define BLOCK	if(*(const udword*)(CurrentBoxListX + Offset)<=MaxLimit)		\
				{const float* box = (const float*)(CurrentBoxListYZ + Offset*2);\
					if(SIMD_OVERLAP_TEST)										\
						goto BeforeLoop;										\
					Offset += 8;

		goto StartLoop;
		_asm	align 16
BeforeLoop:
		{const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
		outputPair(RIndex0, Index, pairs, Remap);
		Offset += 8;
		}
		_asm	align 16
StartLoop:
		BLOCK
			BLOCK
				BLOCK
					BLOCK
						BLOCK
						}
					}
				}
			}
			goto StartLoop;
		}
#endif

		Index0++;
	}

	_aligned_free(BoxListYZ);
	DELETEARRAY(BoxListX);
	return true;
}
