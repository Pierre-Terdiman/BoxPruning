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

struct SIMD_AABB_X
{
	__forceinline	SIMD_AABB_X()	{}
	__forceinline	~SIMD_AABB_X()	{}

	void	InitFrom(const AABB& b)
	{
		mMinX	= b.mMin.x;
		mMaxX	= b.mMax.x;
	}

    float mMinX;
    float mMaxX;
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

static void /*__cdecl*/ /*__forceinline*/ outputPair(udword id0, udword id1, Container& pairs, const udword* remap)
{
	pairs.Add(id0).Add(remap[id1]);
}

static void outputPair(udword id0, Container& pairs, const udword* remap, const char* const CurrentBoxListX, SIMD_AABB_X* BoxListX, udword Offset)
{
	const udword id1 = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
	pairs.Add(id0).Add(remap[id1]);
}

#define BIP_VERSION4

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
bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB* list0, udword nb1, const AABB* list1, Container& pairs)
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
		BoxListX0[nb0].mMinX = FLT_MAX;
#ifdef BIP_VERSION4
		BoxListX0[nb0+1].mMinX = FLT_MAX;
		BoxListX0[nb0+2].mMinX = FLT_MAX;
		BoxListX0[nb0+3].mMinX = FLT_MAX;
		BoxListX0[nb0+4].mMinX = FLT_MAX;
		BoxListX0[nb0+5].mMinX = FLT_MAX;
#endif

		for(udword i=0;i<nb1;i++)
		{
			const udword SortedIndex = Remap1[i];
			BoxListX1[i].InitFrom(list1[SortedIndex]);
			BoxListYZ1[i].InitFrom(list1[SortedIndex]);
		}
		BoxListX1[nb1].mMinX = FLT_MAX;
#ifdef BIP_VERSION4
		BoxListX1[nb1+1].mMinX = FLT_MAX;
		BoxListX1[nb1+2].mMinX = FLT_MAX;
		BoxListX1[nb1+3].mMinX = FLT_MAX;
		BoxListX1[nb1+4].mMinX = FLT_MAX;
		BoxListX1[nb1+5].mMinX = FLT_MAX;
#endif

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}

#define BIP_VERSION2_UNROLL

	// 3) Prune the lists
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const SIMD_AABB_X& Box0X = BoxListX0[Index0];

		const float MinLimit = Box0X.mMinX;
		while(BoxListX1[RunningAddress1].mMinX<MinLimit)
			RunningAddress1++;

		const SIMD_AABB_YZ& Box0YZ = BoxListYZ0[Index0];
		SIMD_OVERLAP_INIT(Box0YZ)

		const udword RIndex0 = Remap0[Index0];
		const float MaxLimit = Box0X.mMaxX;

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
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
		pairs.Add(RIndex0).Add(Remap1[Index1]);
		}
		_asm	align 16
StartLoop4:
		while(*(const float*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
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
#define BLOCK	if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)			\
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

	////

	Index0 = 0;
	udword RunningAddress0 = 0;
	while(RunningAddress0<nb0 && Index0<nb1)
	{
		const SIMD_AABB_X& Box1X = BoxListX1[Index0];

		const float MinLimit = Box1X.mMinX;
		while(BoxListX0[RunningAddress0].mMinX<=MinLimit)
			RunningAddress0++;

		const SIMD_AABB_YZ& Box1YZ = BoxListYZ1[Index0];
		SIMD_OVERLAP_INIT(Box1YZ)

		const udword RIndex1 = Remap1[Index0];
		const float MaxLimit = Box1X.mMaxX;

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ0[RunningAddress0];
		const char* const CurrentBoxListX = (const char*)&BoxListX0[RunningAddress0];

#ifdef BIP_VERSION4
#define BLOCK4(x, label)	{const float* box = (const float*)(CurrentBoxListYZ + Offset*2 + x*2);	\
							if(SIMD_OVERLAP_TEST)													\
								goto label;	}
		goto StartLoop4b;
		_asm	align 16
FoundOverlap3b:
		Offset += 8;
		_asm	align 16
FoundOverlap2b:
		Offset += 8;
		_asm	align 16
FoundOverlap1b:
		Offset += 8;
		_asm	align 16
FoundOverlap0b:
		Offset += 8;
		_asm	align 16
FoundOverlapb:
		{
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
		pairs.Add(RIndex1).Add(Remap0[Index1]);
		}
		_asm	align 16
StartLoop4b:
		while(*(const float*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0b)
			BLOCK4(8, FoundOverlap1b)
			BLOCK4(16, FoundOverlap2b)
			BLOCK4(24, FoundOverlap3b)
			Offset += 40;
			BLOCK4(-8, FoundOverlapb)
		}
#undef BLOCK4
#endif

#ifdef BIP_VERSION2_UNROLL
#define BLOCK	if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)			\
				{const float* box = (const float*)(CurrentBoxListYZ + Offset*2);\
					if(SIMD_OVERLAP_TEST)										\
						goto OverlapFound2;										\
					Offset += 8;

		goto LoopStart2;
		_asm	align	16
OverlapFound2:
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
		outputPair(RIndex1, Index1, pairs, Remap0);
		Offset += 8;
		_asm	align	16
LoopStart2:
		BLOCK
			BLOCK
				BLOCK
				}
			}
			goto LoopStart2;
		}
#undef BLOCK
#endif

		Index0++;
	}

	_aligned_free(BoxListYZ1);
	DELETEARRAY(BoxListX1);
	_aligned_free(BoxListYZ0);
	DELETEARRAY(BoxListX0);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Complete box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to the same set.
 *	\param		nb		[in] number of boxes
 *	\param		list	[in] list of boxes
 *	\param		pairs	[out] list of overlapping pairs
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	SIMD_AABB_X* BoxListX = new SIMD_AABB_X[nb+1+5];	// MODIFIED
	SIMD_AABB_YZ* BoxListYZ = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*(nb+1), 16);

	udword* Remap;
//	{
		// Allocate some temporary data
		float* PosList = new float[nb+1];

		// 1) Build main list using the primary axis
		for(udword i=0;i<nb;i++)
			PosList[i] = list[i].mMin.x;
		PosList[nb] = FLT_MAX;

		// 2) Sort the list
		static PRUNING_SORTER RS;	// Static for coherence
		Remap = RS.Sort(PosList, nb+1).GetRanks();

		for(udword i=0;i<nb;i++)
		{
			const udword SortedIndex = Remap[i];
			BoxListX[i].InitFrom(list[SortedIndex]);
			BoxListYZ[i].InitFrom(list[SortedIndex]);
		}
		BoxListX[nb].mMinX = FLT_MAX;
		// MODIFIED
		BoxListX[nb+1].mMinX = FLT_MAX;
		BoxListX[nb+2].mMinX = FLT_MAX;
		BoxListX[nb+3].mMinX = FLT_MAX;
		BoxListX[nb+4].mMinX = FLT_MAX;
		BoxListX[nb+5].mMinX = FLT_MAX;
		DELETEARRAY(PosList);
//	}

	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const SIMD_AABB_X& Box0X = BoxListX[Index0];

		const float MinLimit = Box0X.mMinX;
		while(BoxListX[RunningAddress++].mMinX<MinLimit);

		const SIMD_AABB_YZ& Box0YZ = BoxListYZ[Index0];
		SIMD_OVERLAP_INIT(Box0YZ)

		const float MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ[RunningAddress];
		const char* const CurrentBoxListX = (const char*)&BoxListX[RunningAddress];

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
		pairs.Add(RIndex0).Add(Remap[Index]);
//		outputPair(RIndex0, Index, pairs, Remap);
//		outputPair(RIndex0, pairs, Remap, CurrentBoxListX, BoxListX, Offset);
		}
		_asm	align 16
StartLoop4:
		while(*(const float*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0)
			BLOCK4(8, FoundOverlap1)
			BLOCK4(16, FoundOverlap2)
			BLOCK4(24, FoundOverlap3)

			Offset += 40;
			BLOCK4(-8, FoundOverlap)
		}

// The following code from Version 14a becomes our safe loop
//#define VERSION3
#ifdef VERSION3
#define BLOCK	if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)			\
				{const float* box = (const float*)(CurrentBoxListYZ + Offset*2);\
					if(SIMD_OVERLAP_TEST)										\
						goto BeforeLoop;										\
					Offset += 8;

		goto StartLoop;
		_asm	align 16
BeforeLoop:
		{const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
//		pairs.Add(RIndex0).Add(Remap[Index]);
		outputPair(RIndex0, Index, pairs, Remap);
//		outputPair(RIndex0, pairs, Remap, CurrentBoxListX, BoxListX, Offset);
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
