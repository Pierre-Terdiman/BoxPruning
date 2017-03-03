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

	SIMD_AABB_X* BoxListX0 = new SIMD_AABB_X[nb0+1];
	SIMD_AABB_YZ* BoxListYZ0 = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*(nb0+1), 16);
	SIMD_AABB_X* BoxListX1 = new SIMD_AABB_X[nb1+1];
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

		for(udword i=0;i<nb1;i++)
		{
			const udword SortedIndex = Remap1[i];
			BoxListX1[i].InitFrom(list1[SortedIndex]);
			BoxListYZ1[i].InitFrom(list1[SortedIndex]);
		}
		BoxListX1[nb1].mMinX = FLT_MAX;

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}

//#define BIP_VERSION1_UNROLL2
//#define BIP_VERSION2_UNROLL3
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

/*		_asm	align	16
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
				outputPair(RIndex0, Index1, pairs, Remap1);
			}
			Offset += 8;
		}*/

#ifdef BIP_VERSION1_UNROLL2
		_asm	align	16
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
				outputPair(RIndex0, Index1, pairs, Remap1);
			}
			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
				{
					const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
					outputPair(RIndex0, Index1, pairs, Remap1);
				}
				Offset += 8;
			}
		}
#endif // BIP_VERSION1_UNROLL2


#ifdef BIP_VERSION2_UNROLL3
		goto LoopStart;
		_asm	align	16
OverlapFound:
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX1)>>3;
		outputPair(RIndex0, Index1, pairs, Remap1);
		Offset += 8;
		_asm	align	16
LoopStart:
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
				goto OverlapFound;

			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
					goto OverlapFound;

				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
						goto OverlapFound;

					Offset += 8;
				}
			}
		}
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

/*		_asm	align	16
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
				outputPair(RIndex1, Index1, pairs, Remap0);
			}
			Offset += 8;
		}*/
#ifdef BIP_VERSION1_UNROLL2
		_asm	align	16
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
				outputPair(RIndex1, Index1, pairs, Remap0);
			}
			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
				{
					const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
					outputPair(RIndex1, Index1, pairs, Remap0);
				}
				Offset += 8;
			}
		}
#endif

#ifdef BIP_VERSION2_UNROLL3
		goto LoopStart2;
		_asm	align	16
OverlapFound2:
		const udword Index1 = (CurrentBoxListX + Offset - (const char*)BoxListX0)>>3;
		outputPair(RIndex1, Index1, pairs, Remap0);
		Offset += 8;
		_asm	align	16
LoopStart2:
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
				goto OverlapFound2;

			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
					goto OverlapFound2;

				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
						goto OverlapFound2;

					Offset += 8;
				}
			}
		}
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

	SIMD_AABB_X* BoxListX = new SIMD_AABB_X[nb+1];
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

//#define VERSION1_UNROLL2
#ifdef VERSION1_UNROLL2
		_asm	align	16

		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
				outputPair(RIndex0, Index, pairs, Remap);
			}
			Offset += 8;

			// START OF ADDED BLOCK
			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
				{
					const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
					outputPair(RIndex0, Index, pairs, Remap);
				}
				Offset += 8;
			}
			// END OF ADDED BLOCK
		}
#endif	// VERSION1_UNROLL2

//#define VERSION1_UNROLL3
#ifdef VERSION1_UNROLL3
		_asm	align	16

		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
			{
				const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
				outputPair(RIndex0, Index, pairs, Remap);
			}
			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
				{
					const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
					outputPair(RIndex0, Index, pairs, Remap);
				}
				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
					{
						const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
						outputPair(RIndex0, Index, pairs, Remap);
					}
					Offset += 8;
				}
			}
		}
#endif	// VERSION1_UNROLL3

//#define VERSION2_UNROLL2
#ifdef VERSION2_UNROLL2

			goto LoopStart;
	OverlapFound:
			const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
			outputPair(RIndex0, Index, pairs, Remap);
			Offset += 8;

			_asm	align	16
	LoopStart:
			while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
					goto OverlapFound;

				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
						goto OverlapFound;

					Offset += 8;
				}
			}

#endif	// VERSION2_UNROLL2

//#define VERSION2_UNROLL3
#ifdef VERSION2_UNROLL3
		goto LoopStart;
OverlapFound:
		const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
		outputPair(RIndex0, Index, pairs, Remap);
		Offset += 8;

		_asm	align	16
LoopStart:
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
				goto OverlapFound;

			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
					goto OverlapFound;

				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
						goto OverlapFound;

					Offset += 8;
				}
			}
		}
#endif	// VERSION2_UNROLL3

//#define VERSION2_UNROLL4
#ifdef VERSION2_UNROLL4
		goto LoopStart;
OverlapFound:
		const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
		outputPair(RIndex0, Index, pairs, Remap);
		Offset += 8;

		_asm	align	16
LoopStart:
		while(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
		{
			const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
			if(SIMD_OVERLAP_TEST)
				goto OverlapFound;

			Offset += 8;

			if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
			{
				const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
				if(SIMD_OVERLAP_TEST)
					goto OverlapFound;

				Offset += 8;

				if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
				{
					const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
					if(SIMD_OVERLAP_TEST)
						goto OverlapFound;

					Offset += 8;

					if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)
					{
						const float* box = (const float*)(CurrentBoxListYZ + Offset*2);
						if(SIMD_OVERLAP_TEST)
							goto OverlapFound;

						Offset += 8;
					}
				}
			}
		}
#endif	// VERSION2_UNROLL4

//#define VERSION2_UNROLL
#ifdef VERSION2_UNROLL
#define BLOCK	if(*(const float*)(CurrentBoxListX + Offset)<=MaxLimit)			\
				{const float* box = (const float*)(CurrentBoxListYZ + Offset*2);\
					if(SIMD_OVERLAP_TEST)										\
						goto OverlapFound;										\
					Offset += 8;
		goto LoopStart;
		_asm	align 16
OverlapFound:
		const udword Index = (CurrentBoxListX + Offset - (const char*)BoxListX)>>3;
		outputPair(RIndex0, Index, pairs, Remap);
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
#endif	// VERSION2_UNROLL

#define VERSION3
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
