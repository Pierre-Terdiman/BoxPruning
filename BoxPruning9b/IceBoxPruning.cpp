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

using namespace Meshmerizer;

// InsertionSort has better coherence, RadixSort is better for one-shot queries.
#define PRUNING_SORTER	RadixSort
//#define PRUNING_SORTER	InsertionSort

//#define NAIVE_VERSION
#ifndef NAIVE_VERSION
	#define SAFE_VERSION

	// NEW
	//__declspec(align(8))
	struct SIMD_AABB
	{
		__forceinline	SIMD_AABB()		{}
		__forceinline	~SIMD_AABB()	{}

		void	InitFrom(const AABB& b)
		{
			mMinX	= b.mMin.x;
			mMaxX	= b.mMax.x;
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

		float mMinX;
		float mMaxX;
		float mMinY;
		float mMinZ;
		float mMaxY;
		float mMaxZ;
	};

	// The numbers are different because the SIMD code implements intersects2D2
	// instead of the initial intersects2D.

	// NEW
	static __forceinline int intersects2D(const SIMD_AABB& a, const SIMD_AABB& b)
	{
		   if(    b.mMaxY < a.mMinY || a.mMaxY < b.mMinY
		   ||     b.mMaxZ < a.mMinZ || a.mMaxZ < b.mMinZ)
				  return 0;
		   return 1;
	}

	// NEW
	static __forceinline int intersects2D2(const SIMD_AABB& a, const SIMD_AABB& b)
	{
		   if(    b.mMaxY <= a.mMinY || a.mMaxY < b.mMinY
		   ||     b.mMaxZ <= a.mMinZ || a.mMaxZ < b.mMinZ)
				  return 0;
		   return 1;
	}

		#include <xmmintrin.h>
		#include <emmintrin.h>

	/*#define SIMD_OVERLAP_INIT(box)    \
		   const __m128i b = _mm_shuffle_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(&box.mMinY)), 78);

	#define SIMD_OVERLAP_TEST(box)    \
		   const __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&box.mMinY));  \
		   const __m128i d = _mm_cmpgt_epi32(a, b);                                                                      \
		   if(_mm_movemask_epi8(d)==0x0000ff00)*/

	// NEW
	#ifdef SAFE_VERSION
		#define SIMD_OVERLAP_INIT(box)	\
			   __m128 b = _mm_shuffle_ps(_mm_loadu_ps(&box.mMinY), _mm_loadu_ps(&box.mMinY), 78);\
				const float Coeff = -1.0f;\
				b = _mm_mul_ps(b, _mm_load1_ps(&Coeff));

		#define SIMD_OVERLAP_TEST(box)						\
			   const __m128 a = _mm_loadu_ps(&box.mMinY);	\
			   const __m128 d = _mm_cmpge_ps(a, b);			\
			   if(_mm_movemask_ps(d)==15)
	#else
		#define SIMD_OVERLAP_INIT(box)	\
			   const __m128 b = _mm_shuffle_ps(_mm_loadu_ps(&box.mMinY), _mm_loadu_ps(&box.mMinY), 78);

		#define SIMD_OVERLAP_TEST(box)						\
			   const __m128 a = _mm_loadu_ps(&box.mMinY);	\
			   const __m128 d = _mm_cmpge_ps(a, b);			\
			   if(_mm_movemask_ps(d)==12)
	#endif
#endif

/*static __forceinline int intersects2D(const AABB& a, const AABB& b)
{
       if(    b.mMax.y < a.mMin.y || a.mMax.y < b.mMin.y
       ||     b.mMax.z < a.mMin.z || a.mMax.z < b.mMin.z)
              return 0;
       return 1;
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Bipartite box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to a different set.
 *	\param		nb0		[in] number of boxes in the first set
 *	\param		list0	[in] list of boxes for the first set
 *	\param		nb1		[in] number of boxes in the second set
 *	\param		list1	[in] list of boxes for the second set
 *	\param		pairs	[out] list of overlapping pairs
 *	\param		axes	[in] projection order (0,2,1 is often best)
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB* list0, udword nb1, const AABB* list1, Container& pairs)
{
	// Checkings
	if(!nb0 || !list0 || !nb1 || !list1)
		return false;

#ifdef NAIVE_VERSION
	AABB* BoxList0 = new AABB[nb0+1];
	AABB* BoxList1 = new AABB[nb1+1];
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
			BoxList0[i] = list0[Remap0[i]];
		BoxList0[nb0].mMin.x = FLT_MAX;

		for(udword i=0;i<nb1;i++)
			BoxList1[i] = list1[Remap1[i]];
		BoxList1[nb1].mMin.x = FLT_MAX;

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}

	// 3) Prune the lists
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const AABB& Box0 = BoxList0[Index0];

		const float MinLimit = Box0.mMin.x;
		while(BoxList1[RunningAddress1].mMin.x<MinLimit)
			RunningAddress1++;

		const udword RIndex0 = Remap0[Index0];
		udword Index1 = RunningAddress1;
		const float MaxLimit = Box0.mMax.x;
		while(BoxList1[Index1].mMin.x<=MaxLimit)
		{
			const AABB& a = Box0;
			const AABB& b = BoxList1[Index1];
			const __m128 b0 = _mm_set_ps(b.mMax.y, a.mMax.y, b.mMax.z, a.mMax.z);
			const __m128 b1 = _mm_set_ps(a.mMin.y, b.mMin.y, a.mMin.z, b.mMin.z);
			const __m128 d = _mm_cmplt_ps(b0, b1);
			if(!_mm_movemask_ps(d))
				pairs.Add(RIndex0).Add(Remap1[Index1]);

			Index1++;
		}
		Index0++;
	}

	////

	Index0 = 0;
	udword RunningAddress0 = 0;
	while(RunningAddress0<nb0 && Index0<nb1)
	{
		const AABB& Box1 = BoxList1[Index0];

		const float MinLimit = Box1.mMin.x;
		while(BoxList0[RunningAddress0].mMin.x<=MinLimit)
			RunningAddress0++;

		const udword RIndex1 = Remap1[Index0];
		udword Index1 = RunningAddress0;
		const float MaxLimit = Box1.mMax.x;
		while(BoxList0[Index1].mMin.x<=MaxLimit)
		{
			const AABB& a = BoxList0[Index1];
			const AABB& b = Box1;
			const __m128 b0 = _mm_set_ps(b.mMax.y, a.mMax.y, b.mMax.z, a.mMax.z);
			const __m128 b1 = _mm_set_ps(a.mMin.y, b.mMin.y, a.mMin.z, b.mMin.z);
			const __m128 d = _mm_cmplt_ps(b0, b1);
			if(!_mm_movemask_ps(d))
				pairs.Add(Remap0[Index1]).Add(RIndex1);

			Index1++;
		}
		Index0++;
	}

	DELETEARRAY(BoxList1);
	DELETEARRAY(BoxList0);
#else
	SIMD_AABB* BoxList0 = new SIMD_AABB[nb0+1];
	SIMD_AABB* BoxList1 = new SIMD_AABB[nb1+1];
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
			BoxList0[i].InitFrom(list0[Remap0[i]]);
		BoxList0[nb0].mMinX = FLT_MAX;

		for(udword i=0;i<nb1;i++)
			BoxList1[i].InitFrom(list1[Remap1[i]]);
		BoxList1[nb1].mMinX = FLT_MAX;

		DELETEARRAY(PosList1);
		DELETEARRAY(PosList0);
	}

	// 3) Prune the lists
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const SIMD_AABB& Box0 = BoxList0[Index0];

		const float MinLimit = Box0.mMinX;
		while(BoxList1[RunningAddress1].mMinX<MinLimit)
			RunningAddress1++;

		SIMD_OVERLAP_INIT(Box0)

		const udword RIndex0 = Remap0[Index0];
		udword Index1 = RunningAddress1;
		const float MaxLimit = Box0.mMaxX;
		while(BoxList1[Index1].mMinX<=MaxLimit)
		{
			SIMD_OVERLAP_TEST(BoxList1[Index1])
				pairs.Add(RIndex0).Add(Remap1[Index1]);

			Index1++;
		}
		Index0++;
	}

	////

	Index0 = 0;
	udword RunningAddress0 = 0;
	while(RunningAddress0<nb0 && Index0<nb1)
	{
		const SIMD_AABB& Box1 = BoxList1[Index0];

		const float MinLimit = Box1.mMinX;
		while(BoxList0[RunningAddress0].mMinX<=MinLimit)
			RunningAddress0++;

		SIMD_OVERLAP_INIT(Box1)

		const udword RIndex1 = Remap1[Index0];
		udword Index1 = RunningAddress0;
		const float MaxLimit = Box1.mMaxX;
		while(BoxList0[Index1].mMinX<=MaxLimit)
		{
			SIMD_OVERLAP_TEST(BoxList0[Index1])
				pairs.Add(Remap0[Index1]).Add(RIndex1);

			Index1++;
		}
		Index0++;
	}

	DELETEARRAY(BoxList1);
	DELETEARRAY(BoxList0);
#endif

	return true;
}


#ifdef NAIVE_VERSION
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Complete box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to the same set.
 *	\param		nb		[in] number of boxes
 *	\param		list	[in] list of boxes
 *	\param		pairs	[out] list of overlapping pairs
 *	\param		axes	[in] projection order (0,2,1 is often best)
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	AABB* BoxList = new AABB[nb+1];
	udword* Remap;
	{
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
			BoxList[i] = list[Remap[i]];
		BoxList[nb].mMin.x = FLT_MAX;
		DELETEARRAY(PosList);
	}

	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const AABB& Box0 = BoxList[Index0];

		const float MinLimit = Box0.mMin.x;
		while(BoxList[RunningAddress++].mMin.x<MinLimit);

		udword Index1 = RunningAddress;
		const float MaxLimit = Box0.mMax.x;

		const udword RIndex0 = Remap[Index0];

		while(BoxList[Index1].mMin.x<=MaxLimit)
		{
			if(Index0!=Index1)
			{
				// MODIFIED
				const AABB& a = BoxList[Index0];
				const AABB& b = BoxList[Index1];
				const __m128 b0 = _mm_set_ps(b.mMax.y, a.mMax.y, b.mMax.z, a.mMax.z);
				const __m128 b1 = _mm_set_ps(a.mMin.y, b.mMin.y, a.mMin.z, b.mMin.z);
				const __m128 d = _mm_cmplt_ps(b0, b1);
				if(!_mm_movemask_ps(d))
//				if(intersects2D(BoxList[Index0], BoxList[Index1]))
					pairs.Add(RIndex0).Add(Remap[Index1]);
			}
			Index1++;
		}
		Index0++;
	}

	DELETEARRAY(BoxList);
	return true;
}

#else

#define ALLOC_DEFAULT			// 32000
//#define ALLOC_ALIGNED16		// 18000
//#define ALLOC_XP_ALIGNED4		// 33000
//#define ALLOC_XP_ALIGNED8		// 19000

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Complete box pruning. Returns a list of overlapping pairs of boxes, each box of the pair belongs to the same set.
 *	\param		nb		[in] number of boxes
 *	\param		list	[in] list of boxes
 *	\param		pairs	[out] list of overlapping pairs
 *	\param		axes	[in] projection order (0,2,1 is often best)
 *	\return		true if success.
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	// MODIFIED
#ifdef ALLOC_DEFAULT
	SIMD_AABB* BoxList = new SIMD_AABB[nb+1];
#endif
#ifdef ALLOC_ALIGNED16
	SIMD_AABB* BoxList = (SIMD_AABB*)_aligned_malloc(sizeof(SIMD_AABB)*(nb+1), 16);
#endif
#ifdef ALLOC_XP_ALIGNED4
	SIMD_AABB* BoxList = (SIMD_AABB*)(((char*)_aligned_malloc(sizeof(SIMD_AABB)*(nb+1+1), 16))+4);
#endif
#ifdef ALLOC_XP_ALIGNED8
	SIMD_AABB* BoxList = (SIMD_AABB*)(((char*)_aligned_malloc(sizeof(SIMD_AABB)*(nb+1+1), 16))+8);
#endif

	udword* Remap;
	{
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
			BoxList[i].InitFrom(list[Remap[i]]);	// MODIFIED
		BoxList[nb].mMinX = FLT_MAX;
		DELETEARRAY(PosList);
	}
//CPUID
	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const SIMD_AABB& Box0 = BoxList[Index0];	// MODIFIED

		const float MinLimit = Box0.mMinX;
		while(BoxList[RunningAddress++].mMinX<MinLimit);

		SIMD_OVERLAP_INIT(Box0)	// NEW

		udword Index1 = RunningAddress;
		const float MaxLimit = Box0.mMaxX;

		const udword RIndex0 = Remap[Index0];

		while(BoxList[Index1].mMinX<=MaxLimit)
		{
			if(Index0!=Index1)
			{
				ASSERT(Index0<Index1);

				// MODIFIED
				SIMD_OVERLAP_TEST(BoxList[Index1])
//					udword id0 = Index0;
//					udword id1 = Index1;
//					if(Remap[Index1]<Remap[Index0])
//						TSwap(id0, id1);
//				if(intersects2D2(BoxList[id0], BoxList[id1]))
//				if(intersects2D2(BoxList[Index0], BoxList[Index1]))
//				if(intersects2D(BoxList[Index0], BoxList[Index1]))
					pairs.Add(RIndex0).Add(Remap[Index1]);
			}
			Index1++;
		}
		Index0++;
	}
//CPUID
#ifdef ALLOC_DEFAULT
	DELETEARRAY(BoxList);
#endif
#ifdef ALLOC_ALIGNED16
	_aligned_free(BoxList);
#endif
#ifdef ALLOC_XP_ALIGNED4
	_aligned_free(BoxList);
#endif
#ifdef ALLOC_XP_ALIGNED8
	_aligned_free(BoxList);
#endif
	return true;
}
#endif
