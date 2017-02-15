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

static __forceinline int intersects2D(const AABB& a, const AABB& b)
{
       if(    b.mMax.y < a.mMin.y || a.mMax.y < b.mMin.y
       ||     b.mMax.z < a.mMin.z || a.mMax.z < b.mMin.z)
              return 0;
       return 1;
}

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
			if(intersects2D(Box0, BoxList1[Index1]))
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
			if(intersects2D(BoxList0[Index1], Box1))
				pairs.Add(Remap0[Index1]).Add(RIndex1);

			Index1++;
		}
		Index0++;
	}

	DELETEARRAY(BoxList1);
	DELETEARRAY(BoxList0);

	return true;
}



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
		mMinY	= b.mMin.y;
		mMinZ	= b.mMin.z;
		mMaxY	= b.mMax.y;
		mMaxZ	= b.mMax.z;
	}

    float mMinY;
    float mMinZ;
    float mMaxY;
    float mMaxZ;
};

/*static __forceinline int intersects2D(const SIMD_AABB_YZ& a, const SIMD_AABB_YZ& b)
{
       if(    b.mMaxY <= a.mMinY || a.mMaxY < b.mMinY
       ||     b.mMaxZ <= a.mMinZ || a.mMaxZ < b.mMinZ)
              return 0;
       return 1;
}*/

	#include <xmmintrin.h>
	#include <emmintrin.h>

/*#define SIMD_OVERLAP_INIT(box)	\
       const __m128 b = _mm_shuffle_ps(_mm_load_ps(&box.mMinY), _mm_load_ps(&box.mMinY), 78);

#define SIMD_OVERLAP_TEST(box)						\
       const __m128 a = _mm_load_ps(&box.mMinY);	\
       const __m128 d = _mm_cmpgt_ps(a, b);			\
       if(_mm_movemask_ps(d)==12)*/

static void __cdecl outputPair(udword id0, udword id1, Container& pairs, const udword* remap)
{
	pairs.Add(remap[id0]).Add(remap[id1]);
}

static void __cdecl outputPair2(udword id0, udword id1, Container& pairs, const udword* remap)
{
	pairs.Add(id0).Add(remap[id1]);
}

/*static void debug(udword id0, udword id1, Container& pairs, const udword* remap)
{
	if(id0==9 && remap[id1]==6)
		int FoundIt=0;
}*/


static void Error()
{
	printf("ERROR!\n");
}

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

//#define NAIVE_VERSION
#ifdef NAIVE_VERSION
		const float MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		__m128	SavedXMM1;
		__m128	SavedXMM2;

		_asm
		{
			movss		xmm1, MaxLimit				// xmm1 = MaxLimit
			mov			edx, BoxListYZ				// edx = BoxListYZ

			mov			edi, Index0					// edi = Index0
			add			edi, edi					// edi = Index0 * 2
			lea			esi, dword ptr [edx+edi*8]	// esi = BoxListYZ + Index0*16 = &BoxListYZ[Index0] (*16 because sizeof(SIMD_AABB_YZ)==16)
			movaps		xmm2, xmmword ptr [esi]		// xmm2 = BoxListYZ[Index0] = Box0YZ <= that's the _mm_load_ps in SIMD_OVERLAP_INIT
			shufps		xmm2, xmm2, 4Eh				// that's our _mm_shuffle_ps in SIMD_OVERLAP_INIT

			mov			edi, RunningAddress			// edi = Index1
			mov			eax, edi					// eax = Index1
			shl			eax, 4						// eax = Index1 * 16
			add			edx, eax					// edx = &BoxListYZ[Index1]

			mov			eax, BoxListX				// eax = BoxListX
			lea			esi, dword ptr [eax+edi*8]	// esi = BoxListX + Index1*8 = &BoxListX[Index1] (*8 because sizeof(SIMD_AABB_X)==8)
			comiss		xmm1, xmmword ptr [esi]		// [esi] = BoxListX[Index1].mMinX, compared to MaxLimit
			jb			ExitLoop

			align		16							// Align start of loop on 16-byte boundary for perf
EnterLoop:
			movaps		xmm3, xmm2					// xmm2 = pre-shuffled box, constant for the duration of the loop
			cmpltps		xmm3, xmmword ptr [edx]		// that's _mm_cmpgt_ps in SIMD_OVERLAP_TEST
			movmskps	eax, xmm3					// that's _mm_movemask_ps in SIMD_OVERLAP_TEST

			cmp			eax, 0Ch
			jne			NoOverlap

			movaps		SavedXMM1, xmm1
			movaps		SavedXMM2, xmm2
			pushad
				push		Remap
				push		pairs
				push		edi
				push		RIndex0
				call		outputPair2;
				add         esp, 16
			popad
			movaps		xmm1, SavedXMM1
			movaps		xmm2, SavedXMM2

			align		16
NoOverlap:;
			inc			edi							// Index1++
			add			esi, 8
			add			edx, 16
			comiss		xmm1, xmmword ptr [esi]		// [esi] = BoxListX[Index1].mMinX, compared to MaxLimit
			jae			EnterLoop
ExitLoop:;
		}
#endif

#define VERSION4
#ifdef VERSION4
		const float MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		__m128	SavedXMM1;
		__m128	SavedXMM2;

		_asm
		{
			movss		xmm1, MaxLimit				// xmm1 = MaxLimit
			mov			edx, BoxListYZ				// edx = BoxListYZ

			mov			edi, Index0					// edi = Index0
			add			edi, edi					// edi = Index0 * 2
			lea			esi, dword ptr [edx+edi*8]	// esi = BoxListYZ + Index0*16 = &BoxListYZ[Index0] (*16 because sizeof(SIMD_AABB_YZ)==16)
			movaps		xmm2, xmmword ptr [esi]		// xmm2 = BoxListYZ[Index0] = Box0YZ <= that's the _mm_load_ps in SIMD_OVERLAP_INIT
			shufps		xmm2, xmm2, 4Eh				// that's our _mm_shuffle_ps in SIMD_OVERLAP_INIT

			mov			edi, RunningAddress			// edi = Index1
			mov			eax, edi					// eax = Index1
			shl			eax, 4						// eax = Index1 * 16
			add			edx, eax					// edx = &BoxListYZ[Index1]

			mov			eax, BoxListX				// eax = BoxListX
			lea			esi, dword ptr [eax+edi*8]	// esi = BoxListX + Index1*8 = &BoxListX[Index1] (*8 because sizeof(SIMD_AABB_X)==8)
			comiss		xmm1, xmmword ptr [esi]		// [esi] = BoxListX[Index1].mMinX, compared to MaxLimit
			jb			ExitLoop

			xor			ecx, ecx

			align		16							// Align start of loop on 16-byte boundary for perf
EnterLoop:
			// ~11600 with this:
			movaps		xmm3, xmmword ptr [edx+ecx*2]		// Box1YZ
			cmpnltps	xmm3, xmm2

			// ~11900 with this:
//			movaps		xmm3, xmm2
//			cmpltps		xmm3, xmmword ptr [edx+ecx*2]

//			movaps		xmm3, xmm2
//			subps		xmm3, xmmword ptr [edx+ecx*2]		// Box1YZ

			movmskps	eax, xmm3

			cmp			eax, 0Ch
			jne			NoOverlap

			movaps		SavedXMM1, xmm1
			movaps		SavedXMM2, xmm2
			pushad
				// Recompute Index1
				add			ecx, esi
				sub			ecx, BoxListX
				shr			ecx, 3

				push		Remap
				push		pairs
				push		ecx
//				push		Index0
//				call		outputPair;
				push		RIndex0
				call		outputPair2;
				add         esp, 16

			popad
			movaps		xmm1, SavedXMM1
			movaps		xmm2, SavedXMM2

			align		16
NoOverlap:;
			add			ecx, 8
			comiss		xmm1, xmmword ptr [esi+ecx]		// [esi] = BoxListX[Index1].mMinX, compared to MaxLimit
			jae			EnterLoop
ExitLoop:;
		}
#endif
		
		Index0++;
	}

	_aligned_free(BoxListYZ);
	DELETEARRAY(BoxListX);
	return true;
}

