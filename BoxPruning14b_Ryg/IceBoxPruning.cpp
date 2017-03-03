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


// Munge the float bits to return produce an unsigned order-preserving
// ranking of floating-point numbers.
// (Old trick: http://stereopsis.com/radix.html FloatFlip, with a new
// spin to get rid of -0.0f)

// In /fp:precise, we can just calc "x + 0.0f" and get what we need.
// But fast math optimizes it away. Could use #pragma float_control,
// but that prohibits inlining of MungeFloat. So do this silly thing
// instead.
float g_global_this_always_zero = 0.0f;

static inline udword MungeFloat(float f)
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
//#pragma float_control(pop)

struct SIMD_AABB_X
{
	__forceinline	SIMD_AABB_X()	{}
	__forceinline	~SIMD_AABB_X()	{}

	void	InitFrom(const AABB& b)
	{
		mMinX	= MungeFloat(b.mMin.x);
		mMaxX	= MungeFloat(b.mMax.x);
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

	SIMD_AABB_X* BoxListX = new SIMD_AABB_X[nb+5];
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
		BoxListX[nb+0].mMinX = ~0u;
		BoxListX[nb+1].mMinX = ~0u;
		BoxListX[nb+2].mMinX = ~0u;
		BoxListX[nb+3].mMinX = ~0u;
		BoxListX[nb+4].mMinX = ~0u;
		DELETEARRAY(PosList);
//	}

	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const SIMD_AABB_X& Box0X = BoxListX[Index0];

		const udword MinLimit = Box0X.mMinX;
		while(BoxListX[RunningAddress++].mMinX<MinLimit);

        if(BoxListX[RunningAddress].mMinX > Box0X.mMaxX)
        {
            Index0++;
            continue;
        }

#define VERSION4
#ifdef VERSION4
		const udword MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		__m128	SavedXMM1;
		__m128	SavedXMM2;

		_asm
		{
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

            mov         edi, MaxLimit               // edi = MaxLimit
			xor			ecx, ecx

			align		16							// Align start of loop on 16-byte boundary for perf
FastLoop:
            cmp         edi, [esi+ecx+24]           // [esi] = BoxListX[Index1].mMinX, compared to MaxLimit - can safely do another 4 iters of this?
            jb          CarefulLoop // nope!

            // Unroll 0
            movaps      xmm3, xmmword ptr [edx+ecx*2+0]  // Box1YZ
            cmpnleps    xmm3, xmm2
            movmskps    eax, xmm3
            cmp         eax, 0Ch
            je          FoundSlot0

            // Unroll 1
            movaps      xmm3, xmmword ptr [edx+ecx*2+16] // Box1YZ
            cmpnleps    xmm3, xmm2
            movmskps    eax, xmm3
            cmp         eax, 0Ch
            je          FoundSlot1

            // Unroll 2
            movaps      xmm3, xmmword ptr [edx+ecx*2+32]  // Box1YZ
            cmpnleps    xmm3, xmm2
            movmskps    eax, xmm3
            cmp         eax, 0Ch
            je          FoundSlot2

            // Unroll 3
            movaps      xmm3, xmmword ptr [edx+ecx*2+48]  // Box1YZ
            add         ecx, 32                           // Advance
            cmpnleps    xmm3, xmm2
            movmskps    eax, xmm3
            cmp         eax, 0Ch
            jne         FastLoop

            jmp         FastFoundOne
            // slots 2-0 fall through for increment (I'm being lazy...)
FoundSlot2:
            add         ecx, 8
FoundSlot1:
            add         ecx, 8
FoundSlot0:
            add         ecx, 8
FastFoundOne:
			movaps		SavedXMM1, xmm1
			movaps		SavedXMM2, xmm2
			push        eax
            push        ecx
            push        edx

				// Recompute Index1
                lea         ecx, [ecx+esi-8]
				sub			ecx, BoxListX
				shr			ecx, 3

				push		Remap
				push		pairs
				push		ecx
				push		RIndex0
				call		outputPair2;
				add         esp, 16

			pop         edx
            pop         ecx
            pop         eax
			movaps		xmm1, SavedXMM1
			movaps		xmm2, SavedXMM2
            jmp         FastLoop


CarefulLoop:
            cmp         edi, [esi+ecx]                   // [esi] = BoxListX[Index1].mMinX, compared to MaxLimit - can safely do another 4 iters of this?
            jb          ExitLoop

			// ~11600 with this:
			movaps		xmm3, xmmword ptr [edx+ecx*2]		// Box1YZ
            add         ecx, 8
			cmpnleps    xmm3, xmm2
			movmskps    eax, xmm3
			cmp			eax, 0Ch
			jne         CarefulLoop

         // found one!
			movaps		SavedXMM1, xmm1
			movaps		SavedXMM2, xmm2
			push        eax
            push        ecx
            push        edx

				// Recompute Index1
                lea         ecx, [ecx+esi-8]
				sub			ecx, BoxListX
				shr			ecx, 3

				push		Remap
				push		pairs
				push		ecx
				push		RIndex0
				call		outputPair2;
				add         esp, 16

			pop         edx
            pop         ecx
            pop         eax
			movaps		xmm1, SavedXMM1
			movaps		xmm2, SavedXMM2
            jmp         CarefulLoop

ExitLoop:;
		}
#endif
		
		Index0++;
	}

	_aligned_free(BoxListYZ);
	DELETEARRAY(BoxListX);
	return true;
}

