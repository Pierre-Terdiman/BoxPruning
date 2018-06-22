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
#define USE_INTEGER_XS

#ifdef USE_INTEGER_XS
	typedef udword			XType;
	#define SENTINEL_VALUE	0xffffffff
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
#else
	typedef float			XType;
	#define SENTINEL_VALUE	FLT_MAX
#endif

struct SIMD_AABB_X
{
	__forceinline	SIMD_AABB_X()	{}
	__forceinline	~SIMD_AABB_X()	{}

	void	InitFrom(const AABB& b)
	{
#ifdef USE_INTEGER_XS
		const udword* Binary = reinterpret_cast<const udword*>(&b.mMin.x);
		mMinX	= encodeFloat(Binary[0]);
		mMaxX	= encodeFloat(Binary[3]);
#else
		mMinX	= b.mMin.x;
		mMaxX	= b.mMax.x;
#endif
	}

	// NEW
	__forceinline	void		operator = (const SIMD_AABB_X& box)
	{
		mMinX = box.mMinX;
		mMaxX = box.mMaxX;
	}

    XType mMinX;
    XType mMaxX;
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

	// NEW
	__forceinline	void		operator = (const SIMD_AABB_YZ& box)
	{
		_mm_store_ps(&mMinY, _mm_load_ps(&box.mMinY));
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

static void /*__cdecl*/ outputPair(udword id0, udword id1, Pairs& pairs, const udword* remap)
{
	pairs.AddPair(id0, remap[id1]);
}

// NEW
static const ubyte gCodes[] = {	4, 4, 4, 255, 4, 3, 2, 255,
								4, 1, 0, 255, 255, 255, 255, 255 };

static __forceinline udword classifyBoxNew(const AABB& box, const float limitY, const float limitZ)
{
	const bool lowerPart = box.mMax.z < limitZ;
	const bool upperPart = box.mMin.z > limitZ;
	const bool leftPart = box.mMax.y < limitY;
	const bool rightPart = box.mMin.y > limitY;

	// Table-based box classification avoids many branches
	const udword Code = udword(rightPart)|(udword(leftPart)<<1)|(udword(upperPart)<<2)|(udword(lowerPart)<<3);
	assert(gCodes[Code]!=255);
	return gCodes[Code];
}

template<int codepath>
static void bipartiteKernel(udword nb0, udword nb1,
							const SIMD_AABB_X* BoxListX0, const SIMD_AABB_X* BoxListX1,
							const SIMD_AABB_YZ* BoxListYZ0, const SIMD_AABB_YZ* BoxListYZ1,
							const udword* Remap0, const udword* Remap1, Pairs& pairs)
{
	udword Index0 = 0;
	udword RunningAddress1 = 0;
	while(RunningAddress1<nb1 && Index0<nb0)
	{
		const SIMD_AABB_X& Box0X = BoxListX0[Index0];

		const XType MinLimit = Box0X.mMinX;
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
		const XType MaxLimit = Box0X.mMaxX;

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ1[RunningAddress1];
		const char* const CurrentBoxListX = (const char*)&BoxListX1[RunningAddress1];

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
		pairs.AddPair(RIndex0, Remap1[Index1]);
		}
		_asm	align 16
StartLoop4:
		while(*(const XType*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0)
			BLOCK4(8, FoundOverlap1)
			BLOCK4(16, FoundOverlap2)
			BLOCK4(24, FoundOverlap3)
			Offset += 40;
			BLOCK4(-8, FoundOverlap)
		}
#undef BLOCK4

#define BLOCK	if(*(const XType*)(CurrentBoxListX + Offset)<=MaxLimit)			\
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

		Index0++;
	}
}

// MODIFIED (moved to separate function)
static void DoBipartiteBoxPruning(	udword nb0, const SIMD_AABB_X* BoxListX0, const SIMD_AABB_YZ* BoxListYZ0, const udword* Remap0,
									udword nb1, const SIMD_AABB_X* BoxListX1, const SIMD_AABB_YZ* BoxListYZ1, const udword* Remap1, Pairs& pairs)
{
	if(!nb0 || !nb1)
		return;

	bipartiteKernel<0>(nb0, nb1, BoxListX0, BoxListX1, BoxListYZ0, BoxListYZ1, Remap0, Remap1, pairs);
	bipartiteKernel<1>(nb1, nb0, BoxListX1, BoxListX0, BoxListYZ1, BoxListYZ0, Remap1, Remap0, pairs);
}

// MODIFIED (moved to separate function)
static void DoCompleteBoxPruning(udword nb, const SIMD_AABB_X* BoxListX, const SIMD_AABB_YZ* BoxListYZ, const udword* Remap, Pairs& pairs)
{
	if(!nb)
		return;

	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const SIMD_AABB_X& Box0X = BoxListX[Index0];

		const XType MinLimit = Box0X.mMinX;
		while(BoxListX[RunningAddress++].mMinX<MinLimit);

		const SIMD_AABB_YZ& Box0YZ = BoxListYZ[Index0];
		SIMD_OVERLAP_INIT(Box0YZ)

		const XType MaxLimit = Box0X.mMaxX;
		const udword RIndex0 = Remap[Index0];

		udword Offset = 0;
		const char* const CurrentBoxListYZ = (const char*)&BoxListYZ[RunningAddress];
		const char* const CurrentBoxListX = (const char*)&BoxListX[RunningAddress];

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
		pairs.AddPair(RIndex0, Remap[Index]);
		}
		_asm	align 16
StartLoop4:
		while(*(const XType*)(CurrentBoxListX + Offset + 8*5)<=MaxLimit)
		{
			BLOCK4(0, FoundOverlap0)
			BLOCK4(8, FoundOverlap1)
			BLOCK4(16, FoundOverlap2)
			BLOCK4(24, FoundOverlap3)
			Offset += 40;
			BLOCK4(-8, FoundOverlap)
		}

#define BLOCK	if(*(const XType*)(CurrentBoxListX + Offset)<=MaxLimit)			\
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

		Index0++;
	}
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
#define NB_BUCKETS				5
#define NB_SENTINEL_PER_BUCKET	6
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs_)
{
	// Checkings
	if(!nb || !list)
		return false;

	// MODIFIED
	SIMD_AABB_X* BoxListXBuffer = new SIMD_AABB_X[nb+NB_SENTINEL_PER_BUCKET*NB_BUCKETS];
	SIMD_AABB_YZ* BoxListYZBuffer = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*nb, 16);

	float* PosList = new float[nb+1];

	// NEW
	__m128 minV = _mm_loadu_ps(&list[0].mMin.x);
	__m128 maxV = _mm_loadu_ps(&list[0].mMax.x);
	PosList[0] = list[0].mMin.x;
	// MODIFIED
	for(udword i=1;i<nb;i++)
	{
		PosList[i] = list[i].mMin.x;
		minV = _mm_min_ps(minV, _mm_loadu_ps(&list[i].mMin.x));
		maxV = _mm_max_ps(maxV, _mm_loadu_ps(&list[i].mMax.x));
	}
	PosList[nb] = FLT_MAX;

	// NEW
	__declspec(align(16)) float mergedMin[4];
	__declspec(align(16)) float mergedMax[4];
	_mm_store_ps(mergedMin, minV);
	_mm_store_ps(mergedMax, maxV);

	static PRUNING_SORTER RS;	// Static for coherence
	const udword* Sorted = RS.Sort(PosList, nb+1).GetRanks();

	// NEW
	const float limitY = (mergedMax[1] + mergedMin[1]) * 0.5f;
	const float limitZ = (mergedMax[2] + mergedMin[2]) * 0.5f;

	udword Counters[NB_BUCKETS];
	for(udword i=0;i<NB_BUCKETS;i++)
		Counters[i] = 0;

	udword* Indices = reinterpret_cast<udword*>(PosList);
	{
		const AABB* boxes = list;
		udword* dst = Indices;
		udword NbToGo = nb;
		while(NbToGo--)
		{
			const udword index = classifyBoxNew(*boxes++, limitY, limitZ);
			*dst++ = index;
			Counters[index]++;
		}
	}

	udword* Remap = RS.GetRecyclable();
	SIMD_AABB_X* BoxListX[NB_BUCKETS];
	SIMD_AABB_YZ* BoxListYZ[NB_BUCKETS];
	udword* RemapBase[NB_BUCKETS];
	{
		SIMD_AABB_X* CurrentBoxListXBuffer = BoxListXBuffer;
		SIMD_AABB_YZ* CurrentBoxListYZBuffer = BoxListYZBuffer;
		udword* CurrentRemap = Remap;
		for(udword i=0;i<NB_BUCKETS;i++)
		{
			const udword Nb = Counters[i];
			BoxListX[i] = CurrentBoxListXBuffer;
			BoxListYZ[i] = CurrentBoxListYZBuffer;
			RemapBase[i] = CurrentRemap;
			CurrentBoxListXBuffer += Nb+NB_SENTINEL_PER_BUCKET;
			CurrentBoxListYZBuffer += Nb;
			CurrentRemap += Nb;
		}
		assert(CurrentBoxListXBuffer == BoxListXBuffer + nb + NB_SENTINEL_PER_BUCKET*NB_BUCKETS);
		assert(CurrentBoxListYZBuffer == BoxListYZBuffer + nb);
		assert(CurrentRemap == Remap + nb);
	}

	for(udword i=0;i<NB_BUCKETS;i++)
		Counters[i] = 0;

	udword NbToGo = nb;
	while(NbToGo--)
	{
		const udword SortedIndex = *Sorted++;
		const udword TargetBucket = Indices[SortedIndex];
		const udword IndexInTarget = Counters[TargetBucket]++;

		SIMD_AABB_X* TargetBoxListX = BoxListX[TargetBucket];
		SIMD_AABB_YZ* TargetBoxListYZ = BoxListYZ[TargetBucket];
		udword* TargetRemap = RemapBase[TargetBucket];

		const AABB& SrcBox = list[SortedIndex];
		TargetBoxListX[IndexInTarget].InitFrom(SrcBox);
		TargetBoxListYZ[IndexInTarget].InitFrom(SrcBox);
		TargetRemap[IndexInTarget] = SortedIndex;
	}

	for(udword i=0;i<NB_BUCKETS;i++)
	{
		SIMD_AABB_X* TargetBoxListX = BoxListX[i];
		const udword IndexInTarget = Counters[i];
		for(udword j=0;j<NB_SENTINEL_PER_BUCKET;j++)
			TargetBoxListX[IndexInTarget+j].mMinX = SENTINEL_VALUE;
	}

	DELETEARRAY(PosList);

	Pairs pairs(pairs_);

	for(udword i=0;i<NB_BUCKETS;i++)
	{
		DoCompleteBoxPruning(Counters[i], BoxListX[i], BoxListYZ[i], RemapBase[i], pairs);
	}

	const udword LastBucket = NB_BUCKETS-1;
	for(udword i=0;i<LastBucket;i++)
	{
		DoBipartiteBoxPruning(	Counters[i], BoxListX[i], BoxListYZ[i], RemapBase[i],
								Counters[LastBucket], BoxListX[LastBucket], BoxListYZ[LastBucket], RemapBase[LastBucket],
								pairs);
	}

	_aligned_free(BoxListYZBuffer);
	DELETEARRAY(BoxListXBuffer);

	return true;
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

#define NB_BUCKETS	5
#define USE_BUCKET_BOUNDS

	class Buckets
	{
		public:
			Buckets(udword nb, const AABB* list, RadixSort& rs);
			~Buckets();

			SIMD_AABB_X*	mBoxListXBuffer;
			SIMD_AABB_YZ*	mBoxListYZBuffer;
			SIMD_AABB_X*	mBoxListX[NB_BUCKETS];
			SIMD_AABB_YZ*	mBoxListYZ[NB_BUCKETS];
			udword			mOffsets[NB_BUCKETS];
			udword			mCounters[NB_BUCKETS];
#ifdef USE_BUCKET_BOUNDS
			struct Bounds
			{
				__m128	mMinV;
				__m128	mMaxV;
			};
			__declspec(align(16)) Bounds	mBucketBounds[NB_BUCKETS];
#endif
	};

Buckets::Buckets(udword nb, const AABB* list, RadixSort& rs)
{
	udword* Remap;

	// Counters will add up to 'nb', then we need +1 for each bucket, so that's +5 overall
	SIMD_AABB_X* BoxListXBuffer = new SIMD_AABB_X[nb+NB_SENTINEL_PER_BUCKET*NB_BUCKETS];
	SIMD_AABB_YZ* BoxListYZBuffer = (SIMD_AABB_YZ*)_aligned_malloc(sizeof(SIMD_AABB_YZ)*nb, 16);

	float* PosList = new float[nb+1];

	__m128 minV = _mm_loadu_ps(&list[0].mMin.x);
	__m128 maxV = _mm_loadu_ps(&list[0].mMax.x);
	PosList[0] = list[0].mMin.x;
	for(udword i=1;i<nb;i++)
	{
		PosList[i] = list[i].mMin.x;
		minV = _mm_min_ps(minV, _mm_loadu_ps(&list[i].mMin.x));
		maxV = _mm_max_ps(maxV, _mm_loadu_ps(&list[i].mMax.x));
	}
	PosList[nb] = FLT_MAX;

	__declspec(align(16)) float mergedMin[4];
	__declspec(align(16)) float mergedMax[4];
	_mm_store_ps(mergedMin, minV);
	_mm_store_ps(mergedMax, maxV);

	const udword* Sorted = rs.Sort(PosList, nb+1).GetRanks();

	const float limitY = (mergedMax[1] + mergedMin[1]) * 0.5f;
	const float limitZ = (mergedMax[2] + mergedMin[2]) * 0.5f;

	udword Counters[NB_BUCKETS];
	for(udword i=0;i<NB_BUCKETS;i++)
		Counters[i] = 0;

	udword* Indices = reinterpret_cast<udword*>(PosList);
	for(udword i=0;i<nb;i++)
	{
		const udword index = classifyBoxNew(list[i], limitY, limitZ);
		Indices[i] = index;
		Counters[index]++;
	}

	SIMD_AABB_X* CurrentBoxListXBuffer = BoxListXBuffer;
	SIMD_AABB_YZ* CurrentBoxListYZBuffer = BoxListYZBuffer;
	SIMD_AABB_X* BoxListX[NB_BUCKETS];
	SIMD_AABB_YZ* BoxListYZ[NB_BUCKETS];
	for(udword i=0;i<NB_BUCKETS;i++)
	{
		const udword Nb = Counters[i];
		BoxListX[i] = mBoxListX[i] = CurrentBoxListXBuffer;
		BoxListYZ[i] = mBoxListYZ[i] = CurrentBoxListYZBuffer;
		CurrentBoxListXBuffer += Nb+NB_SENTINEL_PER_BUCKET;
		CurrentBoxListYZBuffer += Nb;
	}
	assert(CurrentBoxListXBuffer == BoxListXBuffer + nb + NB_SENTINEL_PER_BUCKET*NB_BUCKETS);
	assert(CurrentBoxListYZBuffer == BoxListYZBuffer + nb);

	udword Offsets[NB_BUCKETS];
	Offsets[0]=0;
	for(udword i=0;i<NB_BUCKETS-1;i++)
		Offsets[i+1] = Offsets[i] + Counters[i];

	Remap = rs.GetRecyclable();

	for(udword i=0;i<NB_BUCKETS;i++)
		Counters[i] = 0;

#ifdef USE_BUCKET_BOUNDS
	const float maxFloat = FLT_MAX;
	const float minFloat = -FLT_MAX;
	const __m128 minV0 = _mm_load1_ps(&maxFloat);
	const __m128 maxV0 = _mm_load1_ps(&minFloat);
	for(udword i=0;i<NB_BUCKETS;i++)
	{
		mBucketBounds[i].mMinV = minV0;
		mBucketBounds[i].mMaxV = maxV0;
	}
#endif

	for(udword i=0;i<nb;i++)
	{
		const udword SortedIndex = Sorted[i];
		const udword TargetBucket = Indices[SortedIndex];

		const udword CurrentOffset = Offsets[TargetBucket]++;
		Remap[CurrentOffset] = SortedIndex;

		SIMD_AABB_X* TargetBoxListX = BoxListX[TargetBucket];
		SIMD_AABB_YZ* TargetBoxListYZ = BoxListYZ[TargetBucket];

		const udword IndexInTarget = Counters[TargetBucket]++;

		const AABB& SrcBox = list[SortedIndex];
		TargetBoxListX[IndexInTarget].InitFrom(SrcBox);
		TargetBoxListYZ[IndexInTarget].InitFrom(SrcBox);
#ifdef USE_BUCKET_BOUNDS
		mBucketBounds[TargetBucket].mMinV = _mm_min_ps(mBucketBounds[TargetBucket].mMinV, _mm_loadu_ps(&SrcBox.mMin.x));
		mBucketBounds[TargetBucket].mMaxV = _mm_max_ps(mBucketBounds[TargetBucket].mMaxV, _mm_loadu_ps(&SrcBox.mMax.x));
#endif
	}

	for(udword i=0;i<NB_BUCKETS;i++)
	{
		SIMD_AABB_X* TargetBoxListX = BoxListX[i];
		const udword IndexInTarget = Counters[i];
		for(udword j=0;j<NB_SENTINEL_PER_BUCKET;j++)
			TargetBoxListX[IndexInTarget+j].mMinX = SENTINEL_VALUE;
	}
	Offsets[0]=0;
	for(udword i=0;i<NB_BUCKETS-1;i++)
		Offsets[i+1] = Offsets[i] + Counters[i];

	DELETEARRAY(PosList);

	mBoxListXBuffer = BoxListXBuffer;
	mBoxListYZBuffer = BoxListYZBuffer;

	for(udword i=0;i<NB_BUCKETS;i++)
	{
		mOffsets[i] = Offsets[i];
		mCounters[i] = Counters[i];
	}
}

Buckets::~Buckets()
{
	_aligned_free(mBoxListYZBuffer);
	DELETEARRAY(mBoxListXBuffer);
}

bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB* list0, udword nb1, const AABB* list1, Container& pairs_)
{
	// Checkings
	if(!nb0 || !list0 || !nb1 || !list1)
		return false;

	static PRUNING_SORTER RS0, RS1;	// Static for coherence.
	Buckets Buckets0(nb0, list0, RS0);
	Buckets Buckets1(nb1, list1, RS1);

	Pairs pairs(pairs_);

	const udword* Remap0 = RS0.GetRecyclable();
	const udword* Remap1 = RS1.GetRecyclable();
	for(udword i=0;i<NB_BUCKETS;i++)
	{
		const udword Nb0 = Buckets0.mCounters[i];
		const SIMD_AABB_X* BoxListX0 = Buckets0.mBoxListX[i];
		const SIMD_AABB_YZ* BoxListYZ0 = Buckets0.mBoxListYZ[i];
		const udword* Remap0_ = Remap0 + Buckets0.mOffsets[i];

#ifdef USE_BUCKET_BOUNDS
		__declspec(align(16)) float mergedMin0[4];
		__declspec(align(16)) float mergedMax0[4];
		_mm_store_ps(mergedMin0, Buckets0.mBucketBounds[i].mMinV);
		_mm_store_ps(mergedMax0, Buckets0.mBucketBounds[i].mMaxV);

		AABB Bounds0;
		Bounds0.mMin.x = mergedMin0[0];
		Bounds0.mMin.y = mergedMin0[1];
		Bounds0.mMin.z = mergedMin0[2];
		Bounds0.mMax.x = mergedMax0[0];
		Bounds0.mMax.y = mergedMax0[1];
		Bounds0.mMax.z = mergedMax0[2];
#endif

		for(udword j=0;j<NB_BUCKETS;j++)
		{
#ifdef USE_BUCKET_BOUNDS
			__declspec(align(16)) float mergedMin1[4];
			__declspec(align(16)) float mergedMax1[4];
			_mm_store_ps(mergedMin1, Buckets1.mBucketBounds[j].mMinV);
			_mm_store_ps(mergedMax1, Buckets1.mBucketBounds[j].mMaxV);

			AABB Bounds1;
			Bounds1.mMin.x = mergedMin1[0];
			Bounds1.mMin.y = mergedMin1[1];
			Bounds1.mMin.z = mergedMin1[2];
			Bounds1.mMax.x = mergedMax1[0];
			Bounds1.mMax.y = mergedMax1[1];
			Bounds1.mMax.z = mergedMax1[2];

			if(Bounds0.Intersect(Bounds1))
#endif
			{
				const udword Nb1 = Buckets1.mCounters[j];
				const SIMD_AABB_X* BoxListX1 = Buckets1.mBoxListX[j];
				const SIMD_AABB_YZ* BoxListYZ1 = Buckets1.mBoxListYZ[j];
				const udword* Remap1_ = Remap1 + Buckets1.mOffsets[j];

				DoBipartiteBoxPruning(	Nb0, BoxListX0, BoxListYZ0, Remap0_,
										Nb1, BoxListX1, BoxListYZ1, Remap1_,
										pairs);
			}
		}
	}

	return true;
}
