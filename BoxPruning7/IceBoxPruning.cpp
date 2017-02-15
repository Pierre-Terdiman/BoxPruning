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
/*       if(    b.mMax.y <= a.mMin.y || a.mMax.y < b.mMin.y
       ||     b.mMax.z <= a.mMin.z || a.mMax.z < b.mMin.z)
              return 0;*/
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

// NEW
static __forceinline udword encodeFloat(udword ir)
{
	if(ir & 0x80000000) //negative?
		return ~ir;//reverse sequence of negative numbers
	else
		return ir | 0x80000000; // flip sign
}

// NEW
struct IntegerAABB
{
	__forceinline	IntegerAABB()	{}
	__forceinline	~IntegerAABB()	{}

	void	InitFrom(const AABB& b)
	{
		const udword* Binary = reinterpret_cast<const udword*>(&b.mMin.x);
		mMinX	= encodeFloat(Binary[0]);
		mMinY	= encodeFloat(Binary[1]);
		mMinZ	= encodeFloat(Binary[2]);
		mMaxX	= encodeFloat(Binary[3]);
		mMaxY	= encodeFloat(Binary[4]);
		mMaxZ	= encodeFloat(Binary[5]);
	}

	udword	mMinX;
	udword	mMinY;
	udword	mMinZ;
	udword	mMaxX;
	udword	mMaxY;
	udword	mMaxZ;
};

// NEW
static __forceinline int intersects2D(const IntegerAABB& a, const IntegerAABB& b)
{
	if(1)
	{
		if(    b.mMaxY < a.mMinY || a.mMaxY < b.mMinY
		||     b.mMaxZ < a.mMinZ || a.mMaxZ < b.mMinZ)
			return 0;
		return 1;
	}
	else
	{
		const udword bits0 = ((b.mMaxY>>1) - (a.mMinY>>1))&0x80000000;
		const udword bits1 = ((b.mMaxZ>>1) - (a.mMinZ>>1))&0x80000000;
		const udword bits2 = ((a.mMaxY>>1) - (b.mMinY>>1))&0x80000000;
		const udword bits3 = ((a.mMaxZ>>1) - (b.mMinZ>>1))&0x80000000;
		const udword mask = bits0|(bits1>>1)|(bits2>>2)|(bits3>>3);
		return !mask;
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
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB* list, Container& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	IntegerAABB* BoxList = new IntegerAABB[nb+1];	// MODIFIED
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
		const float FltMax = FLT_MAX;
		BoxList[nb].mMinX = encodeFloat(IR(FltMax));	// MODIFIED
		DELETEARRAY(PosList);
	}

	// 3) Prune the list
	udword RunningAddress = 0;
	udword Index0 = 0;
	while(RunningAddress<nb && Index0<nb)
	{
		const IntegerAABB& Box0 = BoxList[Index0];	// MODIFIED

		const udword MinLimit = Box0.mMinX;	// MODIFIED
		while(BoxList[RunningAddress++].mMinX<MinLimit);	// NOW INTEGER CMP

		udword Index1 = RunningAddress;
		const udword MaxLimit = Box0.mMaxX;	// MODIFIED

		const udword RIndex0 = Remap[Index0];

		while(BoxList[Index1].mMinX<=MaxLimit)	// NOW INTEGER CMP
		{
			if(Index0!=Index1)
			{
				if(intersects2D(BoxList[Index0], BoxList[Index1]))	// NOW INTEGER CMP
					pairs.Add(RIndex0).Add(Remap[Index1]);
			}
			Index1++;
		}
		Index0++;
	}
	DELETEARRAY(BoxList);
	return true;
}

