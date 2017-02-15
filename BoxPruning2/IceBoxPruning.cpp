///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for box pruning.
 *	\file		IceBoxPruning.cpp
 *	\author		Pierre Terdiman
 *	\date		January, 29, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	You could use a complex sweep-and-prune as implemented in I-Collide.
	You could use a complex hashing scheme as implemented in V-Clip or recently in ODE it seems.
	You could use a "Recursive Dimensional Clustering" algorithm as implemented in GPG2.

	Or you could use this.
	Faster ? I don't know. Probably not. It would be a shame. But who knows ?
	Easier ? Definitely. Enjoy the sheer simplicity.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Precompiled Header
#include "Stdafx.h"

using namespace Meshmerizer;

// InsertionSort has better coherence, RadixSort is better for one-shot queries.
#define PRUNING_SORTER	RadixSort
//#define PRUNING_SORTER	InsertionSort

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
bool Meshmerizer::BipartiteBoxPruning(udword nb0, const AABB** list0, udword nb1, const AABB** list1, Container& pairs, const Axes& axes)
{
	// Checkings
	if(!nb0 || !list0 || !nb1 || !list1)	return false;

	// Catch axes
	udword Axis0 = axes.Axis0;
	udword Axis1 = axes.Axis1;
	udword Axis2 = axes.Axis2;

	// Allocate some temporary data
	float* MinPosList0 = new float[nb0];
	float* MinPosList1 = new float[nb1];

	// 1) Build main lists using the primary axis
	for(udword i=0;i<nb0;i++)	MinPosList0[i] = list0[i]->GetMin(Axis0);
	for(udword i=0;i<nb1;i++)	MinPosList1[i] = list1[i]->GetMin(Axis0);

	// 2) Sort the lists
	static PRUNING_SORTER RS0, RS1;	// Static for coherence.
	udword* Sorted0 = RS0.Sort(MinPosList0, nb0).GetRanks();
	udword* Sorted1 = RS1.Sort(MinPosList1, nb1).GetRanks();

	// 3) Prune the lists
	udword Index0, Index1;

	const udword* const LastSorted0 = &Sorted0[nb0];
	const udword* const LastSorted1 = &Sorted1[nb1];
	const udword* RunningAddress0 = Sorted0;
	const udword* RunningAddress1 = Sorted1;

	while(RunningAddress1<LastSorted1 && Sorted0<LastSorted0)
	{
		Index0 = *Sorted0++;

		while(RunningAddress1<LastSorted1 && MinPosList1[*RunningAddress1]<MinPosList0[Index0])	RunningAddress1++;

		const udword* RunningAddress2_1 = RunningAddress1;

		while(RunningAddress2_1<LastSorted1 && MinPosList1[Index1 = *RunningAddress2_1++]<=list0[Index0]->GetMax(Axis0))
		{
			if(list0[Index0]->Intersect(*list1[Index1], Axis1))
			{
				if(list0[Index0]->Intersect(*list1[Index1], Axis2))
				{
					pairs.Add(Index0).Add(Index1);
				}
			}
		}
	}

	////

	while(RunningAddress0<LastSorted0 && Sorted1<LastSorted1)
	{
		Index0 = *Sorted1++;

		while(RunningAddress0<LastSorted0 && MinPosList0[*RunningAddress0]<=MinPosList1[Index0])	RunningAddress0++;

		const udword* RunningAddress2_0 = RunningAddress0;

		while(RunningAddress2_0<LastSorted0 && MinPosList0[Index1 = *RunningAddress2_0++]<=list1[Index0]->GetMax(Axis0))
		{
			if(list0[Index1]->Intersect(*list1[Index0], Axis1))
			{
				if(list0[Index1]->Intersect(*list1[Index0], Axis2))
				{
					pairs.Add(Index1).Add(Index0);
				}
			}

		}
	}

	DELETEARRAY(MinPosList1);
	DELETEARRAY(MinPosList0);

	return true;
}

//#include <intrin.h>

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
bool Meshmerizer::CompleteBoxPruning(udword nb, const AABB** list, Container& pairs, const Axes& axes)
{
	// Checkings
	if(!nb || !list)
		return false;

	// Catch axes
	udword Axis0 = axes.Axis0;
	udword Axis1 = axes.Axis1;
	udword Axis2 = axes.Axis2;

	// Allocate some temporary data
	float* PosList = new float[nb];

	// 1) Build main list using the primary axis
	for(udword i=0;i<nb;i++)
		PosList[i] = list[i]->GetMin(Axis0);

	// 2) Sort the list
	static PRUNING_SORTER RS;	// Static for coherence
	const udword* Sorted = RS.Sort(PosList, nb).GetRanks();

//CPUID
	// 3) Prune the list
	const udword* const LastSorted = &Sorted[nb];
	const udword* RunningAddress = Sorted;
	udword Index0, Index1;
	while(RunningAddress<LastSorted && Sorted<LastSorted)
	{
		Index0 = *Sorted++;

		while(RunningAddress<LastSorted && PosList[*RunningAddress++]<PosList[Index0]);

		const udword* RunningAddress2 = RunningAddress;

//		NOPS
		while(RunningAddress2<LastSorted && PosList[Index1 = *RunningAddress2++]<=list[Index0]->GetMax(Axis0))
		{
//		NOPS

			if(Index0!=Index1)
			{
				if(list[Index0]->Intersect(*list[Index1], Axis1))
				{
					if(list[Index0]->Intersect(*list[Index1], Axis2))
					{
						pairs.Add(Index0).Add(Index1);
					}
				}
			}
		}
	}
//CPUID
//	unsigned long long time = __rdtsc();

	DELETEARRAY(PosList);

//	time = __rdtsc() - time;
//	printf("time: %d\n", time/1024);

	return true;
}
