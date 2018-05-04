#include "Stdafx.h"
#include <conio.h>

static const bool gPresortBounds = false;
static const bool gTestCompleteBoxPruning = true;
static const bool gTestBipartiteBoxPruning = false;

#ifdef USE_STL
static bool BruteForceCompleteBoxTest(udword nb, const AABB** list, std::vector<udword>& pairs)
{
	// Checkings
	if(!nb || !list)
		return false;

	// Brute-force n(n-1)/2 overlap tests
	for(udword i=0;i<nb;i++)
	{
		for(udword j=i+1;j<nb;j++)
		{
			if(list[i]->Intersect(*list[j]))
			{
				pairs.push_back(i);
				pairs.push_back(j);
			}
		}
	}
	return true;
}
#endif

static void RunPerformanceTest()
{
	const udword NbBoxes = 10000;

	// Create random boxes
	AABB* Boxes = new AABB[NbBoxes];
	const AABB** List = new const AABB*[NbBoxes];

	srand(42);
    const udword BoxSize = 128;
    for(udword i=0;i<NbBoxes;i++)
    {
        const float x = float(rand() & 4095) - 2048.0f;// + UnitRandomFloat();
        const float y = float(rand() & 4095) - 2048.0f;// + UnitRandomFloat();
        const float z = float(rand() & 4095) - 2048.0f;// + UnitRandomFloat();
        float ex = float(rand() & (BoxSize-1));
        float ey = float(rand() & (BoxSize-1));
        float ez = float(rand() & (BoxSize-1));
		// Avoid zero-area boxes for testing, as they're not handled consistently by all implementations
/*		if(ex==0.0f)
			ex = 0.01f;
		if(ey==0.0f)
			ey = 0.01f;
		if(ez==0.0f)
			ez = 0.01f;*/
		Boxes[i].mMin.x = x - ex;
		Boxes[i].mMin.y = y - ey;
		Boxes[i].mMin.z = z - ez;
		Boxes[i].mMax.x = x + ex;
		Boxes[i].mMax.y = y + ey;
		Boxes[i].mMax.z = z + ez;

		List[i] = &Boxes[i];
    }

	if(gPresortBounds)
	{
		float* PosList = new float[NbBoxes];
		for(udword i=0;i<NbBoxes;i++)
			PosList[i] = Boxes[i].mMin.x;

		RadixSort RS;
		const udword* Sorted = RS.Sort(PosList, NbBoxes).GetRanks();
		
		AABB* Copy = new AABB[NbBoxes];
		memcpy(Copy, Boxes, sizeof(AABB)*NbBoxes);

		for(udword i=0;i<NbBoxes;i++)
			Boxes[i] = Copy[Sorted[i]];

		DELETEARRAY(Copy);
		DELETEARRAY(PosList);
	}


#ifndef USE_HARDCODED_AXES
	Axes axes;
	axes.Axis0 = 0;
	axes.Axis1 = 2;
	axes.Axis2 = 1;
#endif

#ifdef USE_STL
	std::vector<udword> Pairs;
#else
	Container Pairs;
//	Pairs.SetSize(25000*2);
//	Pairs.Reset();
#endif
	const udword NB = 160;

	// Test "complete" pruning
	if(gTestCompleteBoxPruning)
	{
		udword MinTime = 0xffffffff;
		udword Time;

		if(0)
		{
			// Brute-force
			for(udword i=0;i<NB;i++)
			{
#ifdef USE_STL
				Pairs.clear();
#else
				Pairs.Reset();
#endif
				StartProfile(Time);
					BruteForceCompleteBoxTest(NbBoxes, List, Pairs);
				EndProfile(Time);
				if(Time<MinTime)
					MinTime = Time;
			}
#ifdef USE_STL
			printf("Complete test (brute force): found %d intersections in %d K-cycles.\n", Pairs.size()>>1, MinTime/1024);
#else
			printf("Complete test (brute force): found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);
#endif
		}

		// Optimized
		MinTime = 0xffffffff;
		for(udword i=0;i<NB;i++)
		{
#ifdef USE_STL
				Pairs.clear();
//				Pairs.shrink_to_fit();
//				printf("capa: %d\n", Pairs.capacity());
#else
				Pairs.Reset();
#endif

			StartProfile(Time);
#ifdef USE_HARDCODED_AXES
	#ifdef USE_DIRECT_BOUNDS
//				CompleteBoxPruningSTL(NbBoxes, Boxes, STLPairs);
				CompleteBoxPruning(NbBoxes, Boxes, Pairs);
	#else
				CompleteBoxPruning(NbBoxes, List, Pairs);
	#endif
#else
				CompleteBoxPruning(NbBoxes, List, Pairs, axes);
#endif
			EndProfile(Time);
			if(Time<MinTime)
				MinTime = Time;
			printf(" %d K-cycles.\n", Time/1024);
		}
#ifdef USE_STL
		printf("Complete test (box pruning): found %d intersections in %d K-cycles.\n", Pairs.size()>>1, MinTime/1024);
#else
		printf("Complete test (box pruning): found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);
#endif
	}

#ifndef USE_STL
	// Test "bipartite" pruning
	if(gTestBipartiteBoxPruning)
	{
		const udword NbBoxes0 = NbBoxes/2;
		const udword NbBoxes1 = NbBoxes - NbBoxes0;
		const AABB** List1 = List + NbBoxes0;
		const AABB* Boxes1 = Boxes + NbBoxes0;

		udword MinTime = 0xffffffff;
		udword Time;

		if(0)
		{
			// Brute-force
			for(udword i=0;i<NB;i++)
			{
				Pairs.Reset();
				StartProfile(Time);
					BruteForceBipartiteBoxTest(NbBoxes0, List, NbBoxes1, List1, Pairs);
				EndProfile(Time);
				if(Time<MinTime)
					MinTime = Time;
			}
			printf("Bipartite test (brute force): found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);
		}

		// Optimized
		MinTime = 0xffffffff;
		for(udword i=0;i<NB;i++)
		{
			Pairs.Reset();
			StartProfile(Time);
#ifdef USE_HARDCODED_AXES
	#ifdef USE_DIRECT_BOUNDS
				BipartiteBoxPruning(NbBoxes0, Boxes, NbBoxes1, Boxes1, Pairs);
	#else
				BipartiteBoxPruning(NbBoxes0, List, NbBoxes1, List1, Pairs);
	#endif
#else
				BipartiteBoxPruning(NbBoxes0, List, NbBoxes1, List1, Pairs, axes);
#endif
			EndProfile(Time);
			if(Time<MinTime)
				MinTime = Time;
			printf(" %d K-cycles.\n", Time/1024);
		}
		printf("Bipartite test (box pruning): found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);
	}
#endif

	DELETEARRAY(List);
	DELETEARRAY(Boxes);
}

#ifndef USE_STL
//#define VERBOSE
static void TestCompleteBoxPruning(	udword NbBoxes, const AABB** List, const AABB* Boxes, 
									Container& Pairs0, Container& Pairs1,
									RadixSort& RS0, RadixSort& RS1,
									Container& Keys0, Container& Keys1,
									udword TestIndex)
{
#ifndef USE_HARDCODED_AXES
	Axes axes;
	axes.Axis0 = 0;
	axes.Axis1 = 2;
	axes.Axis2 = 1;
#endif

	Pairs0.Reset();
	BruteForceCompleteBoxTest(NbBoxes, List, Pairs0);

	Pairs1.Reset();
#ifdef USE_HARDCODED_AXES
#ifdef USE_DIRECT_BOUNDS
//	CompleteBoxPruningSTL(NbBoxes, Boxes, STLPairs);
	CompleteBoxPruning(NbBoxes, Boxes, Pairs1);
#else
	CompleteBoxPruning(NbBoxes, List, Pairs1);
#endif
#else
	CompleteBoxPruning(NbBoxes, List, Pairs1, axes);
#endif

	if(Pairs0.GetNbEntries() != Pairs1.GetNbEntries())
	{
		printf("ERROR: different nb pairs!\n");

		const udword NbPairs0 = Pairs0.GetNbEntries()>>1;
#ifdef VERBOSE
		printf(" NbPairs0: %d\n", NbPairs0);
#endif
		const udword NbPairs1 = Pairs1.GetNbEntries()>>1;
#ifdef VERBOSE
		printf(" NbPairs1: %d\n", NbPairs1);
#endif
		Pair* Entries0 = (Pair*)Pairs0.GetEntries();
		Pair* Entries1 = (Pair*)Pairs1.GetEntries();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs0;i++)
		{
			udword id0 = Entries0[i].id0;
			udword id1 = Entries0[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries0[i].id0 = id0;
			Entries0[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS0.Sort(Keys1.GetEntries(), NbPairs0, false);
		RS0.Sort(Keys0.GetEntries(), NbPairs0, false);
		const udword* Sorted0 = RS0.GetRanks();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs1;i++)
		{
			udword id0 = Entries1[i].id0;
			udword id1 = Entries1[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries1[i].id0 = id0;
			Entries1[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS1.Sort(Keys1.GetEntries(), NbPairs1, false);
		RS1.Sort(Keys0.GetEntries(), NbPairs1, false);
		const udword* Sorted1 = RS1.GetRanks();

		const udword NbPairs = TMin(NbPairs0, NbPairs1);
		for(udword i=0;i<NbPairs;i++)
		{
			const udword SortedIndex0 = Sorted0[i];
			const udword SortedIndex1 = Sorted1[i];
			if(
				(Entries0[SortedIndex0].id0!=Entries1[SortedIndex1].id0)
				||
				(Entries0[SortedIndex0].id1!=Entries1[SortedIndex1].id1))
			{
#ifdef VERBOSE
				printf("Pair %d: ERROR!\n", i);
				printf(" Entries0: %d | %d\n", Entries0[SortedIndex0].id0, Entries0[SortedIndex0].id1);
				printf(" Entries1: %d | %d\n", Entries1[SortedIndex1].id0, Entries1[SortedIndex1].id1);

				for(udword j=0;j<NbPairs1;j++)
				{
					if(Entries1[j].id0 == Entries0[SortedIndex0].id0
					&&	Entries1[j].id1 == Entries0[SortedIndex0].id1)
						printf("FOUND!\n");
				}

				const AABB& Box0 = Boxes[Entries0[SortedIndex0].id0];
				const AABB& Box1 = Boxes[Entries0[SortedIndex0].id1];
				printf("  %f %f\n", Box0.mMin.x, Box0.mMax.x);
				printf("  %f %f\n", Box0.mMin.y, Box0.mMax.y);
				printf("  %f %f\n", Box0.mMin.z, Box0.mMax.z);
				printf("  %f %f\n", Box1.mMin.x, Box1.mMax.x);
				printf("  %f %f\n", Box1.mMin.y, Box1.mMax.y);
				printf("  %f %f\n", Box1.mMin.z, Box1.mMax.z);
				printf("Intersect: %d\n", Box0.Intersect(Box1));

				AABB Tmp[2];
				Tmp[0] = Box0;
				Tmp[1] = Box1;
				Container TestPairs;
				CompleteBoxPruning(2, Tmp, TestPairs);
				printf("Intersect: %d\n", TestPairs.GetNbEntries());
#endif
				exit(0);
			}
		}


		printf("\n\nTest index: %d\n", TestIndex);
		exit(0);
	}
	else
	{
		const udword NbPairs = Pairs0.GetNbEntries()>>1;

		Pair* Entries0 = (Pair*)Pairs0.GetEntries();
		Pair* Entries1 = (Pair*)Pairs1.GetEntries();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs;i++)
		{
			udword id0 = Entries0[i].id0;
			udword id1 = Entries0[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries0[i].id0 = id0;
			Entries0[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS0.Sort(Keys1.GetEntries(), NbPairs, false);
		RS0.Sort(Keys0.GetEntries(), NbPairs, false);
		const udword* Sorted0 = RS0.GetRanks();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs;i++)
		{
			udword id0 = Entries1[i].id0;
			udword id1 = Entries1[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries1[i].id0 = id0;
			Entries1[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS1.Sort(Keys1.GetEntries(), NbPairs, false);
		RS1.Sort(Keys0.GetEntries(), NbPairs, false);
		const udword* Sorted1 = RS1.GetRanks();

		for(udword i=0;i<NbPairs;i++)
		{
			const udword SortedIndex0 = Sorted0[i];
			const udword SortedIndex1 = Sorted1[i];
			if(
				(Entries0[SortedIndex0].id0!=Entries1[SortedIndex1].id0)
				||
				(Entries0[SortedIndex0].id1!=Entries1[SortedIndex1].id1))
			{
#ifdef VERBOSE
				printf("Pair %d: ERROR!\n", i);
				printf(" Entries0: %d | %d\n", Entries0[SortedIndex0].id0, Entries0[SortedIndex0].id1);
				printf(" Entries1: %d | %d\n", Entries1[SortedIndex1].id0, Entries1[SortedIndex1].id1);

				for(udword j=0;j<NbPairs;j++)
				{
					if(Entries1[j].id0 == Entries0[SortedIndex0].id0
					&&	Entries1[j].id1 == Entries0[SortedIndex0].id1)
						printf("FOUND!\n");
				}

				const AABB& Box0 = Boxes[Entries0[SortedIndex0].id0];
				const AABB& Box1 = Boxes[Entries0[SortedIndex0].id1];
				printf("  %f %f\n", Box0.mMin.x, Box0.mMax.x);
				printf("  %f %f\n", Box0.mMin.y, Box0.mMax.y);
				printf("  %f %f\n", Box0.mMin.z, Box0.mMax.z);
				printf("  %f %f\n", Box1.mMin.x, Box1.mMax.x);
				printf("  %f %f\n", Box1.mMin.y, Box1.mMax.y);
				printf("  %f %f\n", Box1.mMin.z, Box1.mMax.z);
				printf("Intersect: %d\n", Box0.Intersect(Box1));

				AABB Tmp[2];
				Tmp[0] = Box0;
				Tmp[1] = Box1;
				Container TestPairs;
				CompleteBoxPruning(2, Tmp, TestPairs);
				printf("Intersect: %d\n", TestPairs.GetNbEntries());
#endif
				printf("\n\nTest index: %d\n", TestIndex);
				exit(0);
			}
		}
	}
//	printf("%d\n", Pairs0.GetNbEntries());
}

static void TestBipartiteBoxPruning(udword NbBoxes, const AABB** List, const AABB* Boxes, 
									Container& Pairs0, Container& Pairs1,
									RadixSort& RS0, RadixSort& RS1,
									Container& Keys0, Container& Keys1,
									udword TestIndex)
{
	const udword NbBoxes0 = NbBoxes/2;
	const udword NbBoxes1 = NbBoxes - NbBoxes0;
	const AABB** List1 = List + NbBoxes0;
	const AABB* Boxes1 = Boxes + NbBoxes0;

#ifndef USE_HARDCODED_AXES
	Axes axes;
	axes.Axis0 = 0;
	axes.Axis1 = 2;
	axes.Axis2 = 1;
#endif

	Pairs0.Reset();
	BruteForceBipartiteBoxTest(NbBoxes0, List, NbBoxes1, List1, Pairs0);

	Pairs1.Reset();
#ifdef USE_HARDCODED_AXES
#ifdef USE_DIRECT_BOUNDS
	BipartiteBoxPruning(NbBoxes0, Boxes, NbBoxes1, Boxes1, Pairs1);
#else
#error fixme
	CompleteBoxPruning(NbBoxes, List, Pairs1);
#endif
#else
#error fixme
	CompleteBoxPruning(NbBoxes, List, Pairs1, axes);
#endif

	if(Pairs0.GetNbEntries() != Pairs1.GetNbEntries())
	{
		printf("ERROR: different nb pairs!\n");

#ifdef FIX_ME
		const udword NbPairs0 = Pairs0.GetNbEntries()>>1;
#ifdef VERBOSE
		printf(" NbPairs0: %d\n", NbPairs0);
#endif
		const udword NbPairs1 = Pairs1.GetNbEntries()>>1;
#ifdef VERBOSE
		printf(" NbPairs1: %d\n", NbPairs1);
#endif
		Pair* Entries0 = (Pair*)Pairs0.GetEntries();
		Pair* Entries1 = (Pair*)Pairs1.GetEntries();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs0;i++)
		{
			udword id0 = Entries0[i].id0;
			udword id1 = Entries0[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries0[i].id0 = id0;
			Entries0[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS0.Sort(Keys1.GetEntries(), NbPairs0, false);
		RS0.Sort(Keys0.GetEntries(), NbPairs0, false);
		const udword* Sorted0 = RS0.GetRanks();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs1;i++)
		{
			udword id0 = Entries1[i].id0;
			udword id1 = Entries1[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries1[i].id0 = id0;
			Entries1[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS1.Sort(Keys1.GetEntries(), NbPairs1, false);
		RS1.Sort(Keys0.GetEntries(), NbPairs1, false);
		const udword* Sorted1 = RS1.GetRanks();

		const udword NbPairs = TMin(NbPairs0, NbPairs1);
		for(udword i=0;i<NbPairs;i++)
		{
			const udword SortedIndex0 = Sorted0[i];
			const udword SortedIndex1 = Sorted1[i];
			if(
				(Entries0[SortedIndex0].id0!=Entries1[SortedIndex1].id0)
				||
				(Entries0[SortedIndex0].id1!=Entries1[SortedIndex1].id1))
			{
#ifdef VERBOSE
				printf("Pair %d: ERROR!\n", i);
				printf(" Entries0: %d | %d\n", Entries0[SortedIndex0].id0, Entries0[SortedIndex0].id1);
				printf(" Entries1: %d | %d\n", Entries1[SortedIndex1].id0, Entries1[SortedIndex1].id1);

				for(udword j=0;j<NbPairs1;j++)
				{
					if(Entries1[j].id0 == Entries0[SortedIndex0].id0
					&&	Entries1[j].id1 == Entries0[SortedIndex0].id1)
						printf("FOUND!\n");
				}

				const AABB& Box0 = Boxes[Entries0[SortedIndex0].id0];
				const AABB& Box1 = Boxes[Entries0[SortedIndex0].id1];
				printf("  %f %f\n", Box0.mMin.x, Box0.mMax.x);
				printf("  %f %f\n", Box0.mMin.y, Box0.mMax.y);
				printf("  %f %f\n", Box0.mMin.z, Box0.mMax.z);
				printf("  %f %f\n", Box1.mMin.x, Box1.mMax.x);
				printf("  %f %f\n", Box1.mMin.y, Box1.mMax.y);
				printf("  %f %f\n", Box1.mMin.z, Box1.mMax.z);
				printf("Intersect: %d\n", Box0.Intersect(Box1));

				AABB Tmp[2];
				Tmp[0] = Box0;
				Tmp[1] = Box1;
				Container TestPairs;
				CompleteBoxPruning(2, Tmp, TestPairs);
				printf("Intersect: %d\n", TestPairs.GetNbEntries());
#endif
				exit(0);
			}
		}
#endif

		printf("\n\nTest index: %d\n", TestIndex);
		exit(0);
	}
	else
	{
		const udword NbPairs = Pairs0.GetNbEntries()>>1;

		Pair* Entries0 = (Pair*)Pairs0.GetEntries();
		Pair* Entries1 = (Pair*)Pairs1.GetEntries();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs;i++)
		{
			udword id0 = Entries0[i].id0;
			udword id1 = Entries0[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries0[i].id0 = id0;
			Entries0[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS0.Sort(Keys1.GetEntries(), NbPairs, false);
		RS0.Sort(Keys0.GetEntries(), NbPairs, false);
		const udword* Sorted0 = RS0.GetRanks();

		Keys0.Reset();
		Keys1.Reset();
		for(udword i=0;i<NbPairs;i++)
		{
			udword id0 = Entries1[i].id0;
			udword id1 = Entries1[i].id1;
			if(id1<id0)
				TSwap(id0, id1);
			Entries1[i].id0 = id0;
			Entries1[i].id1 = id1;

			Keys0.Add(id0);
			Keys1.Add(id1);
		}

		RS1.Sort(Keys1.GetEntries(), NbPairs, false);
		RS1.Sort(Keys0.GetEntries(), NbPairs, false);
		const udword* Sorted1 = RS1.GetRanks();

		for(udword i=0;i<NbPairs;i++)
		{
			const udword SortedIndex0 = Sorted0[i];
			const udword SortedIndex1 = Sorted1[i];
			if(
				(Entries0[SortedIndex0].id0!=Entries1[SortedIndex1].id0)
				||
				(Entries0[SortedIndex0].id1!=Entries1[SortedIndex1].id1))
			{
#ifdef VERBOSE
				printf("Pair %d: ERROR!\n", i);
				printf(" Entries0: %d | %d\n", Entries0[SortedIndex0].id0, Entries0[SortedIndex0].id1);
				printf(" Entries1: %d | %d\n", Entries1[SortedIndex1].id0, Entries1[SortedIndex1].id1);

				for(udword j=0;j<NbPairs;j++)
				{
					if(Entries1[j].id0 == Entries0[SortedIndex0].id0
					&&	Entries1[j].id1 == Entries0[SortedIndex0].id1)
						printf("FOUND!\n");
				}

				const AABB& Box0 = Boxes[Entries0[SortedIndex0].id0];
				const AABB& Box1 = Boxes[Entries0[SortedIndex0].id1];
				printf("  %f %f\n", Box0.mMin.x, Box0.mMax.x);
				printf("  %f %f\n", Box0.mMin.y, Box0.mMax.y);
				printf("  %f %f\n", Box0.mMin.z, Box0.mMax.z);
				printf("  %f %f\n", Box1.mMin.x, Box1.mMax.x);
				printf("  %f %f\n", Box1.mMin.y, Box1.mMax.y);
				printf("  %f %f\n", Box1.mMin.z, Box1.mMax.z);
				printf("Intersect: %d\n", Box0.Intersect(Box1));

				AABB Tmp[2];
				Tmp[0] = Box0;
				Tmp[1] = Box1;
				Container TestPairs;
				CompleteBoxPruning(2, Tmp, TestPairs);
				printf("Intersect: %d\n", TestPairs.GetNbEntries());
#endif
				printf("\n\nTest index: %d\n", TestIndex);
				exit(0);
			}
		}
	}
//	printf("%d\n", Pairs0.GetNbEntries());
}

static void RunValidityTest()
{
	srand(42);
	Container Pairs0;
	Container Pairs1;

	RadixSort RS0;
	RadixSort RS1;

	Container Keys0;
	Container Keys1;

	const udword NbTests = 1000;
//	const udword NbTests = 1000000000;
	for(udword TestIndex=0;TestIndex<NbTests;TestIndex++)
	{
#ifdef VERBOSE
		printf("%d\n", TestIndex);
		const udword NbBoxes = rand() & 1023;
		printf("NbBoxes: %d\n", NbBoxes);
//		if(TestIndex==1)
//			NbBoxes = 12;
#else
		printf(".");
		const udword NbBoxes = rand() & 1023;
#endif

		AABB* Boxes = new AABB[NbBoxes];
		const AABB** List = new const AABB*[NbBoxes];

		const udword BoxSize = rand() & 255;
		const udword Range = rand() & 4095;
		for(udword i=0;i<NbBoxes;i++)
		{
			const float x = float(rand() & Range) - float(Range/2);
			const float y = float(rand() & Range) - float(Range/2);
			const float z = float(rand() & Range) - float(Range/2);
			float ex = float(rand() & (BoxSize-1));
			float ey = float(rand() & (BoxSize-1)) + UnitRandomFloat();
			float ez = float(rand() & (BoxSize-1)) + UnitRandomFloat();
			// Avoid zero-area boxes for testing, as they're not handled consistently by all implementations
/*			if(ex==0.0f)
				ex = 0.01f;
			if(ey==0.0f)
				ey = 0.01f;
			if(ez==0.0f)
				ez = 0.01f;*/
			Boxes[i].mMin.x = x - ex;
			Boxes[i].mMin.y = y - ey;
			Boxes[i].mMin.z = z - ez;
			Boxes[i].mMax.x = x + ex;
			Boxes[i].mMax.y = y + ey;
			Boxes[i].mMax.z = z + ez;

			List[i] = &Boxes[i];
		}

		if(0)
		{
			float* PosList = new float[NbBoxes];
			for(udword i=0;i<NbBoxes;i++)
				PosList[i] = Boxes[i].mMin.x;

			RadixSort RS;
			const udword* Sorted = RS.Sort(PosList, NbBoxes).GetRanks();
		
			AABB* Copy = new AABB[NbBoxes];
			memcpy(Copy, Boxes, sizeof(AABB)*NbBoxes);

			for(udword i=0;i<NbBoxes;i++)
				Boxes[i] = Copy[Sorted[i]];

			DELETEARRAY(Copy);
			DELETEARRAY(PosList);
		}

		// Test "complete" pruning
		if(gTestCompleteBoxPruning)
		{
//			if(TestIndex==35)
//				int stop=0;
//			if(TestIndex==1)
				TestCompleteBoxPruning(NbBoxes, List, Boxes, Pairs0, Pairs1, RS0, RS1, Keys0, Keys1, TestIndex);
		}

		// Test "bipartite" pruning
		if(gTestBipartiteBoxPruning)
			TestBipartiteBoxPruning(NbBoxes, List, Boxes, Pairs0, Pairs1, RS0, RS1, Keys0, Keys1, TestIndex);

		DELETEARRAY(List);
		DELETEARRAY(Boxes);
	}
	printf("\nFinished.\n");
}
#endif

/*static void RunEdgeCase()
{
	Container Pairs;

	AABB Boxes[2];
	Boxes[0].mMin = Point(-1.0f, -1.0f, -1.0f);
	Boxes[0].mMax = Point(1.0f, 1.0f, 1.0f);
	Boxes[1].mMin = Point(-1.0f, -1.0f, -1.0f);
	Boxes[1].mMax = Point(1.0f, 1.0f, 1.0f);

	CompleteBoxPruning(2, Boxes, Pairs);
}*/

int main(int argc, char* argv[])
{
	RunPerformanceTest();
//	RunValidityTest();
//	RunEdgeCase();

	while(!_kbhit());

	return 0;
}
