#include "stdafx.h"

#ifdef REMOVED

void RunTest();

// Play with the settings!
static const udword NbBoxes = 5000;
static const float BoxSize = 20.0f;
static const float Spread = 2000.0f;

int main(int argc, char* argv[])
{
	if(1)
	{
		RunTest();
		return 0;
	}

	// 1) Create random boxes
	AABB* Boxes = new AABB[NbBoxes];
	const AABB** List = new const AABB*[NbBoxes];
	for(udword i=0;i<NbBoxes;i++)
	{
		// Get a random center in a cube
		Point Center;
		Center.x = (UnitRandomFloat()-0.5f)*Spread;
		Center.y = (UnitRandomFloat()-0.5f)*Spread;
		Center.z = (UnitRandomFloat()-0.5f)*Spread;

		// Get random extents
		Point Extents;
		Extents.x = UnitRandomFloat()*BoxSize;
		Extents.y = UnitRandomFloat()*BoxSize;
		Extents.z = UnitRandomFloat()*BoxSize;

		// Setup random box
		Boxes[i].mMin.x = Center.x - Extents.x;
		Boxes[i].mMin.y = Center.y - Extents.y;
		Boxes[i].mMin.z = Center.z - Extents.z;
		Boxes[i].mMax.x = Center.x + Extents.x;
		Boxes[i].mMax.y = Center.y + Extents.y;
		Boxes[i].mMax.z = Center.z + Extents.z;
		List[i] = &Boxes[i];
	}

	const udword NB = 16;

	// 2) Do queries
	udword Time;
	Axes axes;
	axes.Axis0 = 0;
	axes.Axis1 = 2;
	axes.Axis2 = 1;
	// 2-1) Brute-force
	Container Pairs;
	udword MinTime = 0xffffffff;

//	StartProfile(Time);
	for(udword i=0;i<NB;i++)
	{
		Pairs.Reset();
		StartProfile(Time);
			BruteForceCompleteBoxTest(NbBoxes, List, Pairs);
		EndProfile(Time);
		if(Time<MinTime)
			MinTime = Time;
	}
//	EndProfile(Time);
//	printf("Brute force: found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, Time/(1024*NB));
	printf("Brute force: found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);

	// 2-1) Sweep&prune
//	Pairs.Reset();
	MinTime = 0xffffffff;
//	StartProfile(Time);
	for(udword i=0;i<NB;i++)
	{
		Pairs.Reset();
		StartProfile(Time);
			CompleteBoxPruning(NbBoxes, List, Pairs, axes);
		EndProfile(Time);
		if(Time<MinTime)
			MinTime = Time;
	}
//	EndProfile(Time);
	printf("Sweep&Prune: found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, MinTime/1024);
//	printf("Sweep&Prune: found %d intersections in %d K-cycles.\n", Pairs.GetNbEntries()>>1, Time/(1024*NB));

	// 3) Free & exit
	DELETEARRAY(Boxes);

	return 0;
}

#endif