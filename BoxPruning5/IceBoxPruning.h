///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains code for the "box pruning revisited" project.
 *	\file		IceBoxPruning.cpp
 *	\author		Pierre Terdiman
 *	\date		February 2017
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEBOXPRUNING_H
#define ICEBOXPRUNING_H

	// Optimized versions
	FUNCTION MESHMERIZER_API bool CompleteBoxPruning(udword nb, const AABB** list, Container& pairs);
	FUNCTION MESHMERIZER_API bool BipartiteBoxPruning(udword nb0, const AABB** list0, udword nb1, const AABB** list1, Container& pairs);

#endif // ICEBOXPRUNING_H
