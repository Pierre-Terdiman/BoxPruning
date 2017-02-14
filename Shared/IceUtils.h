///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains misc. useful macros & defines.
 *	\file		IceUtils.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef __ICEUTILS_H__
#define __ICEUTILS_H__

	#define IS_ALIGNED_2(x)		((x&1)==0)
	#define IS_ALIGNED_4(x)		((x&3)==0)
	#define IS_ALIGNED_8(x)		((x&7)==0)

	#define START_RUNONCE		{ static bool __RunOnce__ = false;	if(!__RunOnce__){
	#define END_RUNONCE			__RunOnce__ = true;}}

	//! Reverse all the bits in a 32 bit word (from Steve Baker's Cute Code Collection)
	inline_ void ReverseBits(udword& n)
	{
		n = ((n >>  1) & 0x55555555) | ((n <<  1) & 0xaaaaaaaa) ;
		n = ((n >>  2) & 0x33333333) | ((n <<  2) & 0xcccccccc) ;
		n = ((n >>  4) & 0x0f0f0f0f) | ((n <<  4) & 0xf0f0f0f0) ;
		n = ((n >>  8) & 0x00ff00ff) | ((n <<  8) & 0xff00ff00) ;
		n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000) ;
	}

	//! Count the number of '1' bits in a 32 bit word (from Steve Baker's Cute Code Collection)
	inline_ void CountBits(udword& n)
	{
		n = (n & 0x55555555) + ((n & 0xaaaaaaaa) >> 1);
		n = (n & 0x33333333) + ((n & 0xcccccccc) >> 2);
		n = (n & 0x0f0f0f0f) + ((n & 0xf0f0f0f0) >> 4);
		n = (n & 0x00ff00ff) + ((n & 0xff00ff00) >> 8);
		n = (n & 0x0000ffff) + ((n & 0xffff0000) >> 16);
	}

	//! Test to see if a number is an exact power of two (from Steve Baker's Cute Code Collection)
	inline_ bool IsPowerOfTwo(udword n)				{ return ((n&(n-1))==0);					}

	//! Zero the least significant '1' bit in a word. (from Steve Baker's Cute Code Collection)
	inline_ void ZeroLeastSetBit(udword& n)			{ n&=(n-1);									}

	//! Set the least significant N bits in a word. (from Steve Baker's Cute Code Collection)
	inline_ void SetLeastNBits(udword& x, udword n)	{ x|=~(~0<<n);								}

	//! Classic XOR swap (from Steve Baker's Cute Code Collection)
	inline_ void Swap(udword& x, udword& y)			{ x ^= y; y ^= x; x ^= y;					}

	//! Little/Big endian (from Steve Baker's Cute Code Collection)
	inline_ char LittleEndian()						{ int i = 1; return *((char*)&i);			}

	//!< Alternative abs function
	inline_ udword abs_(sdword x)					{ sdword y= x >> 31; return (x^y)-y;		}

	/*
	"Just call it repeatedly with various input values and always with the same variable as "memory".
	The sharpness determines the degree of filtering, where 0 completely filters out the input, and 1
	does no filtering at all.

	I seem to recall from college that this is called an IIR (Infinite Impulse Response) filter. As opposed
	to the more typical FIR (Finite Impulse Response).

	Also, I'd say that you can make more intelligent and interesting filters than this, for example filters
	that remove wrong responses from the mouse because it's being moved too fast. You'd want such a filter
	to be applied before this one, of course."

	(JCAB on Flipcode)
	*/
	inline_ float FeedbackFilter(float val, float& memory, float sharpness)
	{
		return memory = val * sharpness + memory * (1.0f - sharpness);
	}

	// Generic functions
	template<class Type> inline_ void TSwap(Type &a, Type &b)								{ const Type c = a; a = b; b = c;			}
	template<class Type> inline_ Type TClamp(const Type &x, const Type &lo, const Type &hi)	{ return ((x<lo) ? lo : (x>hi) ? hi : x);	}

	// Prevent nasty user-manipulations (strategy borrowed from Charles Bloom)
//	#define PREVENT_COPY(curclass)	void operator = (const curclass& object)	{	ASSERT(!"Bad use of operator =");	}
	// ... actually this is better !
	#define PREVENT_COPY(curclass)	private: curclass(const curclass& object);	curclass& operator=(const curclass& object);

	//! TO BE DOCUMENTED
	#define OFFSET_OF(Class, Member)	(size_t)&(((Class*)0)->Member)
	//! TO BE DOCUMENTED
//	#define ARRAYSIZE(p)				(sizeof(p)/sizeof(p[0]))

	ICECORE_API udword Alignment(udword address);

#endif // __ICEUTILS_H__
