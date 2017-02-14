///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains all memory macros.
 *	\file		IceMemoryMacros.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef __ICEMEMORYMACROS_H__
#define __ICEMEMORYMACROS_H__

#undef ZeroMemory
#undef CopyMemory
#undef MoveMemory
#undef FillMemory

	//!	A function to clear a buffer.
	//!	\param	addr		[in] buffer address
	//!	\param	size		[in] buffer length
	//!	\see	FillMemory
	//!	\see	CopyMemory
	//!	\see	MoveMemory
	inline_ void ZeroMemory(void* addr, udword size)					{ memset(addr, 0, size);		}

	//!	A function to fill a buffer with a given byte.
	//!	\param	addr		[in] buffer address
	//!	\param	size		[in] buffer length
	//!	\param	val			[in] the byte value
	//!	\see	ZeroMemory
	//!	\see	CopyMemory
	//!	\see	MoveMemory
	inline_ void FillMemory(void* dest, udword size, ubyte val)			{ memset(dest, val, size);		}

	//!	A function to copy a buffer.
	//!	\param	addr		[in] destination buffer address
	//!	\param	addr		[in] source buffer address
	//!	\param	size		[in] buffer length
	//!	\see	ZeroMemory
	//!	\see	FillMemory
	//!	\see	MoveMemory
	inline_ void CopyMemory(void* dest, const void* src, udword size)	{ memcpy(dest, src, size);		}

	//!	A function to move a buffer.
	//!	\param	addr		[in] destination buffer address
	//!	\param	addr		[in] source buffer address
	//!	\param	size		[in] buffer length
	//!	\see	ZeroMemory
	//!	\see	FillMemory
	//!	\see	CopyMemory
	inline_ void MoveMemory(void* dest, const void* src, udword size)	{ memmove(dest, src, size);		}

	#define SIZEOFOBJECT		sizeof(*this)									//!< Gives the size of current object. Avoid some mistakes (e.g. "sizeof(this)").
	//#define CLEAROBJECT		{ memset(this, 0, SIZEOFOBJECT);	}			//!< Clears current object. Laziness is my business. HANDLE WITH CARE.
	#define DELETESINGLE(x)		if (x) { delete x;				x = null; }		//!< Deletes an instance of a class.
	#define DELETEARRAY(x)		if (x) { delete []x;			x = null; }		//!< Deletes an array.
	#define SAFE_RELEASE(x)		if (x) { (x)->Release();		(x) = null; }	//!< Safe D3D-style release
	#define SAFE_DESTRUCT(x)	if (x) { (x)->SelfDestruct();	(x) = null; }	//!< Safe ICE-style release

#ifdef __ICEERROR_H__
	#define CHECKALLOC(x)		if(!x) return SetIceError("Out of memory.", EC_OUTOFMEMORY);	//!< Standard alloc checking. HANDLE WITH CARE.
#else
	#define CHECKALLOC(x)		if(!x) return false;
#endif

	//! Standard allocation cycle
	#define SAFE_ALLOC(ptr, type, count)	DELETEARRAY(ptr);	ptr = new type[count];	CHECKALLOC(ptr);

#endif // __ICEMEMORYMACROS_H__
