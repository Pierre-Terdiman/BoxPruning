// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__7739E2A3_1E44_11D6_8B0F_0050BAC83302__INCLUDED_)
#define AFX_STDAFX_H__7739E2A3_1E44_11D6_8B0F_0050BAC83302__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "..\Shared\StdAfx.h"

#define USE_HARDCODED_AXES
#define USE_DIRECT_BOUNDS
#define USE_STL

namespace Meshmerizer
{
	#include "IceBoxPruning.h"
}
using namespace Meshmerizer;

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__7739E2A3_1E44_11D6_8B0F_0050BAC83302__INCLUDED_)
