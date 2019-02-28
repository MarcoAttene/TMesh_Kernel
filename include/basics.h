/****************************************************************************
* basics.h                                                                  *
* This file is part of the TMesh_Kernel Library                             *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
* IMATI-GE / CNR is Consiglio Nazionale delle Ricerche                      *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Genova (Italy)                                                            *
*                                                                           *
* TMesh_Kernel is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* TMesh_Kernel is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with TMesh_Kernel.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

#ifndef _BASICS_H
#define _BASICS_H

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include "coordinates.h"
#include <stdint.h>

namespace T_MESH
{

#ifdef EXTENSIBLE_TMESH
#define TMESH_VIRTUAL virtual
#else
#define TMESH_VIRTUAL
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DISPMSG_ACTION_SETWIDGET	1
#define DISPMSG_ACTION_PUTNEWLINE	2
#define DISPMSG_ACTION_PUTPROGRESS	3
#define DISPMSG_ACTION_PUTMESSAGE	4
#define DISPMSG_ACTION_ERRORDIALOG	5

typedef unsigned char	UBYTE;
typedef   signed char	 BYTE;
#ifndef _INC_WINDOWS
typedef unsigned short UINT16;
typedef   signed short	INT16;
#endif

#if defined (IS64BITPLATFORM) || defined(_WIN64) || defined(__x86_64__) || defined(__ppc64__)
typedef int64_t j_voidint;
#else
typedef int32_t	 j_voidint;
#endif 


#define FABS(a) (((a)<0)?(-(a)):(a))
#define LOG2(a) (log(a)/log(2))
#define PI2	(M_PI/2.0)

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif


//////// Swaps two pointers. ///////////////////////////////

inline void p_swap(void **a, void **b) { void *t = *a; *a = *b; *b = t; }

/////////////////////////////////////////////////////////////////////////////////////////////


#define TMESH_VERSION	"5.0"
#define TMESH_YEAR		2019


struct thread_initializer
{
	explicit thread_initializer() {}

	//Copy constructor that does the init
	thread_initializer(thread_initializer& _it)
	{
#ifdef USE_HYBRID_KERNEL
		setFPUModeToRoundUP();
#endif
	}
};

extern thread_initializer tmesh_thread_initializer;

class TMesh
{
 public:

 static void (*display_message)(const char *, int);
 static char *app_name;
 static char *app_version;
 static char *app_year;
 static char *app_authors;
 static char *app_url;
 static char *app_maillist;

 static const char *filename; // This might be null. If not, it represents the file we are currently working with.

 static bool quiet;

 static void init(void (*)(const char *, int) = NULL);

 static void info(const char *, ...);
 static void warning(const char *, ...);
 static void error(const char *, ...);
 static void begin_progress();
 static void report_progress(const char *, ...);
 static void end_progress();

 //! When called with a nonzero argument 'ts', launches a chronometer with a timeout of 'ts' seconds.
 //! Later calls without arguments check the chronometer, and if it is over 'ts' the program exits.
 static void exitOnTimeout(clock_t ts = 0);

 //! Appends a line to the file "tmesh.log"
 static void addMessageToLogFile(const char *msg);

 //! Formats a message headed with date/time/filename, appends it to "tmesh.log", and exits with error
 static void logToFileAndExit(const char *msg);

 //! When called without arguments prints the elapsed time from the latest reset.
 static void printElapsedTime(bool reset = false);

 //! Turn and check the use of rational arithmetic
 static void useRationals(bool u);
 static bool isUsingRationals();

 //! Swith to exact arithmetic and return the status before the switch
 static bool useRationals() { bool t = isUsingRationals(); useRationals(true); return t; }

 //! Set the working project's filename
 static void setFilename(const char *fname) { filename = fname; }

 static double tri_orientation(double *pa, double *pb, double *pc);
 static double tet_orientation(double *pa, double *pb, double *pc, double *pd);
 static double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
 static double incircle3D(double *pa, double *pb, double *pc, double *pd);
};

} //namespace T_MESH

#endif //_BASICS_H

