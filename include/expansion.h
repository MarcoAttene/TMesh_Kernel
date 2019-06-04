/****************************************************************************
* expansion.h                                                               *
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

#ifndef EXPANSION_H
#define EXPANSION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

#ifndef FABS
#define FABS(a)  fabs(a)
#endif
#define TO_NONNEGATIVE(d) ((d) = (FABS(d)))


// The following macros are fast implementations of basic expansion arithmetic due
// to Dekker, Knuth, Priest, Shewchuk, and others.

// See Y. Hida, X. S. Li,  D. H. Bailey "Algorithms for Quad-Double Precision Floating Point Arithmetic"

// Sums
#define Quick_Two_Sum(a, b, x, y) x = a + b; y = b - (x - a)
#define Two_Sum(a, b, x, y) x = a + b; _bv = x - a; y = (a - (x - _bv)) + (b - _bv)
#define Two_One_Sum(a1, a0, b, x2, x1, x0) Two_Sum(a0, b , _i, x0); Two_Sum(a1, _i, x2, x1)

// Differences
#define Two_Diff(a, b, x, y) x = a - b; _bv = a - x; y = (a - (x + _bv)) + (_bv - b)
#define Two_One_Diff(a1, a0, b, x2, x1, x0) Two_Diff(a0, b , _i, x0); Two_Sum( a1, _i, x2, x1)

// Products
#define Split(a, _ah, _al) _c = 1.3421772800000003e+008 * a; _ah = _c - (_c - a); _al = a - _ah
#define Two_Prod_PreSplit(a, b, _bh, _bl, x, y) x = a * b; Split(a, _ah, _al); y = (_al * _bl) - (((x - (_ah * _bh)) - (_al * _bh)) - (_ah * _bl))
#define Two_Product_2Presplit(a, _ah, _al, b, _bh, _bl, x, y) x = a * b; y = (_al * _bl) - (((x - _ah * _bh) - (_al * _bh)) - (_ah * _bl))


// An instance of the following must be created to access functions for expansion arithmetic
class expansionObject
{
	// Temporary vars used in low-level arithmetic
	double _bv, _c, _ah, _al, _bh, _bl, _ch, _cl, _i, _j, _k, _l, _m, _n, _0, _1, _2, _u3;

public:
	expansionObject() {}

	inline void two_Sum(const double& a, const double&b, double& x, double& y) { Two_Sum(a, b, x, y); }
	inline void two_Diff(const double& a, const double&b, double& x, double& y) { Two_Diff(a, b, x, y); }

	// [x,y] = [a]*[b]		 Multiplies two expansions [a] and [b] of length one
	inline void Two_Prod(const double& a, const double&b, double& x, double& y) 
	{
		_u3 = a * b; 
		Split(a, _ah, _al); Split(b, _bh, _bl); 
		y = ((_ah*_bh - _u3) + _ah*_bl + _al*_bh) + _al*_bl;
		x = _u3;
	}
	inline void Two_Prod(const double& a, const double&b, double*xy) { Two_Prod(a, b, xy[1], xy[0]); }

	// [x,y] = [a]^2		Squares an expansion of length one
	inline void Square(const double& a, double& x, double& y)
	{
		x = a * a;
		Split(a, _ah, _al);
		y = (_al * _al) - ((x - (_ah * _ah)) - ((_ah + _ah) * _al));
	}

	// [x3,x2,x1,x0] = [a1,a0]*[b]		Multiplies an expansion [a1,a0] of length two by an expansion [b] of length one
	inline void Two_One_Prod(const double& a1, const double&a0, const double& b, double& x3, double& x2, double& x1, double& x0)
	{
		Split(b, _bh, _bl); 
		Two_Prod_PreSplit(a0, b, _bh, _bl, _i, x0); Two_Prod_PreSplit(a1, b, _bh, _bl, _j, _0); 
		Two_Sum(_i, _0, _k, x1); Quick_Two_Sum(_j, _k, x3, x2);
	}
	inline void Two_One_Prod(const double *a, const double& b, double *x) { Two_One_Prod(a[1], a[0], b, x[3], x[2], x[1], x[0]); }

	// [x3,x2,x1,x0] = [a1,a0]+[b1,b0]		Calculates the sum of two expansions of length two
	inline void Two_Two_Sum(const double& a1, const double&a0, const double& b1, const double& b0, double& x3, double& x2, double& x1, double& x0)
	{
		Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1);
	}

	// [x3,x2,x1,x0] = [a1,a0]-[b1,b0]		Calculates the difference between two expansions of length two
	inline void Two_Two_Diff(const double& a1, const double&a0, const double& b1, const double& b0, double& x3, double& x2, double& x1, double& x0) 
	{
		Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, _u3, x2, x1); x3 = _u3;
	}
	inline void Two_Two_Diff(const double *a, const double *b, double *x) { Two_Two_Diff(a[1], a[0], b[1], b[0], x[3], x[2], x[1], x[0]); }

	// Calculates the second component 'y' of the expansion [x,y] = [a]-[b] when 'x' is known
	inline void Two_Diff_Back(const double& a, const double&b, double& x, double& y) { _bv = a - x; y = (a - (x + _bv)) + (_bv - b); }
	inline void Two_Diff_Back(const double& a, const double&b, double *xy) { Two_Diff_Back(a, b, xy[1], xy[0]); }

	// [h] = [a1,a0]^2		Squares an expansion of length 2
	// 'h' must be allocated by the caller with 6 components.
	void Two_Square(const double& a1, const double& a0, double *x);

	// [h7,h6,...,h0] = [a1,a0]*[b1,b0]		Calculates the product of two expansions of length two.
	// 'h' must be allocated by the caller with eight components.
	void Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double *h);

	// [x] = 2*[e]		Fast multiplication by two
	// 'h' must be allocated by the caller with at least elen components.
	inline int Gen_Doubleval(const int elen, const double *e, double *h)
	{
		for (int i = 0; i < elen; i++) h[i] = e[i]*2;
		return elen;
	}

	// [e] <- 2*[e]		Inplace inversion
	inline void Gen_Invert(const int elen, double *e) {	for (int i = 0; i < elen; i++) e[i] = -e[i]; }

	// [h] = [e] + [f]		Sums two expansions and returns number of components of result
	// 'h' must be allocated by the caller with at least elen+flen components.
	int Gen_Sum(const int elen, const double *e, const int flen, const double *f, double *h);

	// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
	int Gen_Sum_With_Alloc(const int elen, const double *e, const int flen, const double *f, double **h)
	{
		*h = (double *)malloc((elen + flen) * sizeof(double));
		return Gen_Sum(elen, e, flen, f, *h);
	}

	// [h] = [e] * b		Multiplies an expansion by a scalar
	// 'h' must be allocated by the caller with at least elen*2 components.
	int Gen_Scale(const int elen, const double *e, const double& b, double *h);

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least 2*alen*blen components.
	int Sub_product(const int alen, const double *a, const int blen, const double *b, double *h);

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least MAX(2*alen*blen, 8) components.
	int Gen_Product(const int alen, const double *a, const int blen, const double *b, double *h);

	// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
	int Gen_Product_With_Alloc(const int alen, const double *a, const int blen, const double *b, double **h)
	{
		int h_len = alen * blen * 2;
		if (h_len < 8) h_len = 8;
		*h = (double *)malloc(h_len * sizeof(double));
		return Gen_Product(alen, a, blen, b, *h);
	}

	// Approximates the expansion to a double
	static double To_Double(const int elen, const double *e);
};

#endif // EXPANSION_H
