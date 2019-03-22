/****************************************************************************
* expansion.cpp                                                             *
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

/*****************************************************************************/
/*                                                                           */
/*  Arithmetic Expansion                                                     */
/*    Code in this file is freely inspired from ideas first presented        */
/*    in Jonathan Richard Shewchuk's paper "Adaptive Precision Floating-     */
/*    Point3c Arithmetic and Fast Robust Geometric Predicates."  Discrete &    */
/*    Computational Geometry (1998).                                         */
/*                                                                           */
/*****************************************************************************/

#include "expansion.h"
#include "float.h"

#pragma optimize("", off)

int expansionObject::Gen_Sum(const int elen, const double *e, const int flen, const double *f, double *h)
{
  double Q, Qn, hh, en = e[0], fn = f[0];
  int e_i, f_i, h_i;

  h_i = e_i = f_i = 0;
  if ((fn > en) == (fn > -en)) { Q = en; en = e[++e_i]; } 
  else { Q = fn; fn = f[++f_i]; }

  if ((e_i < elen) && (f_i < flen)) 
  {
    if ((fn > en) == (fn > -en)) { Quick_Two_Sum(en, Q, Qn, hh); en = e[++e_i]; }
	else { Quick_Two_Sum(fn, Q, Qn, hh); fn = f[++f_i]; }
    Q = Qn;
    if (hh != 0.0) h[h_i++] = hh;
    while ((e_i < elen) && (f_i < flen)) 
	{
      if ((fn > en) == (fn > -en)) { Two_Sum(Q, en, Qn, hh); en = e[++e_i]; }
	  else { Two_Sum(Q, fn, Qn, hh); fn = f[++f_i]; }
      Q = Qn;
      if (hh != 0.0) h[h_i++] = hh;
    }
  }

  while (e_i < elen) 
  {
    Two_Sum(Q, en, Qn, hh);
    en = e[++e_i];
    Q = Qn;
    if (hh != 0.0) h[h_i++] = hh;
  }
  while (f_i < flen) 
  {
    Two_Sum(Q, fn, Qn, hh);
    fn = f[++f_i];
    Q = Qn;
    if (hh != 0.0) h[h_i++] = hh;
  }
  if ((Q != 0.0) || (h_i == 0)) h[h_i++] = Q;
  
  return h_i;
}

int expansionObject::Gen_Scale(const int elen, const double *e, const double& b, double *h)
{
  double Q, sum, hh, pr1, pr0, enow;
  int e_i, h_i;

  Split(b, _bh, _bl);
  Two_Prod_PreSplit(e[0], b, _bh, _bl, Q, hh);
  h_i = 0;
  if (hh != 0) h[h_i++] = hh;
  
  for (e_i = 1; e_i < elen; e_i++) 
  {
    enow = e[e_i];
    Two_Prod_PreSplit(enow, b, _bh, _bl, pr1, pr0);
    Two_Sum(Q, pr0, sum, hh);
    if (hh != 0) h[h_i++] = hh;
    Quick_Two_Sum(pr1, sum, Q, hh);
    if (hh != 0) h[h_i++] = hh;
  }
  if ((Q != 0.0) || (h_i == 0)) h[h_i++] = Q;
  
  return h_i;
}


void expansionObject::Two_Square(const double& a1, const double& a0, double *x)
{
	Square(a0, _j, x[0]);
	_0 = a0 + a0;
	Two_Prod(a1, _0, _k, _1);
	Two_One_Sum(_k, _1, _j, _l, _2, x[1]);
	Square(a1, _j, _1);
	Two_Two_Sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
}


void expansionObject::Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double *h)
{
	Split(a0, _ah, _al);
	Split(b0, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b0, _bh, _bl, _i, h[0]);
	Split(a1, _ch, _cl);
	Two_Product_2Presplit(a1, _ch, _cl, b0, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _k, _1);
	Quick_Two_Sum(_j, _k, _l, _2);
	Split(b1, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b1, _bh, _bl, _i, _0);
	Two_Sum(_1, _0, _k, h[1]);
	Two_Sum(_2, _k, _j, _1);
	Two_Sum(_l, _j, _m, _2);
	Two_Product_2Presplit(a1, _ch, _cl, b1, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _n, _0);
	Two_Sum(_1, _0, _i, h[2]);
	Two_Sum(_2, _i, _k, _1);
	Two_Sum(_m, _k, _l, _2);
	Two_Sum(_j, _n, _k, _0);
	Two_Sum(_1, _0, _j, h[3]);
	Two_Sum(_2, _j, _i, _1);
	Two_Sum(_l, _i, _m, _2);
	Two_Sum(_1, _k, _i, h[4]);
	Two_Sum(_2, _i, _k, h[5]);
	Two_Sum(_m, _k, h[7], h[6]);
}


int expansionObject::Sub_product(const int alen, const double *a, const int blen, const double *b, double *h)
{
	if (alen == 1) return Gen_Scale(blen, b, a[0], h);
	else
	{
		const double* a1 = a;
		int a1len = alen / 2;
		const double* a2 = a1 + a1len;
		int a2len = alen - a1len;

		int a1blen, a2blen;
		double *a1b = (double *)malloc(2 * a1len * blen * sizeof(double));
		double *a2b = (double *)malloc(2 * a2len * blen * sizeof(double));
		a1blen = Sub_product(a1len, a1, blen, b, a1b);
		a2blen = Sub_product(a2len, a2, blen, b, a2b);
		int hlen = Gen_Sum(a1blen, a1b, a2blen, a2b, h);
		free(a1b);
		free(a2b);
		return hlen;
	}
}


int expansionObject::Gen_Product(const int alen, const double *a, const int blen, const double *b, double *h)
{
	if (alen == 1)
	{
		if (blen == 1) { Two_Prod(a[0], b[0], h); return 2; } else if (blen == 2) { Two_One_Prod(b, a[0], h); return 4; } else return Gen_Scale(blen, b, a[0], h);
	} else if (alen == 2)
	{
		if (blen == 1) { Two_One_Prod(a, b[0], h); return 4; } else if (blen == 2) { Two_Two_Prod(a[1], a[0], b[1], b[0], h); return 8; } else return Sub_product(alen, a, blen, b, h);
	} else
	{
		if (blen == 1) return Gen_Scale(alen, a, b[0], h);
		else if (alen < blen) return Sub_product(alen, a, blen, b, h);
		else return Sub_product(blen, b, alen, a, h);
	}
}


double expansionObject::To_Double(const int elen, const double *e)
{
  double Q = e[0];
  for (int e_i = 1; e_i < elen; e_i++) Q += e[e_i];
  return Q;
}

