/****************************************************************************
* orientation.cpp                                                           *
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

#include "basics.h"
#include "expansion.h"

//#pragma optimize("", off)

namespace T_MESH
{

inline void supo3d1(expansionObject& o, 
	double *c1, double *c2, double *c3, double *c4, double *c5, double *c6,
	double *a1, double *a2, double& i, 	double *k1, double *k2, double *k3, 
	double *k4, int& l1, int& l2)
{
	if (c1[0] == 0.0) {
		if (c2[0] == 0.0) {
			a1[0] = a2[0] = 0.0;
			l1 = l2 = 1;
		}
		else {
			i = -c2[0];
			o.Two_Prod(i, c3[1], a1);
			o.Two_Prod(c2[0], c4[1], a2);
			l1 = l2 = 2;
		}
	}
	else {
		if (c2[0] == 0.0) {
			i = -c1[0];
			o.Two_Prod(c1[0], c5[1], a1);
			o.Two_Prod(i, c6[1], a2);
			l1 = l2 = 2;
		}
		else {
			o.Two_Prod(c1[0], c5[1], k1);
			o.Two_Prod(c2[0], c3[1], k2);
			o.Two_Two_Diff(k1, k2, a1);
			o.Two_Prod(c2[0], c4[1], k3);
			o.Two_Prod(c1[0], c6[1], k4);
			o.Two_Two_Diff(k3, k4, a2);
			l1 = l2 = 4;
		}
	}
}

inline void supo3d2(expansionObject& o,
	double *c1, double *c2, double *c3, double *c4, double *u, 
	int& fl, double fin[2][192], int& wh,
	double *c5, double& i, double *c6, double *c7)
	
{
	if (c1[0] != 0.0) {
		if (c2[0] != 0.0) {
			o.Two_Prod(c1[0], c2[0], c3);
			o.Two_One_Prod(c3, c4[1], u);
			fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
			wh = !wh;
			if (c4[0] != 0.0) {
				o.Two_One_Prod(c3, c4[0], u);
				fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
				wh = !wh;
			}
		}
		if (c5[0] != 0.0) {
			i = -c1[0];
			o.Two_Prod(i, c5[0], c6);
			o.Two_One_Prod(c6, c7[1], u);
			fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
			wh = !wh;
			if (c7[0] != 0.0) {
				o.Two_One_Prod(c6, c7[0], u);
				fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
				wh = !wh;
			}
		}
	}
}

double orient2d_expansions(double* pa, double* pb, double* pc, double dsm)
{
	double acx[2], acy[2], bcx[2], bcy[2], dtl[2], dtr[2], B[4];
	double s[2], t[2], u[4], C1[8], C2[12], D[16];
	int C1l, C2l, Dl;
	expansionObject o;

	acx[1] = (pa[0] - pc[0]);
	bcx[1] = (pb[0] - pc[0]);
	acy[1] = (pa[1] - pc[1]);
	bcy[1] = (pb[1] - pc[1]);

	o.Two_Prod(acx[1], bcy[1], dtl);
	o.Two_Prod(acy[1], bcx[1], dtr);
	o.Two_Two_Diff(dtl, dtr, B);

	double det = expansionObject::To_Double(4, B);
	double eb = 7.7715611723761027e-016 * dsm;
	if ((det >= eb) || (-det >= eb)) return det;

	o.Two_Diff_Back(pa[0], pc[0], acx);
	o.Two_Diff_Back(pb[0], pc[0], bcx);
	o.Two_Diff_Back(pa[1], pc[1], acy);
	o.Two_Diff_Back(pb[1], pc[1], bcy);

	if ((acx[0] == 0.0) && (acy[0] == 0.0) && (bcx[0] == 0.0) && (bcy[0] == 0.0)) return det;

	eb = 1.1093356479670487e-031 * dsm + 1.1102230246251565e-016 * FABS(det);
	det += (acx[1] * bcy[0] + bcy[1] * acx[0]) - (acy[1] * bcx[0] + bcx[1] * acy[0]);
	if ((det >= eb) || (-det >= eb)) return det;

	o.Two_Prod(acx[0], bcy[1], s);
	o.Two_Prod(acy[0], bcx[1], t);
	o.Two_Two_Diff(s, t, u);
	C1l = o.Gen_Sum(4, B, 4, u, C1);

	o.Two_Prod(acx[1], bcy[0], s);
	o.Two_Prod(acy[1], bcx[0], t);
	o.Two_Two_Diff(s, t, u);
	C2l = o.Gen_Sum(C1l, C1, 4, u, C2);

	o.Two_Prod(acx[0], bcy[0], s);
	o.Two_Prod(acy[0], bcx[0], t);
	o.Two_Two_Diff(s, t, u);
	Dl = o.Gen_Sum(C2l, C2, 4, u, D);

	return(D[Dl - 1]);
}

double orient3d_expansions(double* pa, double* pb, double* pc, double* pd)
{
	double fadx, fbdx, fcdx, fady, fbdy, fcdy, fadz, fbdz, fcdz, pm, eb;
	double fbdxcdy, fcdxbdy, fcdxady, fadxcdy, fadxbdy, fbdxady, det;

	fadx = pa[0] - pd[0]; fbdx = pb[0] - pd[0]; fcdx = pc[0] - pd[0];
	fady = pa[1] - pd[1]; fbdy = pb[1] - pd[1]; fcdy = pc[1] - pd[1];
	fadz = pa[2] - pd[2]; fbdz = pb[2] - pd[2]; fcdz = pc[2] - pd[2];

	fbdxcdy = fbdx * fcdy; fcdxbdy = fcdx * fbdy;
	fcdxady = fcdx * fady; fadxcdy = fadx * fcdy;
	fadxbdy = fadx * fbdy; fbdxady = fbdx * fady;

	det = fadz * (fbdxcdy - fcdxbdy) + fbdz * (fcdxady - fadxcdy) + fcdz * (fadxbdy - fbdxady);
	pm = (FABS(fbdxcdy) + FABS(fcdxbdy)) * FABS(fadz) + (FABS(fcdxady) + FABS(fadxcdy)) * FABS(fbdz) + (FABS(fadxbdy) + FABS(fbdxady)) * FABS(fcdz);
	eb = 3.3306690738754716e-016 * pm;
	if ((det > eb) || (-det > eb)) return det;

	double adx[2], bdx[2], cdx[2], ady[2], bdy[2], cdy[2], adz[2], bdz[2], cdz[2];
	double bdxcdy[2], cdxbdy[2], cdxady[2], adxcdy[2], adxbdy[2], bdxady[2];
	double bc[4], ca[4], ab[4];
	double bdxt_cdy[2], cdxt_bdy[2], cdxt_ady[2];
	double adxt_cdy[2], adxt_bdy[2], bdxt_ady[2];
	double bdyt_cdx[2], cdyt_bdx[2], cdyt_adx[2];
	double adyt_cdx[2], adyt_bdx[2], bdyt_adx[2];
	double bdxt_cdyt[2], cdxt_bdyt[2], cdxt_adyt[2];
	double adxt_cdyt[2], adxt_bdyt[2], bdxt_adyt[2];
	double u[4], v[12], w[16];
	double adet[8], bdet[8], cdet[8], abdet[16];
	double fin[2][192];
	int wh = 0;
	double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
	double bct[8], cat[8], abt[8];
	int alen, blen, clen, finlen, vlen, wlen;
	int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
	int bctlen, catlen, abtlen;
	int ablen;
	double inv;

	expansionObject o;

	adx[1] = fadx;
	bdx[1] = fbdx;
	cdx[1] = fcdx;
	ady[1] = fady;
	bdy[1] = fbdy;
	cdy[1] = fcdy;
	adz[1] = fadz;
	bdz[1] = fbdz;
	cdz[1] = fcdz;

	o.Two_Prod(bdx[1], cdy[1], bdxcdy);
	o.Two_Prod(cdx[1], bdy[1], cdxbdy);
	o.Two_Two_Diff(bdxcdy, cdxbdy, bc);
	alen = o.Gen_Scale(4, bc, adz[1], adet);

	o.Two_Prod(cdx[1], ady[1], cdxady);
	o.Two_Prod(adx[1], cdy[1], adxcdy);
	o.Two_Two_Diff(cdxady, adxcdy, ca);
	blen = o.Gen_Scale(4, ca, bdz[1], bdet);

	o.Two_Prod(adx[1], bdy[1], adxbdy);
	o.Two_Prod(bdx[1], ady[1], bdxady);
	o.Two_Two_Diff(adxbdy, bdxady, ab);
	clen = o.Gen_Scale(4, ab, cdz[1], cdet);

	ablen = o.Gen_Sum(alen, adet, blen, bdet, abdet);
	finlen = o.Gen_Sum(ablen, abdet, clen, cdet, fin[wh]);

	det = expansionObject::To_Double(finlen, fin[wh]);
	eb = 3.3306690738754731e-016 * pm;
	if ((det >= eb) || (-det >= eb)) return det;

	o.Two_Diff_Back(pa[0], pd[0], adx);
	o.Two_Diff_Back(pb[0], pd[0], bdx);
	o.Two_Diff_Back(pc[0], pd[0], cdx);
	o.Two_Diff_Back(pa[1], pd[1], ady);
	o.Two_Diff_Back(pb[1], pd[1], bdy);
	o.Two_Diff_Back(pc[1], pd[1], cdy);
	o.Two_Diff_Back(pa[2], pd[2], adz);
	o.Two_Diff_Back(pb[2], pd[2], bdz);
	o.Two_Diff_Back(pc[2], pd[2], cdz);

	if ((adx[0] == 0.0) && (bdx[0] == 0.0) && (cdx[0] == 0.0) &&
		(ady[0] == 0.0) && (bdy[0] == 0.0) && (cdy[0] == 0.0) &&
		(adz[0] == 0.0) && (bdz[0] == 0.0) && (cdz[0] == 0.0)) return det;

	eb = 3.2047474274603644e-031 * pm + 1.1102230246251565e-016 * FABS(det);
	det += (adz[1] * ((bdx[1] * cdy[0] + cdy[1] * bdx[0])
		- (bdy[1] * cdx[0] + cdx[1] * bdy[0]))
		+ adz[0] * (bdx[1] * cdy[1] - bdy[1] * cdx[1]))
		+ (bdz[1] * ((cdx[1] * ady[0] + ady[1] * cdx[0])
			- (cdy[1] * adx[0] + adx[1] * cdy[0]))
			+ bdz[0] * (cdx[1] * ady[1] - cdy[1] * adx[1]))
		+ (cdz[1] * ((adx[1] * bdy[0] + bdy[1] * adx[0])
			- (ady[1] * bdx[0] + bdx[1] * ady[0]))
			+ cdz[0] * (adx[1] * bdy[1] - ady[1] * bdx[1]));
	if ((det >= eb) || (-det >= eb)) return det;

	// Filters did not work. Compute exactly...
	supo3d1(o, adx, ady, bdx, cdx, bdy, cdy, at_b, at_c, inv,
		adxt_bdy, adyt_bdx, adyt_cdx, adxt_cdy, at_blen, at_clen);

	supo3d1(o, bdx, bdy, cdx, adx, cdy, ady, bt_c, bt_a, inv,
		bdxt_cdy, bdyt_cdx, bdyt_adx, bdxt_ady, bt_alen, bt_clen);

	supo3d1(o, cdx, cdy, adx, bdx, ady, bdy, ct_a, ct_b, inv,
		cdxt_ady, cdyt_adx, cdyt_bdx, cdxt_bdy, ct_alen, ct_blen);

	bctlen = o.Gen_Sum(bt_clen, bt_c, ct_blen, ct_b, bct);
	wlen = o.Gen_Scale(bctlen, bct, adz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]); wh = !wh;

	catlen = o.Gen_Sum(ct_alen, ct_a, at_clen, at_c, cat);
	wlen = o.Gen_Scale(catlen, cat, bdz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]); wh = !wh;

	abtlen = o.Gen_Sum(at_blen, at_b, bt_alen, bt_a, abt);
	wlen = o.Gen_Scale(abtlen, abt, cdz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]); wh = !wh;

	if (adz[0] != 0.0) {
		vlen = o.Gen_Scale(4, bc, adz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]); wh = !wh;
	}
	if (bdz[0] != 0.0) {
		vlen = o.Gen_Scale(4, ca, bdz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]); wh = !wh;
	}
	if (cdz[0] != 0.0) {
		vlen = o.Gen_Scale(4, ab, cdz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]); wh = !wh;
	}

	supo3d2(o, adx, bdy, adxt_bdyt, cdz, u, finlen, fin, wh, cdy, inv, adxt_cdyt, bdz);
	supo3d2(o, bdx, cdy, bdxt_cdyt, adz, u, finlen, fin, wh, ady, inv, bdxt_adyt, cdz);
	supo3d2(o, cdx, ady, cdxt_adyt, bdz, u, finlen, fin, wh, bdy, inv, cdxt_bdyt, adz);

	if (adz[0] != 0.0) {
		wlen = o.Gen_Scale(bctlen, bct, adz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]); wh = !wh;
	}
	if (bdz[0] != 0.0) {
		wlen = o.Gen_Scale(catlen, cat, bdz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]); wh = !wh;
	}
	if (cdz[0] != 0.0) {
		wlen = o.Gen_Scale(abtlen, abt, cdz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]);	wh = !wh;
	}

	return fin[wh][finlen - 1];
}


/*****************************************************************************/
/*                                                                           */
/*  PUBLIC FUNCTIONS                                                         */
/*  tri_orientation()    Computes the orientation of three 2D points.        */
/*  tet_orientation()    Computes the orientation of four 3D points.         */
/*  insphere()    Computes whether a point in the sphere by other four.      */
/*  incircle3D()    Computes whether a point in the circle by other three.   */
/*                                                                           */
/*****************************************************************************/

double TMesh::tri_orientation(double *pa, double *pb, double *pc)
{
 double dte_left, dte_right, det, dsm, eb;

 dte_left = (pa[0]-pc[0])*(pb[1]-pc[1]);
 dte_right = (pa[1]-pc[1])*(pb[0]-pc[0]);
 det = dte_left - dte_right;

 if (dte_left > 0.0) {if (dte_right <= 0.0) return det; else dsm = dte_left + dte_right;}
 else if (dte_left < 0.0) {if (dte_right >= 0.0) return det; else dsm = -dte_left - dte_right;}
 else return det;

 eb = 3.3306690738754706e-016*dsm;
 if ((det>=eb) || (-det>=eb)) return det;

 return orient2d_expansions(pa, pb, pc, dsm);
}

double TMesh::tet_orientation(double *pa, double *pb, double *pc, double *pd)
{
 double fadx, fbdx, fcdx, fady, fbdy, fcdy, fadz, fbdz, fcdz, pm, eb;
 double fbdxcdy, fcdxbdy, fcdxady, fadxcdy, fadxbdy, fbdxady, det;

 fadx = pa[0]-pd[0]; fbdx = pb[0]-pd[0]; fcdx = pc[0]-pd[0];
 fady = pa[1]-pd[1]; fbdy = pb[1]-pd[1]; fcdy = pc[1]-pd[1];
 fadz = pa[2]-pd[2]; fbdz = pb[2]-pd[2]; fcdz = pc[2]-pd[2];

 fbdxcdy = fbdx*fcdy; fcdxbdy = fcdx*fbdy;
 fcdxady = fcdx*fady; fadxcdy = fadx*fcdy;
 fadxbdy = fadx*fbdy; fbdxady = fbdx*fady;

 det = fadz*(fbdxcdy-fcdxbdy)+fbdz*(fcdxady-fadxcdy)+fcdz*(fadxbdy-fbdxady);
 pm=(FABS(fbdxcdy)+FABS(fcdxbdy))*FABS(fadz)+(FABS(fcdxady)+FABS(fadxcdy))*FABS(fbdz)+(FABS(fadxbdy)+FABS(fbdxady))*FABS(fcdz);
 eb = 3.3306690738754716e-016*pm;
 if ((det>eb) || (-det>eb)) return det;

 return orient3d_expansions(pa, pb, pc, pd);
}


double insphere_filtered(double *pa, double *pb, double *pc, double *pd, double *pe)
{
	double aex = pa[0] - pe[0], bex = pb[0] - pe[0], cex = pc[0] - pe[0], dex = pd[0] - pe[0];
	double aey = pa[1] - pe[1], bey = pb[1] - pe[1], cey = pc[1] - pe[1], dey = pd[1] - pe[1];
	double aez = pa[2] - pe[2], bez = pb[2] - pe[2], cez = pc[2] - pe[2], dez = pd[2] - pe[2];

	double aexbey = aex * bey, bexaey = bex * aey, ab = aexbey - bexaey;
	double bexcey = bex * cey, cexbey = cex * bey, bc = bexcey - cexbey;
	double cexdey = cex * dey, dexcey = dex * cey, cd = cexdey - dexcey;
	double dexaey = dex * aey, aexdey = aex * dey, da = dexaey - aexdey;
	double aexcey = aex * cey, cexaey = cex * aey, ac = aexcey - cexaey;
	double bexdey = bex * dey, dexbey = dex * bey, bd = bexdey - dexbey;

	double abc = aez * bc - bez * ac + cez * ab;
	double bcd = bez * cd - cez * bd + dez * bc;
	double cda = cez * da + dez * ac + aez * cd;
	double dab = dez * ab + aez * bd + bez * da;

	double alift = aex * aex + aey * aey + aez * aez;
	double blift = bex * bex + bey * bey + bez * bez;
	double clift = cex * cex + cey * cey + cez * cez;
	double dlift = dex * dex + dey * dey + dez * dez;

	double det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

	TO_NONNEGATIVE(aez); TO_NONNEGATIVE(bez); TO_NONNEGATIVE(cez); TO_NONNEGATIVE(dez);
	TO_NONNEGATIVE(aexbey); TO_NONNEGATIVE(bexaey);	TO_NONNEGATIVE(bexcey);	TO_NONNEGATIVE(cexbey);
	TO_NONNEGATIVE(cexdey);	TO_NONNEGATIVE(dexcey);	TO_NONNEGATIVE(dexaey);	TO_NONNEGATIVE(aexdey);
	TO_NONNEGATIVE(aexcey);	TO_NONNEGATIVE(cexaey);	TO_NONNEGATIVE(bexdey);	TO_NONNEGATIVE(dexbey);
	double sz =   ((cexdey + dexcey) * bez + (dexbey + bexdey) * cez + (bexcey + cexbey) * dez) * alift +
					((dexaey + aexdey) * cez + (aexcey + cexaey) * dez + (cexdey + dexcey) * aez) * blift +
					((aexbey + bexaey) * dez + (bexdey + dexbey) * aez + (dexaey + aexdey) * bez) * clift +
					((bexcey + cexbey) * aez + (cexaey + aexcey) * bez + (aexbey + bexaey) * cez) * dlift;

	double errbound = 1.7763568394002532e-015 * sz;
	if ((det > errbound) || (-det > errbound)) return det;
	else return 0;
}

double insphere_exact(double *pa, double *pb, double *pc, double *pd, double *pe)
{
	double axby[2], bxcy[2], cxdy[2], dxey[2], exay[2];
	double bxay[2], cxby[2], dxcy[2], exdy[2], axey[2];
	double axcy[2], bxdy[2], cxey[2], dxay[2], exby[2];
	double cxay[2], dxby[2], excy[2], axdy[2], bxey[2];
	double xab[4], xbc[4], xcd[4], xde[4], xea[4];
	double xac[4], xbd[4], ce[4], xda[4], eb[4];
	double gvec8a[8], gvec8b[8], gvec16[16];
	int gvec8alen, gvec8blen, gvec16len;
	double xabc[24], xbcd[24], xcde[24], dea[24], eab[24];
	double abd[24], bce[24], xcda[24], deb[24], eac[24];
	int abclen, bcdlen, cdelen, dealen, eablen;
	int abdlen, bcelen, cdalen, deblen, eaclen;
	double gvec48a[48], gvec48b[48];
	int gvec48alen, gvec48blen;
	double abcd[96], bcde[96], cdea[96], deab[96], eabc[96];
	int abcdlen, bcdelen, cdealen, deablen, eabclen;
	double gvec192[192];
	double det384x[384], det384y[384], det384z[384];
	int xlen, ylen, zlen;
	double detxy[768];
	int xylen;
	double adet[1152], bdet[1152], cdet[1152], ddet[1152], edet[1152];
	int alen, blen, clen, dlen, elen;
	double abdet[2304], cddet[2304], cdedet[3456];
	int ablen, cdlen;
	double deter[5760];
	int deterlen;

	expansionObject o;

	// Exact
	o.Two_Prod(pa[0], pb[1], axby);
	o.Two_Prod(pb[0], pa[1], bxay);
	o.Two_Two_Diff(axby, bxay, xab);

	o.Two_Prod(pb[0], pc[1], bxcy);
	o.Two_Prod(pc[0], pb[1], cxby);
	o.Two_Two_Diff(bxcy, cxby, xbc);

	o.Two_Prod(pc[0], pd[1], cxdy);
	o.Two_Prod(pd[0], pc[1], dxcy);
	o.Two_Two_Diff(cxdy, dxcy, xcd);

	o.Two_Prod(pd[0], pe[1], dxey);
	o.Two_Prod(pe[0], pd[1], exdy);
	o.Two_Two_Diff(dxey, exdy, xde);

	o.Two_Prod(pe[0], pa[1], exay);
	o.Two_Prod(pa[0], pe[1], axey);
	o.Two_Two_Diff(exay, axey, xea);

	o.Two_Prod(pa[0], pc[1], axcy);
	o.Two_Prod(pc[0], pa[1], cxay);
	o.Two_Two_Diff(axcy, cxay, xac);

	o.Two_Prod(pb[0], pd[1], bxdy);
	o.Two_Prod(pd[0], pb[1], dxby);
	o.Two_Two_Diff(bxdy, dxby, xbd);

	o.Two_Prod(pc[0], pe[1], cxey);
	o.Two_Prod(pe[0], pc[1], excy);
	o.Two_Two_Diff(cxey, excy, ce);

	o.Two_Prod(pd[0], pa[1], dxay);
	o.Two_Prod(pa[0], pd[1], axdy);
	o.Two_Two_Diff(dxay, axdy, xda);

	o.Two_Prod(pe[0], pb[1], exby);
	o.Two_Prod(pb[0], pe[1], bxey);
	o.Two_Two_Diff(exby, bxey, eb);

	gvec8alen = o.Gen_Scale(4, xbc, pa[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xac, -pb[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xab, pc[2], gvec8a);
	abclen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, xabc);

	gvec8alen = o.Gen_Scale(4, xcd, pb[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xbd, -pc[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xbc, pd[2], gvec8a);
	bcdlen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, xbcd);

	gvec8alen = o.Gen_Scale(4, xde, pc[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, ce, -pd[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xcd, pe[2], gvec8a);
	cdelen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, xcde);

	gvec8alen = o.Gen_Scale(4, xea, pd[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xda, -pe[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xde, pa[2], gvec8a);
	dealen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, dea);

	gvec8alen = o.Gen_Scale(4, xab, pe[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, eb, -pa[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xea, pb[2], gvec8a);
	eablen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, eab);

	gvec8alen = o.Gen_Scale(4, xbd, pa[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xda, pb[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xab, pd[2], gvec8a);
	abdlen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, abd);

	gvec8alen = o.Gen_Scale(4, ce, pb[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, eb, pc[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xbc, pe[2], gvec8a);
	bcelen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, bce);

	gvec8alen = o.Gen_Scale(4, xda, pc[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xac, pd[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xcd, pa[2], gvec8a);
	cdalen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, xcda);

	gvec8alen = o.Gen_Scale(4, eb, pd[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, xbd, pe[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xde, pb[2], gvec8a);
	deblen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, deb);

	gvec8alen = o.Gen_Scale(4, xac, pe[2], gvec8a);
	gvec8blen = o.Gen_Scale(4, ce, pa[2], gvec8b);
	gvec16len = o.Gen_Sum(gvec8alen, gvec8a, gvec8blen, gvec8b, gvec16);
	gvec8alen = o.Gen_Scale(4, xea, pc[2], gvec8a);
	eaclen = o.Gen_Sum(gvec8alen, gvec8a, gvec16len, gvec16, eac);

	gvec48alen = o.Gen_Sum(cdelen, xcde, bcelen, bce, gvec48a);
	gvec48blen = o.Gen_Sum(deblen, deb, bcdlen, xbcd, gvec48b);
	o.Gen_Invert(gvec48blen, gvec48b);

	bcdelen = o.Gen_Sum(gvec48alen, gvec48a, gvec48blen, gvec48b, bcde);
	xlen = o.Gen_Scale(bcdelen, bcde, pa[0], gvec192);
	xlen = o.Gen_Scale(xlen, gvec192, pa[0], det384x);
	ylen = o.Gen_Scale(bcdelen, bcde, pa[1], gvec192);
	ylen = o.Gen_Scale(ylen, gvec192, pa[1], det384y);
	zlen = o.Gen_Scale(bcdelen, bcde, pa[2], gvec192);
	zlen = o.Gen_Scale(zlen, gvec192, pa[2], det384z);
	xylen = o.Gen_Sum(xlen, det384x, ylen, det384y, detxy);
	alen = o.Gen_Sum(xylen, detxy, zlen, det384z, adet);

	gvec48alen = o.Gen_Sum(dealen, dea, cdalen, xcda, gvec48a);
	gvec48blen = o.Gen_Sum(eaclen, eac, cdelen, xcde, gvec48b);
	o.Gen_Invert(gvec48blen, gvec48b);

	cdealen = o.Gen_Sum(gvec48alen, gvec48a, gvec48blen, gvec48b, cdea);
	xlen = o.Gen_Scale(cdealen, cdea, pb[0], gvec192);
	xlen = o.Gen_Scale(xlen, gvec192, pb[0], det384x);
	ylen = o.Gen_Scale(cdealen, cdea, pb[1], gvec192);
	ylen = o.Gen_Scale(ylen, gvec192, pb[1], det384y);
	zlen = o.Gen_Scale(cdealen, cdea, pb[2], gvec192);
	zlen = o.Gen_Scale(zlen, gvec192, pb[2], det384z);
	xylen = o.Gen_Sum(xlen, det384x, ylen, det384y, detxy);
	blen = o.Gen_Sum(xylen, detxy, zlen, det384z, bdet);

	gvec48alen = o.Gen_Sum(eablen, eab, deblen, deb, gvec48a);
	gvec48blen = o.Gen_Sum(abdlen, abd, dealen, dea, gvec48b);
	o.Gen_Invert(gvec48blen, gvec48b);

	deablen = o.Gen_Sum(gvec48alen, gvec48a, gvec48blen, gvec48b, deab);
	xlen = o.Gen_Scale(deablen, deab, pc[0], gvec192);
	xlen = o.Gen_Scale(xlen, gvec192, pc[0], det384x);
	ylen = o.Gen_Scale(deablen, deab, pc[1], gvec192);
	ylen = o.Gen_Scale(ylen, gvec192, pc[1], det384y);
	zlen = o.Gen_Scale(deablen, deab, pc[2], gvec192);
	zlen = o.Gen_Scale(zlen, gvec192, pc[2], det384z);
	xylen = o.Gen_Sum(xlen, det384x, ylen, det384y, detxy);
	clen = o.Gen_Sum(xylen, detxy, zlen, det384z, cdet);

	gvec48alen = o.Gen_Sum(abclen, xabc, eaclen, eac, gvec48a);
	gvec48blen = o.Gen_Sum(bcelen, bce, eablen, eab, gvec48b);
	o.Gen_Invert(gvec48blen, gvec48b);

	eabclen = o.Gen_Sum(gvec48alen, gvec48a, gvec48blen, gvec48b, eabc);
	xlen = o.Gen_Scale(eabclen, eabc, pd[0], gvec192);
	xlen = o.Gen_Scale(xlen, gvec192, pd[0], det384x);
	ylen = o.Gen_Scale(eabclen, eabc, pd[1], gvec192);
	ylen = o.Gen_Scale(ylen, gvec192, pd[1], det384y);
	zlen = o.Gen_Scale(eabclen, eabc, pd[2], gvec192);
	zlen = o.Gen_Scale(zlen, gvec192, pd[2], det384z);
	xylen = o.Gen_Sum(xlen, det384x, ylen, det384y, detxy);
	dlen = o.Gen_Sum(xylen, detxy, zlen, det384z, ddet);

	gvec48alen = o.Gen_Sum(bcdlen, xbcd, abdlen, abd, gvec48a);
	gvec48blen = o.Gen_Sum(cdalen, xcda, abclen, xabc, gvec48b);
	o.Gen_Invert(gvec48blen, gvec48b);

	abcdlen = o.Gen_Sum(gvec48alen, gvec48a, gvec48blen, gvec48b, abcd);
	xlen = o.Gen_Scale(abcdlen, abcd, pe[0], gvec192);
	xlen = o.Gen_Scale(xlen, gvec192, pe[0], det384x);
	ylen = o.Gen_Scale(abcdlen, abcd, pe[1], gvec192);
	ylen = o.Gen_Scale(ylen, gvec192, pe[1], det384y);
	zlen = o.Gen_Scale(abcdlen, abcd, pe[2], gvec192);
	zlen = o.Gen_Scale(zlen, gvec192, pe[2], det384z);
	xylen = o.Gen_Sum(xlen, det384x, ylen, det384y, detxy);
	elen = o.Gen_Sum(xylen, detxy, zlen, det384z, edet);

	ablen = o.Gen_Sum(alen, adet, blen, bdet, abdet);
	cdlen = o.Gen_Sum(clen, cdet, dlen, ddet, cddet);
	cdelen = o.Gen_Sum(cdlen, cddet, elen, edet, cdedet);
	deterlen = o.Gen_Sum(ablen, abdet, cdelen, cdedet, deter);

	return deter[deterlen - 1];
}

double TMesh::insphere(double *pa, double *pb, double *pc, double *pd, double *pe)
{
	double result = insphere_filtered(pa, pb, pc, pd, pe);
	if (result) return result;
	return insphere_exact(pa, pb, pc, pd, pe);
}


inline int incircle3D_filtered(const double* p0, const double* p1, const double* p2, const double* p3)
{
	double p10x = p1[0] - p0[0], p10y = p1[1] - p0[1], p10z = p1[2] - p0[2];
	double l1 = p10x * p10x + p10y * p10y + p10z * p10z;
	double p20x = p2[0] - p0[0], p20y = p2[1] - p0[1], p20z = p2[2] - p0[2];
	double l2 = p20x * p20x + p20y * p20y + p20z * p20z;

	double ap10x = FABS(p10x), ap10y = FABS(p10y), ap10z = FABS(p10z);
	double ap20x = FABS(p20x), ap20y = FABS(p20y), ap20z = FABS(p20z);
	double max1 = MAX(ap10x, MAX(ap10y, ap10z));
	double max2 = MAX(max1, MAX(ap20x, MAX(ap20y, ap20z)));
	if ((max1 < 2.22985945097100191780e-74) || (max2 > 2.59614842926741294957e+33)) return 0;

	double p30x = p3[0] - p0[0], p30y = p3[1] - p0[1], p30z = p3[2] - p0[2];
	double l3 = p30x * p30x + p30y * p30y + p30z * p30z;

	double a11 = l1 * 2;
	double a12 = (p10x * p20x + p10y * p20y + p10z * p20z)*2;
	double a22 = l2 * 2;
	double a31 = (p30x * p10x + p30y * p10y + p30z * p10z)*2;
	double a32 = (p30x * p20x + p30y * p20y + p30z * p20z)*2;

	double eps = (8.99983341597279045654e-14 * (((max1 * max2) * max2) * max2));
	double det = (a11 * a22) - (a12 * a12);
	int det_sign;
	if ((det > eps)) det_sign = 1;
	else if ((det < -eps)) det_sign = -1;
	else return 0;
	
	TO_NONNEGATIVE(p30x); TO_NONNEGATIVE(p30y); TO_NONNEGATIVE(p30z);
	double max4 = MAX(max1, MAX(p30x, MAX(p30y, p30z)));
	double max3 = MAX(max2, max4);
	if ((MIN(max2, max4) < 4.84416636653081796592e-50) || (max3 > 2.59614842926741294957e+33)) return 0;

	eps = (1.72198804259438718181e-12 * (((((max4 * max2) * max2) * max2) * max3) * max3));
	double r = ((det * l3) - (((a31 * ((a22 * l1) - (a12 * l2)))) + (a32 * ((a11 * l2) - (a12 * l1)))));
	if ((r > eps)) return -det_sign;
	else if ((r < -eps)) return det_sign;
	else return 0;
}

int incircle3D_exact(const double* p0, const double* p1, const double* p2, const double* p3)
{
	double p1x_p0x[2], p1y_p0y[2], p1z_p0z[2], p2x_p0x[2], p2y_p0y[2], p2z_p0z[2], p3x_p0x[2], p3y_p0y[2], p3z_p0z[2];
	double a_1[8], a_2[8], a_3[8];
	double l_part[16];
	double l1[18], l2[18], l3[18];
	double a_part2[24];
	double a11[24], a22[24], a12[24], a31[24], a32[24];
	double r_1[864], r_3[864];
	double r_2[648], r_4[648], det_1[648];
	double det_2[1152], rr_1[1152], rr_2[1152];
	double det[1800];

	int l_part_len, a_part2_len;

	expansionObject o;

	o.two_Diff(p1[0], p0[0], p1x_p0x);
	o.two_Diff(p1[1], p0[1], p1y_p0y);
	o.two_Diff(p1[2], p0[2], p1z_p0z);
	o.Two_Square(p1x_p0x[1], p1x_p0x[0], a_1);
	o.Two_Square(p1y_p0y[1], p1y_p0y[0], a_2);
	o.Two_Square(p1z_p0z[1], p1z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(6, a_1, 6, a_2, l_part);
	int l1_len = o.Gen_Sum(l_part_len, l_part, 6, a_3, l1);

	o.two_Diff(p2[0], p0[0], p2x_p0x);
	o.two_Diff(p2[1], p0[1], p2y_p0y);
	o.two_Diff(p2[2], p0[2], p2z_p0z);
	o.Two_Square(p2x_p0x[1], p2x_p0x[0], a_1);
	o.Two_Square(p2y_p0y[1], p2y_p0y[0], a_2);
	o.Two_Square(p2z_p0z[1], p2z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(6, a_1, 6, a_2, l_part);
	int l2_len = o.Gen_Sum(l_part_len, l_part, 6, a_3, l2);

	o.two_Diff(p3[0], p0[0], p3x_p0x);
	o.two_Diff(p3[1], p0[1], p3y_p0y);
	o.two_Diff(p3[2], p0[2], p3z_p0z);
	o.Two_Square(p3x_p0x[1], p3x_p0x[0], a_1);
	o.Two_Square(p3y_p0y[1], p3y_p0y[0], a_2);
	o.Two_Square(p3z_p0z[1], p3z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(6, a_1, 6, a_2, l_part);
	int l3_len = o.Gen_Sum(l_part_len, l_part, 6, a_3, l3);

	int a11_len = l1_len; o.Double(l1_len, l1, a11);
	int a22_len = l2_len; o.Double(l2_len, l2, a22);

	o.Two_Two_Prod(p1x_p0x[1], p1x_p0x[0], p2x_p0x[1], p2x_p0x[0], a_1);
	o.Two_Two_Prod(p1y_p0y[1], p1y_p0y[0], p2y_p0y[1], p2y_p0y[0], a_2);
	o.Two_Two_Prod(p1z_p0z[1], p1z_p0z[0], p2z_p0z[1], p2z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(8, a_1, 8, a_2, l_part);
	a_part2_len = o.Gen_Sum(l_part_len, l_part, 8, a_3, a_part2);
	int a12_len = a_part2_len; o.Double(a_part2_len, a_part2, a12);

	o.Two_Two_Prod(p3x_p0x[1], p3x_p0x[0], p1x_p0x[1], p1x_p0x[0], a_1);
	o.Two_Two_Prod(p3y_p0y[1], p3y_p0y[0], p1y_p0y[1], p1y_p0y[0], a_2);
	o.Two_Two_Prod(p3z_p0z[1], p3z_p0z[0], p1z_p0z[1], p1z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(8, a_1, 8, a_2, l_part);
	a_part2_len = o.Gen_Sum(l_part_len, l_part, 8, a_3, a_part2);
	int a31_len = a_part2_len; o.Double(a_part2_len, a_part2, a31);

	o.Two_Two_Prod(p3x_p0x[1], p3x_p0x[0], p2x_p0x[1], p2x_p0x[0], a_1);
	o.Two_Two_Prod(p3y_p0y[1], p3y_p0y[0], p2y_p0y[1], p2y_p0y[0], a_2);
	o.Two_Two_Prod(p3z_p0z[1], p3z_p0z[0], p2z_p0z[1], p2z_p0z[0], a_3);
	l_part_len = o.Gen_Sum(8, a_1, 8, a_2, l_part);
	a_part2_len = o.Gen_Sum(l_part_len, l_part, 8, a_3, a_part2);
	int a32_len = a_part2_len; o.Double(a_part2_len, a_part2, a32);

	int det_1_len = o.Gen_Product(a11_len, a11, a22_len, a22, det_1);

	int det_2_len = o.Gen_Product(a12_len, a12, a12_len, a12, det_2);
	o.Gen_Invert(det_2_len, det_2);

	int det_len = o.Gen_Sum(det_1_len, det_1, det_2_len, det_2, det);

	if (det_len == 0 || det[det_len - 1] == 0) return 0;
	int s = (det[det_len - 1] > 0) ? (1) : (-1);

	int r_1_len = o.Gen_Product(a12_len, a12, l1_len, l1, r_1);
	o.Gen_Invert(r_1_len, r_1);
	int r_2_len = o.Gen_Product(a11_len, a11, l2_len, l2, r_2);
	int r_3_len = o.Gen_Product(a12_len, a12, l2_len, l2, r_3);
	o.Gen_Invert(r_3_len, r_3);
	int r_4_len = o.Gen_Product(a22_len, a22, l1_len, l1, r_4);

	int rr_1_len = o.Gen_Sum(r_4_len, r_4, r_3_len, r_3, rr_1);
	int rr_2_len = o.Gen_Sum(r_2_len, r_2, r_1_len, r_1, rr_2);

	double *fr_1, *fr_2, *fr_3, *sr, *r;
	int fr_1_len = o.Gen_Product_With_Alloc(det_len, det, l3_len, l3, &fr_1);
	int fr_2_len = o.Gen_Product_With_Alloc(a31_len, a31, rr_1_len, rr_1, &fr_2);
	int fr_3_len = o.Gen_Product_With_Alloc(a32_len, a32, rr_2_len, rr_2, &fr_3);
	int sr_len = o.Gen_Sum_With_Alloc(fr_3_len, fr_3, fr_2_len, fr_2, &sr);
	o.Gen_Invert(sr_len, sr);
	int r_len = o.Gen_Sum_With_Alloc(fr_1_len, fr_1, sr_len, sr, &r);
	if (r_len == 0 || r[r_len - 1] == 0) s = 0;
	else if (r[r_len - 1] < 0) s = -s;

	free(r);
	free(sr);
	free(fr_1);
	free(fr_2);
	free(fr_3);

	return -s;
}


double TMesh::incircle3D(double *p0, double *p1, double *p2, double *p3)
{
	double result = incircle3D_filtered(p0, p1, p2, p3);
	if (result) return result;
	return incircle3D_exact(p0, p1, p2, p3);
}

} //namespace T_MESH
