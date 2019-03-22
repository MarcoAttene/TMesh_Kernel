/****************************************************************************
* mixedPredicates.h                                                         *
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

#ifndef MIXED_PREDICATES
#define MIXED_PREDICATES

#include "unprecise_numbers.h"

namespace T_MESH
{
	class unprecise_number : public interval_number
	{
	public:
		inline unprecise_number(double a) : interval_number(a) {}
		inline unprecise_number(double a, double b) : interval_number(a, b) {}

		inline unprecise_number(const interval_number& b) : interval_number(b) {}

#ifdef USE_HYBRID_KERNEL
#ifdef USE_LAZY_KERNEL
		inline unprecise_number(const PM_Rational& a)
		{
			if (a.isOfRationalType()) init(a.getVal().unprecise());
			else init(a.getDVal());
		}
#else
		inline unprecise_number(const PM_Rational& a)
		{
			if (a.isOfRationalType()) init(a.getVal());
			else init(a.getDVal());
		}
#endif
#endif

	};

	// Orient 2D
	template <class T>
	inline T unpreciseOrient2d_T(T& ax, T& ay, T& bx, T& by, T& cx, T& cy)
	{
		return (((ax - cx)*(by - cy)) - ((ay - cy)*(bx - cx)));
	}

#ifdef USE_HYBRID_KERNEL
	template <class T>
	inline T unpreciseOrient2d(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry)
	{
		T ax(px), ay(py), bx(qx), by(qy), cx(rx), cy(ry);
		return unpreciseOrient2d_T<T>(ax, ay, bx, by, cx, cy);
	}

#endif

	template <class T>
	inline T unpreciseOrient2d(const double& px, const double& py, const double& qx, const double& qy, const double& rx, const double& ry)
	{
		T ax(px), ay(py), bx(qx), by(qy), cx(rx), cy(ry);
		return unpreciseOrient2d_T<T>(ax, ay, bx, by, cx, cy);
	}

	// Orient 3D
	template <class T>
	inline T unpreciseOrient3d_T(T& a11, T& a12, T& a13, T& a21, T& a22, T& a23, T& a31, T& a32, T& a33, T& cx, T& cy, T& cz)
	{
		a11 -= cx; a12 -= cy; a13 -= cz;
		a21 -= cx; a22 -= cy; a23 -= cz;
		a31 -= cx; a32 -= cy; a33 -= cz;
		return ((a11)*((a22)*(a33)-(a23)*(a32)) - (a12)*((a21)*(a33)-(a23)*(a31)) + (a13)*((a21)*(a32)-(a22)*(a31)));
	}

	template <class T>
	inline T unpreciseOrient3d(const Point3c *d, const Point3c *a, const Point3c *b, const Point3c *c)
	{
		T a11(d->x), a12(d->y), a13(d->z);
		T a21(a->x), a22(a->y), a23(a->z);
		T a31(b->x), a32(b->y), a33(b->z);
		T cx(c->x), cy(c->y), cz(c->z);
		return unpreciseOrient3d_T<T>(a11, a12, a13, a21, a22, a23, a31, a32, a33, cx, cy, cz);
	}

	template <class T>
	inline T unpreciseOrient3d(const double *d, const double *a, const double *b, const double *c)
	{
		T a11((d)[0]), a12((d)[1]), a13((d)[2]);
		T a21((a)[0]), a22((a)[1]), a23((a)[2]);
		T a31((b)[0]), a32((b)[1]), a33((b)[2]);
		T cx((c)[0]), cy((c)[1]), cz((c)[2]);
		return unpreciseOrient3d_T<T>(a11, a12, a13, a21, a22, a23, a31, a32, a33, cx, cy, cz);
	}

	// Insphere
	template <class T>
	inline T unpreciseInsphere_T(T& pex, T& pey, T& pez, T& aex, T& aey, T& aez, T& bex, T& bey, T& bez, T& cex, T& cey, T& cez, T& dex, T& dey, T& dez)
	{
		aex -= pex; aey -= pey; aez -= pez;
		bex -= pex; bey -= pey; bez -= pez;
		cex -= pex; cey -= pey; cez -= pez;
		dex -= pex; dey -= pey; dez -= pez;

		T ab(aex * bey); ab -= (bex * aey);
		T bc(bex * cey); bc -= (cex * bey);
		T cd(cex * dey); cd -= (dex * cey);
		T da(dex * aey); da -= (aex * dey);
		T ac(aex * cey); ac -= (cex * aey);
		T bd(bex * dey); bd -= (dex * bey);

		T abc(aez * bc - bez * ac + cez * ab);
		T bcd(bez * cd - cez * bd + dez * bc);
		T cda(cez * da + dez * ac + aez * cd);
		T dab(dez * ab + aez * bd + bez * da);

		T alift(aex * aex + aey * aey + aez * aez);
		T blift(bex * bex + bey * bey + bez * bez);
		T clift(cex * cex + cey * cey + cez * cez);
		T dlift(dex * dex + dey * dey + dez * dez);

		return ((dlift * abc - clift * dab) + (blift * cda - alift * bcd));
	}

	template <class T>
	inline T unpreciseInsphere(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd, const Point3c *pe)
	{
		T pex(pe->x), pey(pe->y), pez(pe->z);
		T aex(pa->x), aey(pa->y), aez(pa->z);
		T bex(pb->x), bey(pb->y), bez(pb->z);
		T cex(pc->x), cey(pc->y), cez(pc->z);
		T dex(pd->x), dey(pd->y), dez(pd->z);
		return unpreciseInsphere_T<T>(pex, pey, pez, aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}

	template <class T>
	inline T unpreciseInsphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe)
	{
		T pex((pe)[0]), pey(pe[1]), pez(pe[2]);
		T aex((pa)[0]), aey(pa[1]), aez(pa[2]);
		T bex((pb)[0]), bey(pb[1]), bez(pb[2]);
		T cex((pc)[0]), cey(pc[1]), cez(pc[2]);
		T dex((pd)[0]), dey((pd)[1]), dez((pd)[2]);
		return unpreciseInsphere_T<T>(pex, pey, pez, aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}


	// InCircle3D
	template <class T>
	inline T unpreciseInCircle3D_T(T& aex, T& aey, T& aez, T& bex, T& bey, T& bez, T& cex, T& cey, T& cez, T& dex, T& dey, T& dez)
	{
		T p10x = bex - aex, p10y = bey - aey, p10z = bez - aez;
		T l1 = p10x * p10x + p10y * p10y + p10z * p10z;
		T p20x = cex - aex, p20y = cey - aey, p20z = cez - aez;
		T l2 = p20x * p20x + p20y * p20y + p20z * p20z;
		T p30x = dex - aex, p30y = dey - aey, p30z = dez - aez;
		T l3 = p30x * p30x + p30y * p30y + p30z * p30z;

		T a11 = l1 * T(2);
		T a12 = (p10x * p20x + p10y * p20y + p10z * p20z) * T(2);
		T a22 = l2 * T(2);
		T a31 = (p30x * p10x + p30y * p10y + p30z * p10z) * T(2);
		T a32 = (p30x * p20x + p30y * p20y + p30z * p20z) * T(2);

		T det = (a11 * a22) - (a12 * a12);
		T r = ((((a31 * ((a22 * l1) - (a12 * l2)))) + (a32 * ((a11 * l2) - (a12 * l1)))) - (det * l3));

		return det*r;
	}

	template <class T>
	inline T unpreciseInCircle3D(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd)
	{
		T aex(pa->x), aey(pa->y), aez(pa->z);
		T bex(pb->x), bey(pb->y), bez(pb->z);
		T cex(pc->x), cey(pc->y), cez(pc->z);
		T dex(pd->x), dey(pd->y), dez(pd->z);
		return unpreciseInCircle3D_T<T>(aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}

	template <class T>
	inline T unpreciseInCircle3D(const double *pa, const double *pb, const double *pc, const double *pd)
	{
		T aex((pa)[0]), aey(pa[1]), aez(pa[2]);
		T bex((pb)[0]), bey(pb[1]), bez(pb[2]);
		T cex((pc)[0]), cey(pc[1]), cez(pc[2]);
		T dex((pd)[0]), dey((pd)[1]), dez((pd)[2]);
		return unpreciseInCircle3D_T<T>(aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}

#ifdef USE_LAZY_KERNEL
	inline tmesh_fraction preciseOrient2d(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry)
	{
		tmesh_fraction qpx(px.toRational().exact()), qpy(py.toRational().exact());
		tmesh_fraction qqx(qx.toRational().exact()), qqy(qy.toRational().exact());
		tmesh_fraction qrx(rx.toRational().exact()), qry(ry.toRational().exact());
		return unpreciseOrient2d_T<tmesh_fraction>(qpx, qpy, qqx, qqy, qrx, qry);
	}

	inline tmesh_fraction preciseOrient3d(const Point3c *d, const Point3c *a, const Point3c *b, const Point3c *c)
	{
		tmesh_fraction a11(d->x.toRational().exact()), a12(d->y.toRational().exact()), a13(d->z.toRational().exact());
		tmesh_fraction a21(a->x.toRational().exact()), a22(a->y.toRational().exact()), a23(a->z.toRational().exact());
		tmesh_fraction a31(b->x.toRational().exact()), a32(b->y.toRational().exact()), a33(b->z.toRational().exact());
		tmesh_fraction cx(c->x.toRational().exact()), cy(c->y.toRational().exact()), cz(c->z.toRational().exact());
		return unpreciseOrient3d_T<tmesh_fraction>(a11, a12, a13, a21, a22, a23, a31, a32, a33, cx, cy, cz);
	}

	inline tmesh_fraction preciseInsphere(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd, const Point3c *pe)
	{
		tmesh_fraction pex(pe->x.toRational().exact()), pey(pe->y.toRational().exact()), pez(pe->z.toRational().exact());
		tmesh_fraction aex(pa->x.toRational().exact()), aey(pa->y.toRational().exact()), aez(pa->z.toRational().exact());
		tmesh_fraction bex(pb->x.toRational().exact()), bey(pb->y.toRational().exact()), bez(pb->z.toRational().exact());
		tmesh_fraction cex(pc->x.toRational().exact()), cey(pc->y.toRational().exact()), cez(pc->z.toRational().exact());
		tmesh_fraction dex(pd->x.toRational().exact()), dey(pd->y.toRational().exact()), dez(pd->z.toRational().exact());
		return unpreciseInsphere_T<tmesh_fraction>(pex, pey, pez, aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}

	inline tmesh_fraction preciseInCircle3D(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd)
	{
		tmesh_fraction aex(pa->x.toRational().exact()), aey(pa->y.toRational().exact()), aez(pa->z.toRational().exact());
		tmesh_fraction bex(pb->x.toRational().exact()), bey(pb->y.toRational().exact()), bez(pb->z.toRational().exact());
		tmesh_fraction cex(pc->x.toRational().exact()), cey(pc->y.toRational().exact()), cez(pc->z.toRational().exact());
		tmesh_fraction dex(pd->x.toRational().exact()), dey(pd->y.toRational().exact()), dez(pd->z.toRational().exact());
		return unpreciseInCircle3D_T<tmesh_fraction>(aex, aey, aez, bex, bey, bez, cex, cey, cez, dex, dey, dez);
	}
#endif
}

#endif // MIXED_PREDICATES
