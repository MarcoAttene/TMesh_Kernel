/****************************************************************************
* mixedPredicates.cpp                                                       *
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

#include <iostream>

#include "point.h"
#ifdef USE_HYBRID_KERNEL
#include "mixedPredicates.h"
#endif

namespace T_MESH
{

/****************************************************************************
*                                                                           *
*           MIXED PREDICATES                                                *
*                                                                           *
****************************************************************************/

#ifdef USE_HYBRID_KERNEL

char mixedFilteredOrient2d(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry)
{
	unprecise_number det = unpreciseOrient2d<unprecise_number>(px, py, qx, qy, rx, ry);
	if (det.signIsReliable()) return det.sign();
	else
	{
// Here we already know that lazy evaluation fails. Switch to exact directly
#ifdef USE_LAZY_KERNEL
		tmesh_fraction det = preciseOrient2d(px, py, qx, qy, rx, ry);
#else
		bool wwu = PM_Rational::isUsingRationals();
		PM_Rational::useRationals(true);
		PM_Rational det = ((px - rx)*(qy - ry) - (py - ry)*(qx - rx));
		PM_Rational::useRationals(wwu);
#endif
		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
}

char mixedFilteredOrient3d(const Point *d, const Point *a, const Point *b, const Point *c)
{
	unprecise_number det = unpreciseOrient3d<unprecise_number>(d, a, b, c);
	if (det.signIsReliable()) return det.sign();
	else
	{
#ifdef USE_LAZY_KERNEL
		tmesh_fraction det = preciseOrient3d(d, a, b, c);
#else
		bool wwu = PM_Rational::isUsingRationals();
		PM_Rational::useRationals(true);
		//PM_Rational det = TMESH_DETERMINANT3X3(d->x - c->x, d->y - c->y, d->z - c->z, a->x - c->x, a->y - c->y, a->z - c->z, b->x - c->x, b->y - c->y, b->z - c->z);
		PM_Rational det = unpreciseOrient3d<PM_Rational>(d, a, b, c);
		PM_Rational::useRationals(wwu);
#endif
		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
}

char mixedFilteredInsphere(const Point *pa, const Point *pb, const Point *pc, const Point *pd, const Point *pe)
{
	unprecise_number det = unpreciseInsphere<unprecise_number>(pa, pb, pc, pd, pe);
	if (det.signIsReliable()) return det.sign();
	else
	{
#ifdef USE_LAZY_KERNEL
		tmesh_fraction det = preciseInsphere(pa, pb, pc, pd, pe);
#else
		bool wwu = PM_Rational::isUsingRationals();
		PM_Rational::useRationals(true);
		coord det = unpreciseInsphere<coord>(pa, pb, pc, pd, pe);
		PM_Rational::useRationals(wwu);
#endif
		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
}

char mixedFilteredInCircle3D(const Point *pa, const Point *pb, const Point *pc, const Point *pd)
{
	unprecise_number det = unpreciseInCircle3D<unprecise_number>(pa, pb, pc, pd);
	if (det.signIsReliable()) return det.sign();
	else
	{
#ifdef USE_LAZY_KERNEL
		tmesh_fraction det = preciseInCircle3D(pa, pb, pc, pd);
#else
		bool wwu = PM_Rational::isUsingRationals();
		PM_Rational::useRationals(true);
		coord det = unpreciseInCircle3D<coord>(pa, pb, pc, pd);
		PM_Rational::useRationals(wwu);
#endif
		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
}


char orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry)
{
	if (px.isOfDoubleType() && py.isOfDoubleType() && qx.isOfDoubleType() && qy.isOfDoubleType() && rx.isOfDoubleType() && ry.isOfDoubleType())
	{
		double pqr[6];
		pqr[0] = px.getDVal();  pqr[1] = py.getDVal();
		pqr[2] = qx.getDVal();  pqr[3] = qy.getDVal();
		pqr[4] = rx.getDVal();  pqr[5] = ry.getDVal();

		unsigned int prevround = getFPURoundingMode();
		setFPUModeToRoundNEAR();
		double det = TMesh::tri_orientation(pqr, pqr + 2, pqr + 4);
		setFPURoundingMode(prevround);

		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
	else return mixedFilteredOrient2d(px, py, qx, qy, rx, ry);
}

char orient3D(const Point *t, const Point *a, const Point *b, const Point *c)
{
	if (a->x.isOfDoubleType() && a->y.isOfDoubleType() && a->z.isOfDoubleType() &&
		t->x.isOfDoubleType() && t->y.isOfDoubleType() && t->z.isOfDoubleType() &&
		b->x.isOfDoubleType() && b->y.isOfDoubleType() && b->z.isOfDoubleType() &&
		c->x.isOfDoubleType() && c->y.isOfDoubleType() && c->z.isOfDoubleType())
	{
		double p1[3], p2[3], p3[3], p4[3];
		p1[0] = (t->x).getDVal(); p1[1] = (t->y).getDVal(); p1[2] = (t->z).getDVal();
		p2[0] = (a->x).getDVal(); p2[1] = (a->y).getDVal(); p2[2] = (a->z).getDVal();
		p3[0] = (b->x).getDVal(); p3[1] = (b->y).getDVal(); p3[2] = (b->z).getDVal();
		p4[0] = (c->x).getDVal(); p4[1] = (c->y).getDVal(); p4[2] = (c->z).getDVal();

		unsigned int prevround = getFPURoundingMode();
		setFPUModeToRoundNEAR();
		double det = TMesh::tet_orientation(p1, p2, p3, p4);
		setFPURoundingMode(prevround);

		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	}
	else return mixedFilteredOrient3d(t, a, b, c);
}

char inSphere3D(const Point *pa, const Point *pb, const Point *pc, const Point *pd, const Point *pe)
{
	if (pa->x.isOfDoubleType() && pa->y.isOfDoubleType() && pa->z.isOfDoubleType() &&
		pe->x.isOfDoubleType() && pe->y.isOfDoubleType() && pe->z.isOfDoubleType() &&
		pd->x.isOfDoubleType() && pd->y.isOfDoubleType() && pd->z.isOfDoubleType() &&
		pb->x.isOfDoubleType() && pb->y.isOfDoubleType() && pb->z.isOfDoubleType() &&
		pc->x.isOfDoubleType() && pc->y.isOfDoubleType() && pc->z.isOfDoubleType())
	{
		double p1[3], p2[3], p3[3], p4[3], p5[3];
		p1[0] = (pa->x).getDVal(); p1[1] = (pa->y).getDVal(); p1[2] = (pa->z).getDVal();
		p2[0] = (pb->x).getDVal(); p2[1] = (pb->y).getDVal(); p2[2] = (pb->z).getDVal();
		p3[0] = (pc->x).getDVal(); p3[1] = (pc->y).getDVal(); p3[2] = (pc->z).getDVal();
		p4[0] = (pd->x).getDVal(); p4[1] = (pd->y).getDVal(); p4[2] = (pd->z).getDVal();
		p5[0] = (pe->x).getDVal(); p5[1] = (pe->y).getDVal(); p5[2] = (pe->z).getDVal();

		unsigned int prevround = getFPURoundingMode();
		setFPUModeToRoundNEAR();
		double det = TMesh::insphere(p1, p2, p3, p4, p5);
		setFPURoundingMode(prevround);

		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	} else return mixedFilteredInsphere(pa, pb, pc, pd, pe);
}


char inCircle3D(const Point *pa, const Point *pb, const Point *pc, const Point *pd)
{
	if (pa->x.isOfDoubleType() && pa->y.isOfDoubleType() && pa->z.isOfDoubleType() &&
		pd->x.isOfDoubleType() && pd->y.isOfDoubleType() && pd->z.isOfDoubleType() &&
		pb->x.isOfDoubleType() && pb->y.isOfDoubleType() && pb->z.isOfDoubleType() &&
		pc->x.isOfDoubleType() && pc->y.isOfDoubleType() && pc->z.isOfDoubleType())
	{
		double p1[3], p2[3], p3[3], p4[3];
		p1[0] = (pa->x).getDVal(); p1[1] = (pa->y).getDVal(); p1[2] = (pa->z).getDVal();
		p2[0] = (pb->x).getDVal(); p2[1] = (pb->y).getDVal(); p2[2] = (pb->z).getDVal();
		p3[0] = (pc->x).getDVal(); p3[1] = (pc->y).getDVal(); p3[2] = (pc->z).getDVal();
		p4[0] = (pd->x).getDVal(); p4[1] = (pd->y).getDVal(); p4[2] = (pd->z).getDVal();

		unsigned int prevround = getFPURoundingMode();
		setFPUModeToRoundNEAR();
		double det = TMesh::incircle3D(p1, p2, p3, p4);
		setFPURoundingMode(prevround);

		return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
	} else return mixedFilteredInCircle3D(pa, pb, pc, pd);
}

#else

char orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry)
{
//	return ((px - rx)*(qy - ry) - (py - ry)*(qx - rx));

	double pqr[6];
	pqr[0] = TMESH_TO_DOUBLE(px);  pqr[1] = TMESH_TO_DOUBLE(py);
	pqr[2] = TMESH_TO_DOUBLE(qx);  pqr[3] = TMESH_TO_DOUBLE(qy);
	pqr[4] = TMESH_TO_DOUBLE(rx);  pqr[5] = TMESH_TO_DOUBLE(ry);
	PM_Rational det = TMesh::tri_orientation(pqr, pqr + 2, pqr + 4);
	return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
}

char orient3D(const Point *t, const Point *a, const Point *b, const Point *c)
{
//	return TMESH_DETERMINANT3X3(t->x - c->x, t->y - c->y, t->z - c->z, a->x - c->x, a->y - c->y, a->z - c->z, b->x - c->x, b->y - c->y, b->z - c->z);

	double p1[3], p2[3], p3[3], p4[3];
	p1[0] = TMESH_TO_DOUBLE(t->x); p1[1] = TMESH_TO_DOUBLE(t->y); p1[2] = TMESH_TO_DOUBLE(t->z);
	p2[0] = TMESH_TO_DOUBLE(a->x); p2[1] = TMESH_TO_DOUBLE(a->y); p2[2] = TMESH_TO_DOUBLE(a->z);
	p3[0] = TMESH_TO_DOUBLE(b->x); p3[1] = TMESH_TO_DOUBLE(b->y); p3[2] = TMESH_TO_DOUBLE(b->z);
	p4[0] = TMESH_TO_DOUBLE(c->x); p4[1] = TMESH_TO_DOUBLE(c->y); p4[2] = TMESH_TO_DOUBLE(c->z);
	double det =  TMesh::tet_orientation(p1, p2, p3, p4);
	return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
}

char inSphere3D(const Point *pa, const Point *pb, const Point *pc, const Point *pd, const Point *pe)
{
	double p1[3], p2[3], p3[3], p4[3], p5[3];
	p1[0] = TMESH_TO_DOUBLE(pa->x); p1[1] = TMESH_TO_DOUBLE(pa->y); p1[2] = TMESH_TO_DOUBLE(pa->z);
	p2[0] = TMESH_TO_DOUBLE(pb->x); p2[1] = TMESH_TO_DOUBLE(pb->y); p2[2] = TMESH_TO_DOUBLE(pb->z);
	p3[0] = TMESH_TO_DOUBLE(pc->x); p3[1] = TMESH_TO_DOUBLE(pc->y); p3[2] = TMESH_TO_DOUBLE(pc->z);
	p4[0] = TMESH_TO_DOUBLE(pd->x); p4[1] = TMESH_TO_DOUBLE(pd->y); p4[2] = TMESH_TO_DOUBLE(pd->z);
	p5[0] = TMESH_TO_DOUBLE(pe->x); p5[1] = TMESH_TO_DOUBLE(pe->y); p5[2] = TMESH_TO_DOUBLE(pe->z);
	double det = TMesh::insphere(p1, p2, p3, p4, p5);
	return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
}

char inCircle3D(const Point *pa, const Point *pb, const Point *pc, const Point *pd)
{
	double p1[3], p2[3], p3[3], p4[3];
	p1[0] = TMESH_TO_DOUBLE(pa->x); p1[1] = TMESH_TO_DOUBLE(pa->y); p1[2] = TMESH_TO_DOUBLE(pa->z);
	p2[0] = TMESH_TO_DOUBLE(pb->x); p2[1] = TMESH_TO_DOUBLE(pb->y); p2[2] = TMESH_TO_DOUBLE(pb->z);
	p3[0] = TMESH_TO_DOUBLE(pc->x); p3[1] = TMESH_TO_DOUBLE(pc->y); p3[2] = TMESH_TO_DOUBLE(pc->z);
	p4[0] = TMESH_TO_DOUBLE(pd->x); p4[1] = TMESH_TO_DOUBLE(pd->y); p4[2] = TMESH_TO_DOUBLE(pd->z);
	double det = TMesh::incircle3D(p1, p2, p3, p4);
	return (det > 0) ? (1) : ((det < 0) ? (-1) : (0));
}

#endif

} //namespace T_MESH
