/****************************************************************************
* point.cpp                                                                 *
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

#include "point.h"
#include <stdlib.h>
#include <limits.h>
#include <errno.h>

namespace T_MESH
{

const Point3c INFINITE_POINT(TMESH_INFINITY, TMESH_INFINITY, TMESH_INFINITY);

//////// Lexicographic Point3c comparison //////////

// This can be used with std::sort()
bool Point3c::operator<(const Point3c& s) const
{
	if (x<s.x) return true; else if (x>s.x) return false;
	if (y<s.y) return true; else if (y>s.y) return false;
	if (z<s.z) return true; else return false;
}

/////////// Solution of a linear system 3 x 3    //////////
///// System Ax = d, where A = (a,b,c) rows, d = this /////

Point3c Point3c::linearSystem(const Point3c& a, const Point3c& b, const Point3c& c) const
{
	// Optimize - need to avoid redoing same multiplications
	coord det_A = TMESH_DETERMINANT3X3(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
	if (det_A == 0.0) return INFINITE_POINT;
	return Point3c(
		TMESH_DETERMINANT3X3(x, a.y, a.z, y, b.y, b.z, z, c.y, c.z) / det_A,
		TMESH_DETERMINANT3X3(a.x, x, a.z, b.x, y, b.z, c.x, z, c.z) / det_A,
		TMESH_DETERMINANT3X3(a.x, a.y, x, b.x, b.y, y, c.x, c.y, z) / det_A
		);
}

////////////// Projection on the line passing through A and B ///////////

Point3c Point3c::projection(const Point3c *A, const Point3c *B) const
{
	Point3c BA((*B) - (*A));
	coord l = BA.squaredLength();
	if (l == 0.0) return INFINITE_POINT;

	return ((*A) + (BA*((BA*((*this) - (*A))) / (l))));
}


////////////// Projection on the plane passing through A, B and C ///////////

Point3c Point3c::projection(const Point3c *A, const Point3c *B, const Point3c *C) const
{
	Point3c BA((*B) - (*A)), CA((*C) - (*A));
	Point3c plane_vec(BA&CA);
	Point3c p2((*this) + (plane_vec));
	return Point3c::linePlaneIntersection(*this, p2, *A, *B, *C);
}

////////////// Alignment check /////////////

bool Point3c::exactMisalignment(const Point3c *A, const Point3c *B) const
{
	if (orient2D(x, y, A->x, A->y, B->x, B->y) != 0) return true;
	if (orient2D(y, z, A->y, A->z, B->y, B->z) != 0) return true;
	if (orient2D(z, x, A->z, A->x, B->z, B->x) != 0) return true;

	return false;
}

bool Point3c::exactSameSideOnPlane(const Point3c *Q, const Point3c *A, const Point3c *B) const
{
		char o1, o2;

		o1 = orient2D(x, y, A->x, A->y, B->x, B->y);
		o2 = orient2D(Q->x, Q->y, A->x, A->y, B->x, B->y);
		if (o1 != o2) return false;

		o1 = orient2D(y, z, A->y, A->z, B->y, B->z);
		o2 = orient2D(Q->y, Q->z, A->y, A->z, B->y, B->z);
		if (o1 != o2) return false;

		o1 = orient2D(z, x, A->z, A->x, B->z, B->x);
		o2 = orient2D(Q->z, Q->x, A->z, A->x, B->z, B->x);
		if (o1 != o2) return false;

		return true;
}


//////////////////////////////////////////////////////////////////
//
// Basic predicates of type 'pointIn'
//
//////////////////////////////////////////////////////////////////

// Returns true if 'p' is a point of the segment v1-v2 (endpoints excluded)
bool Point3c::pointInInnerSegment(const Point3c *p, const Point3c *v1, const Point3c *v2)
{
	if (!p->exactMisalignment(v1, v2)) // Segment and point aligned
	{
		if (v1->x < v2->x && v1->x < p->x && p->x < v2->x) return true;
		if (v1->y < v2->y && v1->y < p->y && p->y < v2->y) return true;
		if (v1->z < v2->z && v1->z < p->z && p->z < v2->z) return true;
		if (v1->x > v2->x && v1->x > p->x && p->x > v2->x) return true;
		if (v1->y > v2->y && v1->y > p->y && p->y > v2->y) return true;
		if (v1->z > v2->z && v1->z > p->z && p->z > v2->z) return true;
	}
	return false;
}

// Returns true if 'p' is a point of the segment v1-v2 (endpoints included)
bool Point3c::pointInSegment(const Point3c *p, const Point3c *v1, const Point3c *v2)
{
	return ((*p) == (*(v1)) || (*p) == (*(v2)) || Point3c::pointInInnerSegment(p, v1, v2));
}

// Returns true if the coplanar point 'p' is in the inner area of 't'.
// Undetermined if p and t are not coplanar.
bool Point3c::pointInInnerTriangle(const Point3c *p, const Point3c *v1, const Point3c *v2, const Point3c *v3)
{
	char o1, o2, oo2, oo4, oo6;

	o1 = orient2D(p->x, p->y, v2->x, v2->y, v3->x, v3->y);
	o2 = oo2 = orient2D(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y);
	if (o1 != o2) return false;

	o1 = orient2D(p->y, p->z, v2->y, v2->z, v3->y, v3->z);
	o2 = oo4 = orient2D(v1->y, v1->z, v2->y, v2->z, v3->y, v3->z);
	if (o1 != o2) return false;

	o1 = orient2D(p->z, p->x, v2->z, v2->x, v3->z, v3->x);
	o2 = oo6 = orient2D(v1->z, v1->x, v2->z, v2->x, v3->z, v3->x);
	if (o1 != o2) return false;

	o1 = orient2D(p->x, p->y, v3->x, v3->y, v1->x, v1->y);
	o2 = oo2;
	if (o1 != o2) return false;

	o1 = orient2D(p->y, p->z, v3->y, v3->z, v1->y, v1->z);
	o2 = oo4;
	if (o1 != o2) return false;

	o1 = orient2D(p->z, p->x, v3->z, v3->x, v1->z, v1->x);
	o2 = oo6;
	if (o1 != o2) return false;

	o1 = orient2D(p->x, p->y, v1->x, v1->y, v2->x, v2->y);
	o2 = oo2;
	if (o1 != o2) return false;

	o1 = orient2D(p->y, p->z, v1->y, v1->z, v2->y, v2->z);
	o2 = oo4;
	if (o1 != o2) return false;

	o1 = orient2D(p->z, p->x, v1->z, v1->x, v2->z, v2->x);
	o2 = oo6;
	if (o1 != o2) return false;

	return true;
}

// Returns true if the coplanar point 'p' is either in the inner area of
// 't' or on its border. Undetermined if p and t are not coplanar.
bool Point3c::pointInTriangle(const Point3c *p, const Point3c *v1, const Point3c *v2, const Point3c *v3)
{
	if (Point3c::pointInSegment(p, v1, v2)) return true;
	else if (Point3c::pointInSegment(p, v2, v3)) return true;
	else if (Point3c::pointInSegment(p, v3, v1)) return true;
	else return Point3c::pointInInnerTriangle(p, v1, v2, v3);
}


//////////////////////////////////////////////////////////////////
//
// Basic predicates of type 'segmentIntersects'
//
//////////////////////////////////////////////////////////////////

// Returns true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
// Collinear overlapping segments are not considered to be properly intersecting.
bool Point3c::innerSegmentsCross(const Point3c& p1, const Point3c& p2, const Point3c& sp1, const Point3c& sp2)
{
	if (p1 == sp1 || p1 == sp2 || p2 == sp1 || p2 == sp2) return false;	// Endpoints cannot coincide
	if (p1.exactOrientation(&p2, &sp1, &sp2) != 0) return false;	// Must be coplanar
	if (p1.exactSameSideOnPlane(&p2, &sp1, &sp2) || sp1.exactSameSideOnPlane(&sp2, &p1, &p2)) return false; 	// Cannot be either on one side of the other

	if (orient2D(p1.x, p1.y, p2.x, p2.y, sp1.x, sp1.y) != 0) return true;
	if (orient2D(sp2.x, sp2.y, p2.x, p2.y, sp1.x, sp1.y) != 0) return true;
	if (orient2D(p1.y, p1.z, p2.y, p2.z, sp1.y, sp1.z) != 0) return true;
	if (orient2D(sp2.y, sp2.z, p2.y, p2.z, sp1.y, sp1.z) != 0) return true;
	if (orient2D(p1.z, p1.x, p2.z, p2.x, sp1.z, sp1.x) != 0) return true;
	if (orient2D(sp2.z, sp2.x, p2.z, p2.x, sp1.z, sp1.x) != 0) return true;
	return false;
}

// true if (p1-p2) properly intersects (sp1-sp2) at any point (endpoints included).
// Collinear overlapping segments are not considered to be properly intersecting.
bool Point3c::segmentsIntersect(const Point3c *p1, const Point3c *p2, const Point3c *sp1, const Point3c *sp2)
{
	return (p1->exactOrientation(p2, sp1, sp2) == 0 && !p1->exactSameSideOnPlane(p2, sp1, sp2) && !sp1->exactSameSideOnPlane(sp2, p1, p2));
}

bool Point3c::segmentProperlyIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3)
{
	char o1, o2, o3;

	coord mx = MIN(s1->x, s2->x);
	if (v1->x < mx && v2->x < mx && v3->x < mx) return false;
	mx = MAX(s1->x, s2->x);
	if (v1->x > mx && v2->x > mx && v3->x > mx) return false;
	mx = MIN(s1->y, s2->y);
	if (v1->y < mx && v2->y < mx && v3->y < mx) return false;
	mx = MAX(s1->y, s2->y);
	if (v1->y > mx && v2->y > mx && v3->y > mx) return false;
	mx = MIN(s1->z, s2->z);
	if (v1->z < mx && v2->z < mx && v3->z < mx) return false;
	mx = MAX(s1->z, s2->z);
	if (v1->z > mx && v2->z > mx && v3->z > mx) return false;

	if ((o1 = s1->exactOrientation(v1, v2, v3)) == 0) return false;
	if ((o2 = s2->exactOrientation(v1, v2, v3)) == 0) return false;
	if ((o1>0 && o2>0) || (o1<0 && o2<0)) return false;

	// Only one above and one below here...
	if ((o1 = s1->exactOrientation(s2, v1, v2)) == 0) return false;
	if ((o2 = s1->exactOrientation(s2, v2, v3)) == 0) return false;
	if ((o1>0 && o2<0) || (o1<0 && o2>0)) return false;
	if ((o3 = s1->exactOrientation(s2, v3, v1)) == 0) return false;
	if ((o1>0 && o3<0) || (o1<0 && o3>0)) return false;
	if ((o2>0 && o3<0) || (o2<0 && o3>0)) return false;
	return true;
}

bool Point3c::segmentIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3)
{
	char o1, o2, o3;

	coord mx = MIN(s1->x, s2->x);
	if (v1->x < mx && v2->x < mx && v3->x < mx) return false;
	mx = MAX(s1->x, s2->x);
	if (v1->x > mx && v2->x > mx && v3->x > mx) return false;
	mx = MIN(s1->y, s2->y);
	if (v1->y < mx && v2->y < mx && v3->y < mx) return false;
	mx = MAX(s1->y, s2->y);
	if (v1->y > mx && v2->y > mx && v3->y > mx) return false;
	mx = MIN(s1->z, s2->z);
	if (v1->z < mx && v2->z < mx && v3->z < mx) return false;
	mx = MAX(s1->z, s2->z);
	if (v1->z > mx && v2->z > mx && v3->z > mx) return false;

	o1 = s1->exactOrientation(v1, v2, v3);
	o2 = s2->exactOrientation(v1, v2, v3);
	if (o1 == 0 && o2 == 0)
	{
		if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2)) return true;
		if (Point3c::pointInInnerTriangle(s1, v1, v2, v3) && Point3c::pointInInnerTriangle(s2, v1, v2, v3)) return true;
		return false;
	}

	if ((o1>0 && o2>0) || (o1<0 && o2<0)) return false; // s1 and s2 are both above/below v1,v2,v3
	o1 = s1->exactOrientation(s2, v1, v2);
	o2 = s1->exactOrientation(s2, v2, v3);
	if ((o1>0 && o2<0) || (o1<0 && o2>0)) return false;
	o3 = s1->exactOrientation(s2, v3, v1);
	if ((o1>0 && o3<0) || (o1<0 && o3>0)) return false;
	if ((o2>0 && o3<0) || (o2<0 && o3>0)) return false;
	return true;
}

bool Point3c::segmentIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3, const coord& oo1, const coord& oo2)
{
	// In this case the fast reject by bounding box appears to be a disadvantage ...
	if (oo1 == 0 && oo2 == 0)
	{
		if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2)) return true;
		if (Point3c::pointInInnerTriangle(s1, v1, v2, v3) && Point3c::pointInInnerTriangle(s2, v1, v2, v3)) return true;
		return false;
	}

	if ((oo1>0 && oo2>0) || (oo1<0 && oo2<0)) return false; // s1 and s2 are both above/below v1,v2,v3
	char o1, o2, o3;
	o1 = s1->exactOrientation(s2, v1, v2);
	o2 = s1->exactOrientation(s2, v2, v3);
	if ((o1>0 && o2<0) || (o1<0 && o2>0)) return false;
	o3 = s1->exactOrientation(s2, v3, v1);
	if ((o1>0 && o3<0) || (o1<0 && o3>0)) return false;
	if ((o2>0 && o3<0) || (o2<0 && o3>0)) return false;
	return true;
}



// Returns the point of intersection between the two lines defined by (p,q) and (r,s) respectively
// Return INFINITE_POINT is lines do not intersect or if p==q or r==s
Point3c Point3c::lineLineIntersection(const Point3c& p, const Point3c& q, const Point3c& r, const Point3c& s)
{
	Point3c da(q - p);
	Point3c db(s - r);
	Point3c dc(r - p);
	Point3c dab(da&db);

	if (dab.isNull()) return INFINITE_POINT; // parallel
	if ((dc*dab) != 0.0) return INFINITE_POINT; // skew

	return p + (da*(((dc&db)*dab) / (dab*dab)));
}

// Returns the point of intersection between the line for (p,q) and the plane for (r,s,t)
// Returns INFINITE_POINT in case of parallelism
Point3c Point3c::linePlaneIntersection(const Point3c& p, const Point3c& q, const Point3c& r, const Point3c& s, const Point3c& t)
{
	coord a11(p.x - q.x), a12(p.y - q.y), a13(p.z - q.z);
	coord a21(s.x - r.x), a22(s.y - r.y), a23(s.z - r.z);
	coord a31(t.x - r.x), a32(t.y - r.y), a33(t.z - r.z);
	coord a2233(a22*a33 - a23*a32);
	coord a2133(a21*a33 - a23*a31);
	coord a2132(a21*a32 - a22*a31);
	coord den(a11*a2233 - a12*a2133 + a13*a2132);
	if (TMESH_IS_ZERO(den)) return INFINITE_POINT;
	coord num(((p.y - r.y)*(a2133) - (p.x - r.x)*(a2233) - (p.z - r.z)*(a2132))/den);
	return Point3c(p.x + a11*num, p.y + a12*num, p.z + a13*num);
}

// Returns the point of intersection between the line for (v1,v2) and the plane for 'v0' with directional vector 'd'
// Returns INFINITE_POINT in case of parallelism
Point3c Point3c::linePlaneIntersection(const Point3c& v1, const Point3c& v2, const Point3c& v0, const Point3c& d)
{
	Point3c v21(v2 - v1);
	coord den = d*v21;
	if (den == 0) return INFINITE_POINT;
	else return v1 + (v21*((d*(v0 - v1)) / den));
}


coord Point3c::squaredDistanceFromLine(const Point3c *x1, const Point3c *x2) const
{
	Point3c x21((*x2) - (*x1));
	if (x21.isNull()) return TMESH_INFINITY;
	Point3c x10((*x1) - (*this));
	x10 = x21&x10;
	return (x10*x10) / (x21*x21);
}

coord Point3c::squaredDistanceFromPlane(const Point3c& dirver, const Point3c& app_point) const
{
	coord CA2 = dirver*dirver;

	if (CA2 == 0) return TMESH_INFINITY;
	coord d = (dirver*(*this)) - (dirver*(app_point));

	return (d*d) / CA2;
}


//// Computes the closest points of the two lines 'this'-v1 and p1-p2  ////
//// Returns FALSE if the lines are parallel.                          ////

bool Point3c::closestPoints(const Point3c *v1, const Point3c *p1, const Point3c *p2, Point3c *ptOnThis, Point3c *ptOnLine2) const
{
	Point3c u = (*v1) - (*this);
	if (u.isNull()) return false;
	Point3c v = (*p2) - (*p1);
	if (v.isNull()) return false;

	coord A = u*u;
	coord B = u*v;
	coord C = v*v;
	coord denom = (A*C) - (B*B);
	if (denom == 0.0) return false; // Lines are parallel

	Point3c w0 = (*this) - (*p1);
	coord D = u*w0;
	coord E = v*w0;
	coord s = (B*E - C*D) / denom;
	coord t = (A*E - B*D) / denom;

	*ptOnThis = (*this) + (u*s);
	*ptOnLine2 = (*p1) + (v*t);

	return true;
}





//////////////// Normalization /////////////////////////

Point3c& Point3c::normalize()
{
 coord l = length();
 if (l == 0) TMesh::error("normalize : Trying to normalize a null vector !\n");
 operator/=(l);
 return *this;
}


//////////////////// Point3c rotation ////////////////////
/////////// 'ang' radians CCW around 'axis' ////////////

void Point3c::rotate(const Point3c& a, const double& ang)
{
 double l, q[4], m[3][3];
 if ((l = a.length())==0.0) return;
 l = sin(ang/2.0)/l;

 q[0] = TMESH_TO_DOUBLE(a.x)*l;
 q[1] = TMESH_TO_DOUBLE(a.y)*l;
 q[2] = TMESH_TO_DOUBLE(a.z)*l;
 q[3] = cos(ang/2.0);

 m[0][0] = 1.0 - (q[1]*q[1] + q[2]*q[2])*2.0;
 m[0][1] = (q[0] * q[1] + q[2] * q[3])*2.0;
 m[0][2] = (q[2] * q[0] - q[1] * q[3])*2.0;

 m[1][0] = (q[0] * q[1] - q[2] * q[3])*2.0;
 m[1][1] = 1.0 - (q[2] * q[2] + q[0] * q[0])*2.0;
 m[1][2] = (q[1] * q[2] + q[0] * q[3])*2.0;

 m[2][0] = (q[2] * q[0] + q[1] * q[3])*2.0;
 m[2][1] = (q[1] * q[2] - q[0] * q[3])*2.0;
 m[2][2] = 1.0 - (q[1] * q[1] + q[0] * q[0])*2.0;

 q[0] = TMESH_TO_DOUBLE(x); q[1] = TMESH_TO_DOUBLE(y); q[2] = TMESH_TO_DOUBLE(z);
 x = m[0][0]*q[0] + m[1][0]*q[1] + m[2][0]*q[2];
 y = m[0][1]*q[0] + m[1][1]*q[1] + m[2][1]*q[2];
 z = m[0][2]*q[0] + m[1][2]*q[1] + m[2][2]*q[2];
}


/////////// Distance from the line passing through A and B ////////

double Point3c::distanceFromLine(const Point3c *A, const Point3c *B) const
{
	return sqrt(TMESH_TO_DOUBLE(squaredDistanceFromLine(A, B)));
}


/////////////////// Distance from a line ///////////////////////
//// 'cc' is initialized as the point of the line whose     ////
//// distance from 'this' is minimum.                       ////

double Point3c::distanceFromLine(const Point3c *A, const Point3c *B, Point3c *cc) const
{
 Point3c AP((*A)-(*this));
 if (AP.isNull()) { cc->setValue(A); return 0.0; }
 
 Point3c BP((*B) - (*this));
 if (BP.isNull()) { cc->setValue(B); return 0.0; }

 Point3c AB((*A) - (*B));
 coord t = (AB*AB);
 if (t == 0.0) TMESH_INFINITY;
 else t = (AP*AB)/(-t);

 cc->setValue((AB*t) + A);
 return distanceFromLine(A,B);
}


////////////// Distance from a segment /////////////////

double Point3c::distanceFromEdge(const Point3c *A, const Point3c *B) const
{
	Point3c AP((*A) - (*this));
	Point3c BP((*B) - (*this));
	Point3c AB((*A) - (*B));

	if (AB*AP <= 0) return AP.length();
	if (AB*BP >= 0) return BP.length();

	return distanceFromLine(A,B);
}

/////////////////// Distance from a segment ///////////////////////
//// 'cc' is initialized as the point of the segment whose     ////
//// distance from 'this' is minimum.                          ////

double Point3c::distanceFromEdge(const Point3c *A, const Point3c *B, Point3c *cc) const
{
	Point3c AP((*A) - (*this));
	Point3c BP((*B) - (*this));
	Point3c AB((*A) - (*B));

	if (AB*AP <= 0) { cc->setValue(A); return AP.length(); }
	if (AB*BP >= 0) { cc->setValue(B); return BP.length(); }

	return distanceFromLine(A, B, cc);
}

/////////// Distance of two straight lines ///////////////

double Point3c::distanceLineLine(const Point3c *A, const Point3c *A1, const Point3c *B1) const
{
	Point3c uu1 = ((*this) - (*A))&((*A1) - (*B1));
	coord nom = ((*A) - (*A1))*(uu1);
	return sqrt(TMESH_TO_DOUBLE((nom*nom) / uu1.squaredLength()));
}

///////////////// Angle between two vectors ///////////////

double Point3c::getAngle(const Point3c& p) const
{
	return atan2(((*this)&p).length(), TMESH_TO_DOUBLE(((*this)*p)));
}






// This can be used with jqsort
int xyzCompare(const void *a, const void *b)
{
	if (((((Point3c *)a)->x < ((Point3c *)b)->x))) return -1;
	if (((((Point3c *)a)->x >((Point3c *)b)->x))) return 1;
	if (((((Point3c *)a)->y < ((Point3c *)b)->y))) return -1;
	if (((((Point3c *)a)->y >((Point3c *)b)->y))) return 1;
	if (((((Point3c *)a)->z < ((Point3c *)b)->z))) return -1;
	if (((((Point3c *)a)->z >((Point3c *)b)->z))) return 1;

	return 0;
}

} //namespace T_MESH
