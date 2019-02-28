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

const Point INFINITE_POINT(TMESH_INFINITY, TMESH_INFINITY, TMESH_INFINITY);

//////// Lexicographic Point comparison //////////

// This can be used with std::sort()
bool Point::operator<(const Point& s) const
{
	if (x<s.x) return true; else if (x>s.x) return false;
	if (y<s.y) return true; else if (y>s.y) return false;
	if (z<s.z) return true; else return false;
}

/////////// Solution of a linear system 3 x 3    //////////
///// System Ax = d, where A = (a,b,c) rows, d = this /////

Point Point::linearSystem(const Point& a, const Point& b, const Point& c) const
{
	coord det_A = TMESH_DETERMINANT3X3(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
	if (det_A == 0.0) return INFINITE_POINT;
	return Point(
		TMESH_DETERMINANT3X3(x, a.y, a.z, y, b.y, b.z, z, c.y, c.z) / det_A,
		TMESH_DETERMINANT3X3(a.x, x, a.z, b.x, y, b.z, c.x, z, c.z) / det_A,
		TMESH_DETERMINANT3X3(a.x, a.y, x, b.x, b.y, y, c.x, c.y, z) / det_A
		);
}

////////////// Projection on the line passing through A and B ///////////

Point Point::projection(const Point *A, const Point *B) const
{
	Point BA((*B) - (*A));
	coord l = BA*BA;
	if (l == 0.0) return INFINITE_POINT;

	return ((*A) + (BA*((BA*((*this) - (*A))) / (l))));
}


////////////// Projection on the plane passing through A, B and C ///////////

Point Point::projection(const Point *A, const Point *B, const Point *C) const
{
	Point BA((*B) - (*A)), CA((*C) - (*A));
	Point plane_vec(BA&CA);
	Point p2((*this) + (plane_vec));
	return Point::linePlaneIntersection(*this, p2, *A, *B, *C);
}

////////////// Alignment check /////////////

bool Point::exactMisalignment(const Point *A, const Point *B) const
{
	if (orient2D(x, y, A->x, A->y, B->x, B->y) != 0) return true;
	if (orient2D(y, z, A->y, A->z, B->y, B->z) != 0) return true;
	if (orient2D(z, x, A->z, A->x, B->z, B->x) != 0) return true;

	return false;
}

bool Point::exactSameSideOnPlane(const Point *Q, const Point *A, const Point *B) const
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
bool Point::pointInInnerSegment(const Point *p, const Point *v1, const Point *v2)
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
bool Point::pointInSegment(const Point *p, const Point *v1, const Point *v2)
{
	return ((*p) == (*(v1)) || (*p) == (*(v2)) || Point::pointInInnerSegment(p, v1, v2));
}

// Returns true if the coplanar point 'p' is in the inner area of 't'.
// Undetermined if p and t are not coplanar.
bool Point::pointInInnerTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3)
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
bool Point::pointInTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3)
{
	if (Point::pointInSegment(p, v1, v2)) return true;
	else if (Point::pointInSegment(p, v2, v3)) return true;
	else if (Point::pointInSegment(p, v3, v1)) return true;
	else return Point::pointInInnerTriangle(p, v1, v2, v3);
}


//////////////////////////////////////////////////////////////////
//
// Basic predicates of type 'segmentIntersects'
//
//////////////////////////////////////////////////////////////////

// Returns true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
// Collinear overlapping segments are not considered to be properly intersecting.
bool Point::innerSegmentsCross(const Point& p1, const Point& p2, const Point& sp1, const Point& sp2)
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
bool Point::segmentsIntersect(const Point *p1, const Point *p2, const Point *sp1, const Point *sp2)
{
	return (p1->exactOrientation(p2, sp1, sp2) == 0 && !p1->exactSameSideOnPlane(p2, sp1, sp2) && !sp1->exactSameSideOnPlane(sp2, p1, p2));
}

bool Point::segmentProperlyIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3)
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

bool Point::segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3)
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
		if (Point::pointInInnerTriangle(s1, v1, v2, v3) && Point::pointInInnerTriangle(s2, v1, v2, v3)) return true;
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

bool Point::segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3, const coord& oo1, const coord& oo2)
{
	// In this case the fast reject by bounding box appears to be a disadvantage ...
	if (oo1 == 0 && oo2 == 0)
	{
		if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2)) return true;
		if (Point::pointInInnerTriangle(s1, v1, v2, v3) && Point::pointInInnerTriangle(s2, v1, v2, v3)) return true;
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
Point Point::lineLineIntersection(const Point& p, const Point& q, const Point& r, const Point& s)
{
	Point da(q - p);
	Point db(s - r);
	Point dc(r - p);
	Point dab(da&db);

	if ((dc*dab) != 0.0) return INFINITE_POINT;

	coord k = (((dc&db)*dab) / (dab*dab));
	return p + (da*k);
}

// Returns the point of intersection between the line for (p,q) and the plane for (r,s,t)
// Returns INFINITE_POINT in case of parallelism
Point Point::linePlaneIntersection(const Point& p, const Point& q, const Point& r, const Point& s, const Point& t)
{
	coord den = TMESH_DETERMINANT3X3(p.x - q.x, p.y - q.y, p.z - q.z, s.x - r.x, s.y - r.y, s.z - r.z, t.x - r.x, t.y - r.y, t.z - r.z);
	if (den == 0) return INFINITE_POINT;
	coord num = TMESH_DETERMINANT3X3(p.x - r.x, p.y - r.y, p.z - r.z, s.x - r.x, s.y - r.y, s.z - r.z, t.x - r.x, t.y - r.y, t.z - r.z);
	coord gamma = num / den;
	return p + ((q - p)*gamma);
}

// Returns the point of intersection between the line for (v1,v2) and the plane for 'v0' with directional vector 'd'
// Returns INFINITE_POINT in case of parallelism
Point Point::linePlaneIntersection(const Point& v1, const Point& v2, const Point& v0, const Point& d)
{
	Point v21(v2 - v1);
	coord den = d*v21;
	if (den == 0) return INFINITE_POINT;
	else return v1 + (v21*((d*(v0 - v1)) / den));
}


coord Point::squaredDistanceFromLine(const Point *x1, const Point *x2) const
{
	Point x21((*x2) - (*x1));
	if (x21.isNull()) return TMESH_INFINITY;
	Point x10((*x1) - (*this));
	x10 = x21&x10;
	return (x10*x10) / (x21*x21);
}

coord Point::squaredDistanceFromPlane(const Point& dirver, const Point& app_point) const
{
	coord CA2 = dirver*dirver;

	if (CA2 == 0) return TMESH_INFINITY;
	coord d = (dirver*(*this)) - (dirver*(app_point));

	return (d*d) / CA2;
}


//// Computes the closest points of the two lines 'this'-v1 and p1-p2  ////
//// Returns FALSE if the lines are parallel.                          ////

bool Point::closestPoints(const Point *v1, const Point *p1, const Point *p2, Point *ptOnThis, Point *ptOnLine2) const
{
	Point u = (*v1) - (*this);
	if (u.isNull()) return false;
	Point v = (*p2) - (*p1);
	if (v.isNull()) return false;

	coord A = u*u;
	coord B = u*v;
	coord C = v*v;
	coord denom = (A*C) - (B*B);
	if (denom == 0.0) return false; // Lines are parallel

	Point w0 = (*this) - (*p1);
	coord D = u*w0;
	coord E = v*w0;
	coord s = (B*E - C*D) / denom;
	coord t = (A*E - B*D) / denom;

	*ptOnThis = (*this) + (u*s);
	*ptOnLine2 = (*p1) + (v*t);

	return true;
}





//////////////// Normalization /////////////////////////

Point& Point::normalize()
{
 coord l = length();
 if (l == 0) TMesh::error("normalize : Trying to normalize a null vector !\n");
 operator/=(l);
 return *this;
}


//////////////////// Point rotation ////////////////////
/////////// 'ang' radians CCW around 'axis' ////////////

void Point::rotate(const Point& a, const double& ang)
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

double Point::distanceFromLine(const Point *A, const Point *B) const
{
	return sqrt(TMESH_TO_DOUBLE(squaredDistanceFromLine(A, B)));
}


/////////////////// Distance from a line ///////////////////////
//// 'cc' is initialized as the point of the line whose     ////
//// distance from 'this' is minimum.                       ////

double Point::distanceFromLine(const Point *A, const Point *B, Point *cc) const
{
 Point AP((*A)-(*this));
 if (AP.isNull()) { cc->setValue(A); return 0.0; }
 
 Point BP((*B) - (*this));
 if (BP.isNull()) { cc->setValue(B); return 0.0; }

 Point AB((*A) - (*B));
 coord t = (AB*AB);
 if (t == 0.0) TMESH_INFINITY;
 else t = (AP*AB)/(-t);

 cc->setValue((AB*t) + A);
 return distanceFromLine(A,B);
}


////////////// Distance from a segment /////////////////

double Point::distanceFromEdge(const Point *A, const Point *B) const
{
	Point AP((*A) - (*this));
	Point BP((*B) - (*this));
	Point AB((*A) - (*B));

	if (AB*AP <= 0) return AP.length();
	if (AB*BP >= 0) return BP.length();

	return distanceFromLine(A,B);
}

/////////////////// Distance from a segment ///////////////////////
//// 'cc' is initialized as the point of the segment whose     ////
//// distance from 'this' is minimum.                          ////

double Point::distanceFromEdge(const Point *A, const Point *B, Point *cc) const
{
	Point AP((*A) - (*this));
	Point BP((*B) - (*this));
	Point AB((*A) - (*B));

	if (AB*AP <= 0) { cc->setValue(A); return AP.length(); }
	if (AB*BP >= 0) { cc->setValue(B); return BP.length(); }

	return distanceFromLine(A, B, cc);
}

/////////// Distance of two straight lines ///////////////

double Point::distanceLineLine(const Point *A, const Point *A1, const Point *B1) const
{
	Point uu1 = ((*this) - (*A))&((*A1) - (*B1));
	coord nom = ((*A) - (*A1))*(uu1);
	return sqrt(TMESH_TO_DOUBLE((nom*nom) / uu1.squaredLength()));
}

///////////////// Angle between two vectors ///////////////

double Point::getAngle(const Point& p) const
{
	return atan2(((*this)&p).length(), TMESH_TO_DOUBLE(((*this)*p)));
}






// This can be used with jqsort
int xyzCompare(const void *a, const void *b)
{
	if (((((Point *)a)->x < ((Point *)b)->x))) return -1;
	if (((((Point *)a)->x >((Point *)b)->x))) return 1;
	if (((((Point *)a)->y < ((Point *)b)->y))) return -1;
	if (((((Point *)a)->y >((Point *)b)->y))) return 1;
	if (((((Point *)a)->z < ((Point *)b)->z))) return -1;
	if (((((Point *)a)->z >((Point *)b)->z))) return 1;

	return 0;
}

} //namespace T_MESH
