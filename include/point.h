/****************************************************************************
* point.h                                                                   *
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

#ifndef _POINT_H
#define _POINT_H

#include "basics.h"

namespace T_MESH
{
//! Orientation predicates on PM_Rationals

class Point3c;
char orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry);
char orient3D(const Point3c *t, const Point3c *a, const Point3c *b, const Point3c *c);
char inSphere3D(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd, const Point3c *pe);
char inCircle3D(const Point3c *pa, const Point3c *pb, const Point3c *pc, const Point3c *pd);


//! Geometric point definition

//! This class represents a point in the Euclidean 3D space. It can be used
//! to  represent  3D vectors originating at (0,0,0) and terminating at the
//! corresponding point. Several methods of  this  class  are  intended  to
//! manipulate  vectors  rather  than  points;  for  example, a call of the
//! method normalize is an actual normalization if the object is a  vector,
//! but  it  has  to  be intended as a projection on the unit sphere if the
//! object is intended to be a point. An object of type Point3c is a triplet
//! of coordinates. Each coordinate is a number of type 'coord' which can be
//! often treated as a standard double. Operations on points include addition,
//! subtraction, cross and dot product, and others. This class implements
//! several useful operations using vector arithmethic. For example,
//! the simple piece of code "A = B*C;" assignes to A the value of the  dot
//! product of B and C.

class Point3c
{
 public :
 coord x,y,z;					//!< Coordinates

 //! Creates a new point with coordinates (0,0,0).
 inline Point3c() : x(0), y(0), z(0) { }

 //! Creates a new point with the same coordinates as 's'. 
 inline Point3c(const Point3c *s) : x(s->x), y(s->y), z(s->z) { }

 //! Creates a new point with the same coordinates as 's'. 
 inline Point3c(const Point3c& s) : x(s.x), y(s.y), z(s.z) { }

 //! Creates a new point with coordinates (a,b,c).
 inline Point3c(const coord& a, const coord& b, const coord& c) : x(a), y(b), z(c) { }

 //! Do not remove this. It makes the compiler produce a vtable for this object.
 TMESH_VIRTUAL bool isPoint() const { return true; }

 //! Set the coordinates to (a,b,c).
 inline void	setValue(const coord& a, const coord& b, const coord& c) { x = a; y = b; z = c; }

 //! Set the coordinates as those of 'p'
 inline void	setValue(const Point3c& p) { x = p.x; y = p.y; z = p.z; }

 //! Set the coordinates as those of '*p'
 inline void	setValue(const Point3c *p) { x = p->x; y = p->y; z = p->z; }

 //! Returns the vector difference
 inline Point3c 	operator-(const Point3c& p) const { return Point3c(x - p.x, y - p.y, z - p.z); }

 //! Returns the vector sum
 inline Point3c 	operator+(const Point3c& p) const { return Point3c(x + p.x, y + p.y, z + p.z); }

 //! Sums another point
 inline void 	operator+=(const Point3c& p) { x += p.x; y += p.y; z += p.z; }

 //! Subtracts another point
 inline void 	operator-=(const Point3c& p) { x -= p.x; y -= p.y; z -= p.z; }

 //! Returns the Cross Product
 inline Point3c 	operator&(const Point3c& p) const { return Point3c(y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x); }

 //! Returns the Dot Product
 inline coord operator*(const Point3c& p) const { return (x*p.x + y*p.y + z*p.z); }

 //! Returns the product with a scalar
 inline Point3c  operator*(const coord& d) const { return Point3c(x*d, y*d, z*d); }

 //! Multiplies by a scalar
 inline void 	operator*=(const coord& m) { x *= m; y *= m; z *= m; }

 //! Divides by a scalar
 inline void 	operator/=(const coord& m) { x /= m; y /= m; z /= m; }

 //! Returns the vector divided by the scalar
 inline Point3c 	operator/(const coord& d) const { return Point3c(x / d, y / d, z / d); }

 //! TRUE iff coordinates are equal
 inline bool  	operator==(const Point3c& p) const { return (x == p.x && y == p.y && z == p.z); }

 //! FALSE iff coordinates are equal
 inline bool  	operator!=(const Point3c& p) const { return (x != p.x || y != p.y || z != p.z); }

 //! TRUE iff this is lexycographically smaller than s
 bool operator<(const Point3c& s) const;

 //! Get the i'th coordinate
 inline coord& at(unsigned char i) { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }
 inline const coord& at(unsigned char i) const { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }
 inline coord& operator[](unsigned char i) { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }
 inline const coord& operator[](unsigned char i) const { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }

 //! Returns the inverse vector
 inline Point3c 	inverse() const { return Point3c(-x, -y, -z); }

 //! Inverts the vector
 inline void 	invert() { x = -x; y = -y; z = -z; }

 //! TRUE if vector is (0,0,0)
 inline bool  	isNull() const { return (TMESH_IS_ZERO(x) && TMESH_IS_ZERO(y) && TMESH_IS_ZERO(z)); }

 //! Returns the solution of the linear system Ax = d, where A is a 3x3 matrix whose rows are row1, row2 and row3, d = this
 Point3c  linearSystem(const Point3c& row1, const Point3c& row2, const Point3c& row3) const;

 //! Projects the vector on the plane with normal 'n' passing through the origin.
 inline void   project(const Point3c *n) { setValue((*this) - ((*n)*((*this)*(*n)))); }

 //! Returns the projection of the point on the straight line though 'a' and 'b'.
 Point3c  projection(const Point3c *a, const Point3c *b) const;

 //! Returns the projection of the point on the plane though 'a', 'b' and 'c'.
 Point3c  projection(const Point3c *a, const Point3c *b, const Point3c *c) const;

 //! Exact orientation test.
 //! Return value is positive iff the tetrahedron (this,a,b,c) has a positive volume;
 //! It is negative iff the tetrahedron (this,a,b,c) has a negative volume;
 //! It is zero iff the tetrahedron (this,a,b,c) has a zero volume.
 inline char exactOrientation(const Point3c *a, const Point3c *b, const Point3c *c) const { return orient3D(this, a, b, c); }
 inline char exactOrientation(const Point3c& a, const Point3c& b, const Point3c& c) const { return orient3D(this, &a, &b, &c); }

 //! Return value is positive iff this point belongs to the interior of the sphere
 //! by a,b,c,d. Negative if outside. Zero if on the sphere's surface.
 //! Assumes that a->exactOrientation(b, c, d) is positive. Otherwise the sign is flipped.
 inline char inSphere(const Point3c *a, const Point3c *b, const Point3c *c, const Point3c *d) const { return inSphere3D(a, b, c, d, this); }
 inline char inSphere(const Point3c& a, const Point3c& b, const Point3c& c, const Point3c& d) const { return inSphere3D(&a, &b, &c, &d, this); }

 //! Return value is positive iff this point belongs to the interior of the circle
 //! by a,b,c. Negative if outside. Zero if on the circle.
 //! Assumes that this and a,b,c are coplanar. Result is undetermined otherwise.
 inline char incircle3D(const Point3c *a, const Point3c *b, const Point3c *c) const { return inCircle3D(a, b, c, this); }
 inline char incircle3D(const Point3c& a, const Point3c& b, const Point3c& c) const { return inCircle3D(&a, &b, &c, this); }

 //! Exact misalignment test. Returns TRUE iff points are not aligned.
 bool exactMisalignment(const Point3c *a, const Point3c *b) const;
 inline bool 	notAligned(const Point3c *a, const Point3c *b) const { return exactMisalignment(a, b); }

 //! Exact planar side test. Returns TRUE iff 'this', Q, A and B are coplanar
 //! and 'this' and Q are (properly) on the same side of A-B.
 //! Coplanarity is not checked, result is undetermined if
 //! 'this', Q, A and B are not coplanar.
 bool exactSameSideOnPlane(const Point3c *Q, const Point3c *A, const Point3c *B) const;

 //! true if 'p' is a point of the segment v1-v2 (endpoints excluded)
 static bool pointInInnerSegment(const Point3c *p, const Point3c *v1, const Point3c *v2);

 //! true if 'p' is a point of the segment v1-v2 (endpoints included)
 static bool pointInSegment(const Point3c *p, const Point3c *v1, const Point3c *v2);

 //! true if the coplanar point 'p' is in the inner area of v1-v2-v3.
 //! Undetermined if points are not coplanar.
 static bool pointInInnerTriangle(const Point3c *p, const Point3c *v1, const Point3c *v2, const Point3c *v3);

 //! true if the coplanar point 'p' is either in the inner area of v1-v2-v3 or on its border.
 //! Undetermined if points are not coplanar.
 static bool pointInTriangle(const Point3c *p, const Point3c *v1, const Point3c *v2, const Point3c *v3);

 //! true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
 //! Collinear overlapping segments are not considered to be properly intersecting.
 static bool innerSegmentsCross(const Point3c& p1, const Point3c& p2, const Point3c& sp1, const Point3c& sp2);

 //! true if (p1-p2) properly intersects (sp1-sp2) at any point (endpoints included).
 //! Collinear overlapping segments are not considered to be properly intersecting.
 static bool segmentsIntersect(const Point3c *p1, const Point3c *p2, const Point3c *sp1, const Point3c *sp2);

 //! true if inner segment (s1-s2) intersects the triangle v1-v2-v3 (border excluded) at a single point
 static bool segmentProperlyIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3);

 //! true if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
 static bool segmentIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3);

 //! true if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
 //! Accelerated version - relative orientations are passed as parameters.
 static bool segmentIntersectsTriangle(const Point3c *s1, const Point3c *s2, const Point3c *v1, const Point3c *v2, const Point3c *v3, const coord& o1, const coord& o2);

 //! Itersection point between lines p-q and r-s. Return INFINITE_POINT if lines are either non-intersecting or degenerate.
 static Point3c lineLineIntersection(const Point3c& p, const Point3c& q, const Point3c& r, const Point3c& s);

 //! Itersection point between line p-q and plane r-s-t. Return INFINITE_POINT for parallel/degenerate args.
 static Point3c linePlaneIntersection(const Point3c& p, const Point3c& q, const Point3c& r, const Point3c& s, const Point3c& t);

 //! Itersection point between line p-q and plane for 'v0' with directional vector 'd'. Return INFINITE_POINT for parallel/degenerate args.
 static Point3c linePlaneIntersection(const Point3c& p, const Point3c& q, const Point3c& v0, const Point3c& d);


 //! Squared distance from origin
 inline coord squaredLength() const { return (x*x + y*y + z*z); }

 //! Squared distance from '*b'
 inline coord squaredDistance(const Point3c *b) const { return (((*(this)) - (*b)).squaredLength()); }

 //! Squared area of the triangle p-q-r.
 inline static coord squaredTriangleArea3D(const Point3c& p, const Point3c& q, const Point3c& r) { return ((p - r)&(q - r)).squaredLength()*0.25; }

 //! Squared distance from straight line through 'a' and 'b'
 coord squaredDistanceFromLine(const Point3c *a, const Point3c *b) const;

 //! Squared distance from plane through 'app_point' and 'having directional vector 'dirver'
 coord squaredDistanceFromPlane(const Point3c& dirver, const Point3c& app_point) const;

 //! Line-line closest point computation.
 //! Computes the closest points of the line passing through this and this2,
 //! and the line passing through p1 and p2. The computed points are used to
 //! initialize the  coordinates  of  cpOnThis  and  cpOnOther.  The  method
 //! returns FALSE if the lines are parallel or degenerate, TRUE otherwise.
 bool    closestPoints(const Point3c *this2, const Point3c *p1, const Point3c *p2, Point3c *cpOnThis, Point3c *cpOnOther) const;


 // FUNCTIONS BELOW THIS LINE MAY RETURN APPROXIMATE/NOT ROBUST RESULTS EVEN WHEN USING RATIONALS



 //! Distance from origin
 inline double length() const { return sqrt(TMESH_TO_DOUBLE(squaredLength())); }

 //! Divides the vector by its length. If isNull() the application exits with an error.
 Point3c& 	normalize();

 //! Rotates the vector around 'axis' by 'ang' radians ccw.
 void  	rotate(const Point3c& axis, const double& ang);

 //! Distance from 'b'
 inline double distance(const Point3c& b) const { return (((*(this)) - (b)).length()); }

 //! Distance from '*b'
 inline double distance(const Point3c *b) const { return (((*(this)) - (*b)).length()); }

 //! Distance from straight line through 'a' and 'b'
 double distanceFromLine(const Point3c *a, const Point3c *b) const;

 //! Distance from straight line through 'a' and 'b'. *cc is set to the closest line point.
 double distanceFromLine(const Point3c *a, const Point3c *b, Point3c *cc) const;

 double distanceFromEdge(const Point3c *a, const Point3c *b) const; //!< Distance from segment a-b

 //! Distance from segment a-b. *cc is set to the closest edge point.
 double distanceFromEdge(const Point3c *a, const Point3c *b, Point3c *cc) const;

 //! Distance between the straight lines through (this) - l1_p2 and l2_p1 - l2_p2.
 double distanceLineLine(const Point3c *l1_p2, const Point3c *l2_p1, const Point3c *l2_p2) const;

 //! Angle between this vector and 'v' in radians.
 double getAngle(const Point3c& v) const;

 //! Angle defined by <a, *this, b> in radians.
 inline double getAngle(const Point3c& a, const Point3c& b) const { return (a - (*this)).getAngle(b - (*this)); }

 //! Angle defined by <*a, *this, *b> in radians.
 inline double getAngle(const Point3c *a, const Point3c *b) const { return ((*a) - (*this)).getAngle((*b) - (*this)); }

 //! These functions round the coordinates to the closest floating point representation
 void snapToSinglePrecisionFloat() { x = float(TMESH_TO_NEAREST_DOUBLE(x)); y = float(TMESH_TO_NEAREST_DOUBLE(y)); z = float(TMESH_TO_NEAREST_DOUBLE(z)); }
 void snapToDoublePrecisionFloat() { x = TMESH_TO_NEAREST_DOUBLE(x); y = TMESH_TO_NEAREST_DOUBLE(y); z = TMESH_TO_NEAREST_DOUBLE(z); }

 //! Prints the coordinates of the point to a file handler. stdout is the default.
 void 	printPoint(FILE *fp = stdout) const { fprintf(fp, "%f %f %f,\n", TMESH_TO_FLOAT(x), TMESH_TO_FLOAT(y), TMESH_TO_FLOAT(z)); }		// Debug
};

//! Lexycographic comparison to be used with jqsort() or abstractHeap.
int xyzCompare(const void *p1, const void *p2);

//! Static point with 'infinite' coordinates.
extern const Point3c INFINITE_POINT;

// The following is wrong! -INF returns TRUE. Should be fixed...
//! Checks whether a point is INFINITE_POINT.
#define IS_FINITE_POINT(p) ((p).x < TMESH_INFINITY && (p).y < TMESH_INFINITY && (p).z < TMESH_INFINITY)


//! This is mainly for backward compatibility with older ImatiSTL-based Apps.
//! This class adds a generic 'info' field to store additional information.

class Point : public Point3c
{
public:
	 void *info;					//!< Further information

	//! Creates a new point with coordinates (0,0,0).
	inline Point() : Point3c() { }

	//! Creates a new point with the same coordinates as 's'. The info field is not copied.
	inline Point(const Point3c *s) : Point3c(s), info(NULL) { }

	//! Creates a new point with the same coordinates as 's'. The info field is not copied.
	inline Point(const Point3c& s) : Point3c(s), info(NULL) { }

	//! Creates a new point with the same coordinates as 's'. The info field is not copied.
	inline Point(const Point *s) : Point3c(s), info(NULL) { }

	//! Creates a new point with the same coordinates as 's'. The info field is not copied.
	inline Point(const Point& s) : Point3c(s), info(NULL) { }

	//! Creates a new point with coordinates (a,b,c).
	inline Point(const coord& a, const coord& b, const coord& c) : Point3c(a, b, c), info(NULL) { }
};

} //namespace T_MESH

#endif // _POINT_H

