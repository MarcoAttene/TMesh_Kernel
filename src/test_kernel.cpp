/****************************************************************************
* test_kernel.cpp                                                           *
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

#include "tmesh_kernel.h"
#include <iostream>
#include <time.h>

using namespace T_MESH;

int main()
{
	TMesh::init(); // This is mandatory to correctly set up the FPU for exact computation

	// The type 'coord' represents the basic hybrid number in TMesh.
	// It can be handled as an IEEE 754 standard double in most cases.
	coord q_x = 3.0;
	coord q_y = 4.0;
	coord q_z = 5.0;

	// for example...
	q_x += (q_y - q_z);

	// When necessary, the TMESH_TO_DOUBLE macro approximates a coord to an actual double
	double sqrtex = sqrt( TMESH_TO_DOUBLE(q_y) );

	// And a double can be seamlessly assigned to a coord
	q_z = sqrtex;
	
	// A Point in TMesh is a triplet of coordinates, each having type 'coord'
	// The following example shows how to manipulate points as vectors
	// Here we calculate the projection of 'p' on the plane by A, B and C
	Point q(q_x, q_y, q_z);
	Point A(0, 0, 0);
	Point B(2, 0, 0);
	Point C(0, 2, 0);

	Point BA = B - A;
	Point CA = C - A;
	Point plane_vec = BA & CA; // Operator '&' represents the cross product

	// The previous three lines can be compacted into a single line as follows:
	// Point plane_vec = (B - A) & (C - A);

	// Here we calculate the projection as the intersection of a straight line by
	// q and orthogonal to the plane.
	Point lifted_q = q + plane_vec;
	Point projected_point = Point::linePlaneIntersection(q, lifted_q, A, B, C);

	std::cout << "Projected point coordinates: ";
	projected_point.printPoint();


	// By default, TMesh works in floating point mode

	// Declare four coplanar points
	Point p1(5, 1, 0);
	Point p2(4, 2, 0);
	Point p3(8, 1, 1);
	Point p4(7, 2, 1);

	// Predicates in TMesh are always exact - the following reports coplanarity as expected
	// even if we work with floating point coordinates
	if (p1.exactOrientation(p2, p3, p4) == 0) std::cout << "Points pi are coplanar\n";
	else std::cout << "Points pi are not coplanar\n";

	// We now create scaled copies of the pi's
	Point q1 = p1 / 3;
	Point q2 = p2 / 3;
	Point q3 = p3 / 3;
	Point q4 = p4 / 3;

	// Coordinates in the qi's are not exact scales of those in the pi's because we are
	// approximating them using floating point numbers
	if (q1.exactOrientation(q2, q3, q4) == 0) std::cout << "Points qi are coplanar\n";
	else std::cout << "Points qi are not coplanar\n";

	// Now we switch to rational mode
	TMesh::useRationals();

	// We create other scaled copies of the pi's
	Point r1 = p1 / 3;
	Point r2 = p2 / 3;
	Point r3 = p3 / 3;
	Point r4 = p4 / 3;

	// Coordinates in the ri's are exact scales of those in the pi's because we are
	// representing them using rational numbers
	if (r1.exactOrientation(r2, r3, r4) == 0) std::cout << "Points ri are coplanar\n";
	else std::cout << "Points ri are not coplanar\n";

	// Now we switch back to floating point mode
	TMesh::useRationals(false);

	// Even if we switched back to floating point, coordinates in the ri's are still
	// rational numbers. That is why the following still reports coplanarity.
	if (r1.exactOrientation(r2, r3, r4) == 0) std::cout << "Points ri are coplanar\n";
	else std::cout << "Points ri are not coplanar\n";

	// The following is just an example of hybrid calculation where q1 has floating
	// point coords while the ri's have rational coords.
	if (q1.exactOrientation(r2, r3, r4) == 0) std::cout << "Points (q1, r2, r3, r4) are coplanar\n";
	else std::cout << "Points (q1, r2, r3, r4) are not coplanar\n";


	// Multi-threading example
	
	const int num_tests = 1000000;

	// Switch to rational mode to make Points assume exact coordinates.
	// We also divide by three to have a high probability that the resulting
	// coordinates cannot be represented using floating points.
	// In this way, we are stressing TMesh on purpose by making it use
	// rationals in most places to calculate a large number of insphere predicates.
	TMesh::useRationals(true);

	clock_t c0 = clock();
	for (int i = 0; i < num_tests; i++)
	{
		Point p1(rand(), rand(), rand()); p1 /= 3;
		Point p2(rand(), rand(), rand()); p2 /= 3;
		Point p3(rand(), rand(), rand()); p3 /= 3;
		Point p4(rand(), rand(), rand()); p4 /= 3;
		Point p5(rand(), rand(), rand()); p5 /= 3;
		p1.inSphere(&p2, &p3, &p4, &p5);
	}

	std::cout << "Needed " << ((clock() - c0) / 1000) << " seconds to sequentially calculate " << num_tests << " insphere predicates.\n";

	// Here we use OpenMP to concurrently redo the same calculation.
	// Remember to always use 'firstprivate(tmesh_thread_initializer)' to correctly
	// set the FPU mode in all the threads.

	c0 = clock();
#pragma omp parallel for firstprivate(tmesh_thread_initializer) schedule(dynamic) 
	for (int i = 0; i < num_tests; i++)
	{
		Point p1(rand(), rand(), rand()); p1 /= 3;
		Point p2(rand(), rand(), rand()); p2 /= 3;
		Point p3(rand(), rand(), rand()); p3 /= 3;
		Point p4(rand(), rand(), rand()); p4 /= 3;
		Point p5(rand(), rand(), rand()); p5 /= 3;
		p1.inSphere(&p2, &p3, &p4, &p5);
	}

	std::cout << "Needed " << ((clock() - c0) / 1000) << " seconds to concurrently calculate " << num_tests << " insphere predicates.\n";

	return 0;
}
