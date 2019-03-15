/****************************************************************************
* coordinates.cpp                                                           *
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

#include "coordinates.h"
#include <limits>
#include <cmath>

namespace T_MESH
{
#ifdef USE_HYBRID_KERNEL

// Default behaviour = FILTERED KERNEL
#ifdef WIN32
	char PM_Rational::use_rationals = 0;
#else
	thread_local char PM_Rational::use_rationals = 0;
#endif

	double to_upper_double(const tmesh_fraction& a)
	{
		if (sgn(a) <= 0) return a.get_d(); // Default truncation is equivalent to round-up

		casted_double eps;
		double ret = eps.d = a.get_d();
		eps.u &= EXPONENT_MASK;
		if (eps.u == EXPONENT_MASK) return ret; // Not a finite number
		if (eps.u>MACH_DIFFEPS)
		{
			eps.u -= MACH_DIFFEPS;
			return (ret + eps.d);
		} else return DBL_MIN;
	}

	double to_lower_double(const tmesh_fraction& a)
	{
		if (sgn(a) >= 0) return a.get_d(); // Default truncation is equivalent to round-down

		casted_double eps;
		double ret = eps.d = a.get_d();
		eps.u &= EXPONENT_MASK;
		if (eps.u == EXPONENT_MASK) return ret; // Not a finite number
		if (eps.u>MACH_DIFFEPS)
		{
			eps.u -= MACH_DIFFEPS;
			return (ret - eps.d);
		} else return -DBL_MIN;
	}

int PM_Rational::fget(FILE *fp)
{
	tmesh_fraction a;
	if (gmp_fscanf(fp, "%Qd", a.get_mpq_t()))
	{
		double d = a.get_d();
		if (a == d) setFromDouble(d);
		else setFromRational(EXACT_NT(a));
		return 1;
	}
	return 0;
}

PM_Rational ceil(const PM_Rational& a)
{
	if (a.isOfRationalType())
	{
		mpz_t n, d, f;
		const EXACT_NT& en = a.toRational();
		mpz_init(n); mpz_init(d); mpz_init(f);
		mpz_set(n, EXACT_NT_NUMERATOR(&en));
		mpz_set(d, EXACT_NT_DENOMINATOR(&en));
		mpz_cdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(tmesh_fraction(f));
	}
	else
		return PM_Rational(::ceil(a.getDVal()));
}

PM_Rational floor(const PM_Rational& a)
{
	if (a.isOfRationalType())
	{
		mpz_t n, d, f;
		const EXACT_NT& en = a.toRational();
		mpz_init(n); mpz_init(d); mpz_init(f);
		mpz_set(n, EXACT_NT_NUMERATOR(&en));
		mpz_set(d, EXACT_NT_DENOMINATOR(&en));
		mpz_fdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(tmesh_fraction(f));
	} else
		return PM_Rational(::floor(a.getDVal()));
}

PM_Rational round(const PM_Rational& a)
{
	if (a.isOfRationalType())
	{
		mpz_t n, d, f, c;
		mpz_init(n); mpz_init(d); mpz_init(f); mpz_init(c);
		const EXACT_NT& en = a.toRational();
		mpz_set(n, EXACT_NT_NUMERATOR(&en));
		mpz_set(d, EXACT_NT_DENOMINATOR(&en));
		mpz_cdiv_q(c, n, d);
		mpz_clear(n); mpz_clear(d);
		PM_Rational fr = PM_Rational(tmesh_fraction(f));
		PM_Rational cr = PM_Rational(tmesh_fraction(c));
		mpz_clear(f); mpz_clear(c);
		return ((a - fr) < (cr - a)) ? (fr) : (cr);
	} else
		return PM_Rational(::round(a.getDVal()));
}

#endif

} //namespace T_MESH
