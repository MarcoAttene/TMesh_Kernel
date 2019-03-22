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

#ifdef USE_LAZY_KERNEL

void lazy_num_base::compute_exact()
{
	tmesh_fraction *nv;

	if (approximate_value.isExact())
	{
		nv = new tmesh_fraction(approximate_value.low);
	}
	else
	{
		switch (operation)
		{
		case operation_id::none:
		nv = new tmesh_fraction(approximate_value.low);
		break;
		case operation_id::sum:
		nv = new tmesh_fraction(((lazy_num_binary *)this)->operand_1->exact() + ((lazy_num_binary *)this)->operand_2->exact());
		break;
		case operation_id::difference:
		nv = new tmesh_fraction(((lazy_num_binary *)this)->operand_1->exact() - ((lazy_num_binary *)this)->operand_2->exact());
		break;
		case operation_id::product:
		nv = new tmesh_fraction(((lazy_num_binary *)this)->operand_1->exact() * ((lazy_num_binary *)this)->operand_2->exact());
		break;
		case operation_id::division:
		nv = new tmesh_fraction(((lazy_num_binary *)this)->operand_1->exact() / ((lazy_num_binary *)this)->operand_2->exact());
		}
	}
#ifndef TMESH_SAFE_LOW_MEM_FOOTPRINT
#pragma omp critical
#endif
	{
		if (exact_value == NULL)
		{
			exact_value = nv;
#ifdef TMESH_SAFE_LOW_MEM_FOOTPRINT
			if (operation == operation_id::sum || operation == operation_id::difference || operation == operation_id::product || operation == operation_id::division)
			{
				if (((lazy_num_binary *)this)->operand_1.use_count()) ((lazy_num_binary *)this)->operand_1.reset();
				if (((lazy_num_binary *)this)->operand_2.use_count()) ((lazy_num_binary *)this)->operand_2.reset();
			}
#endif
			operation = operation_id::none;
			if (!approximate_value.isExact()) approximate_value = interval_number(*exact_value);
		} else delete nv;
	}
}

void lazy_num_base::print() const
{
	switch (operation)
	{
	case operation_id::none:
	printf("[ %f , %f]", approximate_value.low, approximate_value.high);
	return;
	case operation_id::sum:
	printf("( "); ((lazy_num_binary *)this)->operand_1->print(); printf(" + "); ((lazy_num_binary *)this)->operand_2->print(); printf(" )");
	return;
	case operation_id::difference:
	printf("( "); ((lazy_num_binary *)this)->operand_1->print(); printf(" - "); ((lazy_num_binary *)this)->operand_2->print(); printf(" )");
	return;
	case operation_id::product:
	printf("( "); ((lazy_num_binary *)this)->operand_1->print(); printf(" * "); ((lazy_num_binary *)this)->operand_2->print(); printf(" )");
	return;
	case operation_id::division:
	printf("( "); ((lazy_num_binary *)this)->operand_1->print(); printf(" / "); ((lazy_num_binary *)this)->operand_2->print(); printf(" )");
	}
}

#endif // USE_LAZY_KERNEL

#endif // USE_HYBRID_KERNEL

} //namespace T_MESH
