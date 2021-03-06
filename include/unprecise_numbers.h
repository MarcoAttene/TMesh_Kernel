/****************************************************************************
* unprecise_numbers.h                                                       *
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

#ifndef UNPRECISE_NUMBERS
#define UNPRECISE_NUMBERS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <memory>
#include "coordinates.h"

namespace T_MESH
{
	// Interval_number
	class interval_number
	{
	public:
		double low, high;

		inline interval_number() {}
		inline interval_number(double a) : low(a), high(a) {}
		inline interval_number(double a, double b) : low(a), high(b) {}
		inline interval_number(const interval_number& b) : low(b.low), high(b.high) {}

#ifdef USE_HYBRID_KERNEL
		inline interval_number(const tmesh_fraction& a) { init(a); }

		inline void init(const tmesh_fraction& a)
		{
			casted_double eps;
			low = high = eps.d = a.get_d();
			eps.u &= EXPONENT_MASK;
			if (eps.u == EXPONENT_MASK) return; // Not a finite number (NaN or INFINITY)
			int s = sgn(a);
			if (eps.u>MACH_DIFFEPS)
			{
				eps.u -= MACH_DIFFEPS;
				if (s>0) high += eps.d;
				if (s<0) low -= eps.d;
			} else
			{
				if (s<0) low = -DBL_MIN;
				if (s>0) high = DBL_MIN;
			}
		}
#endif
		inline void init(const interval_number& b) { low = b.low; high = b.high; }
		inline void init(const double b) { low = high = b; }

		inline bool disjointWith(const interval_number& b) const { return (high < b.low || low > b.high); }

		inline bool isExact() const { return (low == high); }
		inline bool signIsReliable() const { return (low>0 || high <0 || (low==0 && high==0)); }
		inline char sign() const { return (low > 0) ? (1) : ((low < 0) ? (-1) : (0)); }

		inline double getMid() const { return (low + high)*0.5; }

		inline bool operator<(const interval_number& b) const { return (high<b.low); }
		inline bool operator>(const interval_number& b) const { return (low>b.high); }
		inline bool operator<=(const interval_number& b) const { return (high <= b.low); }
		inline bool operator>=(const interval_number& b) const { return (low >= b.high); }

		inline interval_number& operator=(const interval_number& b) { low = b.low; high = b.high; return *this; }

		inline interval_number& operator+=(const interval_number& b) { low = -((-low) - b.low); high += b.high; return *this; }
		inline interval_number& operator-=(const interval_number& b) { low = -(b.high - low); high -= b.low; return *this; }
		inline interval_number& operator*=(const interval_number& b) { return operator=((*this) * b); } // This should be optimized
		inline interval_number& operator/=(const interval_number& b) { return operator=((*this) / b); } // This should be optimized

		inline interval_number operator+(const interval_number& b) const { return interval_number(-((-low)-b.low), high+b.high); }

		inline interval_number operator-(const interval_number& b) const { return interval_number(-(b.high - low), high - b.low); }

		inline interval_number operator*(const interval_number& b) const
		{
			casted_double l1(low), h1(high), l2(b.low), h2(b.high);
			uint64_t conf = (l1.is_negative() << 3) + (h1.is_negative() << 2) + (l2.is_negative() << 1) + (h2.is_negative());
			switch (conf)
			{
			case 0: return interval_number(-((-low)*b.low), high*b.high);
			case 2: return interval_number(-((-high)*b.low), high*b.high);
			case 3: return interval_number(-((-high)*b.low), low*b.high);
			case 8: return interval_number(-((-low)*b.high), high*b.high);
			case 10:
				double ll, lh, hl, hh;
				ll = low*b.low; lh = -((-low)*b.high); hl = -((-high)*b.low); hh = high*b.high;
				if (hl < lh) lh = hl;
				if (ll > hh) hh = ll;
				return interval_number(lh, hh);
			case 11: return interval_number(-((-high)*b.low), low*b.low);
			case 12: return interval_number(-((-low)*b.high), high*b.low);
			case 14: return interval_number(-((-low)*b.high), low*b.low);
			case 15: return interval_number(-((-high)*b.high), low*b.low);
			};
			//TMesh::error("interval_number: inconsistent interval.");
			return interval_number(TMESH_NAN);
		}

		inline interval_number operator/(const interval_number& b) const
		{
			casted_double l1(low), h1(high), l2(b.low), h2(b.high);
			uint64_t conf = (l1.is_negative() << 3) + (h1.is_negative() << 2) + (l2.is_negative() << 1) + (h2.is_negative());
			switch (conf)
			{
			case 0: return interval_number(-((-low)/b.high), high/b.low);
			case 2: return interval_number(TMESH_NAN);
			case 3: return interval_number(-((-high) / b.high), low / b.low);
			case 8: return interval_number(-((-low) / b.low), high / b.low);
			case 10: return interval_number(TMESH_NAN);
			case 11: return interval_number(-((-high) / b.high), low / b.high);
			case 12: return interval_number(-((-low) / b.low), high / b.high);
			case 14: return interval_number(-TMESH_NAN);
			case 15: return interval_number(-((-high) / b.low), low / b.high);
			};
			//TMesh::error("interval_number: inconsistent interval.");
			return interval_number(TMESH_NAN);
		}
	};




	//////////////////////////////////////////////// LAZY NUMBERS //////////////////////////////////////////////

#ifdef USE_LAZY_KERNEL

typedef enum
{
	none, 
	sum, difference, product, division
} operation_id;


// Two possibilities to parallelize the exact computation of a lazy number:
// 1) Safe and lower memory footprint
//    > the subdags are released immediately upon exact calculation
// 2) Higher memory footprint but much faster concurrency (seems to be safe too)
//    > the subdags are lazily released upon destruction of the last number that used them

//#define TMESH_SAFE_LOW_MEM_FOOTPRINT


// This is the number's DAG


class lazy_num_base
{
protected:
	interval_number approximate_value;
	tmesh_fraction *exact_value; // NULL, unless previously computed
	operation_id operation; // If 'none', this is a leaf

public:
	inline lazy_num_base() : exact_value(NULL) {}
		
	inline lazy_num_base(const operation_id o) : exact_value(NULL), operation(o) {}

	inline lazy_num_base(double d) :
		approximate_value(d),
		exact_value(NULL),
		operation(operation_id::none)
	{
	}

	inline lazy_num_base(const tmesh_fraction& a) :
		approximate_value(a),
		exact_value(new tmesh_fraction(a)),
		operation(operation_id::none)
	{
	}

	inline ~lazy_num_base()
	{
		if (exact_value) delete exact_value;
	}

	inline const interval_number& unprecise() const { return approximate_value; }

	void print() const;

	inline const tmesh_fraction& exact()
	{
#ifdef TMESH_SAFE_LOW_MEM_FOOTPRINT
#pragma omp critical
#endif
		if (!exact_value) compute_exact();
		return *exact_value;
	}

private:
	void compute_exact();
};

class lazy_num_binary : public lazy_num_base
{
public:
	std::shared_ptr<lazy_num_base> operand_1, operand_2;

	inline lazy_num_binary(const std::shared_ptr<lazy_num_base>& a, const std::shared_ptr<lazy_num_base>& b, const operation_id op) : lazy_num_base(op),
		operand_1(a),
		operand_2(b)
	{
		switch (op)
		{
		case operation_id::sum:
		approximate_value = operand_1->unprecise() + operand_2->unprecise();
		return;
		case operation_id::difference:
		approximate_value = operand_1->unprecise() - operand_2->unprecise();
		return;
		case operation_id::product:
		approximate_value = operand_1->unprecise() * operand_2->unprecise();
		return;
		case operation_id::division:
		approximate_value = operand_1->unprecise() / operand_2->unprecise();
		return;
		default:
		return;
		}
	}
};


class lazy_num
{
	std::shared_ptr<lazy_num_base> root;

public:
	inline lazy_num(double d) { root = std::make_shared<lazy_num_base>(d); }
	inline lazy_num(const tmesh_fraction& a) { root = std::make_shared<lazy_num_base>(a); }
	inline lazy_num(const lazy_num& n) : root(n.root) { }
	inline lazy_num(const lazy_num& o1, const lazy_num& o2, const operation_id& op) { root = std::make_shared<lazy_num_binary>(o1.root, o2.root, op); }
	
	inline lazy_num& operator=(const lazy_num& n) {	root = n.root; return *this; }

	inline const interval_number& unprecise() const { return root->unprecise(); }
	inline const double approx() const { return root->unprecise().getMid(); }
	inline const tmesh_fraction& exact() const { return root->exact(); }
	inline bool isPrecise() const { return (root->unprecise().isExact()); }

	inline lazy_num operator+(const lazy_num& n) const { return lazy_num(*this, n, operation_id::sum); }
	inline lazy_num operator-(const lazy_num& n) const { return lazy_num(*this, n, operation_id::difference); }
	inline lazy_num operator*(const lazy_num& n) const { return lazy_num(*this, n, operation_id::product); }
	inline lazy_num operator/(const lazy_num& n) const { return lazy_num(*this, n, operation_id::division); }

	inline lazy_num& operator+=(const lazy_num& n) { root = std::make_shared<lazy_num_binary>(root, n.root, operation_id::sum); return *this; }
	inline lazy_num& operator-=(const lazy_num& n) { root = std::make_shared<lazy_num_binary>(root, n.root, operation_id::difference); return *this; }
	inline lazy_num& operator*=(const lazy_num& n) { root = std::make_shared<lazy_num_binary>(root, n.root, operation_id::product); return *this; }
	inline lazy_num& operator/=(const lazy_num& n) { root = std::make_shared<lazy_num_binary>(root, n.root, operation_id::division); return *this; }

	void print() const { root->print(); }

	void fput(FILE *fp) const { gmp_fprintf(fp, "%Qd", exact().get_mpq_t()); }

	int fget(FILE *fp)
	{
		root = std::make_shared<lazy_num_base>(0.0);
		return gmp_fscanf(fp, "%Qd", root->exact().get_mpq_t());
	}

	inline int sign() const
	{
		if (unprecise().signIsReliable()) return unprecise().sign();
		else return sgn(exact());
	}

	inline bool operator==(const lazy_num& n) const
	{
		if (unprecise().disjointWith(n.unprecise())) return false;
		if (isPrecise() && n.isPrecise()) return (unprecise().low == n.unprecise().low);
		return (exact() == n.exact());
	}

	inline bool operator!=(const lazy_num& n) const
	{
		if (unprecise().disjointWith(n.unprecise())) return true;
		if (isPrecise() && n.isPrecise()) return (unprecise().low != n.unprecise().low);
		return (exact() != n.exact());
	}

	inline bool operator<(const lazy_num& n) const
	{
		if (unprecise() < n.unprecise()) return true;
		if (unprecise() >= n.unprecise()) return false;
		return (exact()<n.exact());
	}

	inline bool operator<=(const lazy_num& n) const
	{
		if (unprecise() <= n.unprecise()) return true;
		if (unprecise() > n.unprecise()) return false;
		return (exact() <= n.exact());
	}

	inline bool operator>(const lazy_num& n) const
	{
		if (unprecise() > n.unprecise()) return true;
		if (unprecise() <= n.unprecise()) return false;
		return (exact()>n.exact());
	}

	inline bool operator>=(const lazy_num& n) const
	{
		if (unprecise() >= n.unprecise()) return true;
		if (unprecise() < n.unprecise()) return false;
		return (exact() >= n.exact());
	}
};



#endif // USE_LAZY_KERNEL

}

#endif // UNPRECISE_NUMBERS
