/****************************************************************************
* coordinates.h                                                             *
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

#ifndef _COORDINATES_H
#define _COORDINATES_H

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

namespace T_MESH
{

#define TMESH_NAN				NAN
#define TMESH_INFINITY			INFINITY

// This is the maximum value that can be represented in all the possible configurations
#define TMESH_MAX_COORDINATE	DBL_MAX

#ifdef WIN32
	typedef unsigned int fpu_status_t;
	inline fpu_status_t getFPURoundingMode() { return _controlfp(0, 0); }
	inline void setFPU53BitsPrecision() { _control87(_PC_53, _MCW_PC); }
	inline void setFPURoundingMode(fpu_status_t r) { _controlfp(r, _MCW_RC); }
	inline void setFPUModeToRoundUP() { _controlfp(_RC_UP, _MCW_RC); }
	inline void setFPUModeToRoundNEAR() { _controlfp(_RC_NEAR, _MCW_RC); }
#else
}

#include <fenv.h>
#if __APPLE__
#define _FPU_SETCW(cw) // nothing https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/float.3.html
#else
#include <fpu_control.h>
#endif

namespace T_MESH
{
	typedef int fpu_status_t;
	inline fpu_status_t getFPURoundingMode() { return fegetround(); }
	inline void setFPU53BitsPrecision() { int cw = 4722; _FPU_SETCW(cw); }
	inline void setFPURoundingMode(fpu_status_t r) { fesetround(r); }
	inline void setFPUModeToRoundUP() { fesetround(FE_UPWARD); }
	inline void setFPUModeToRoundNEAR() { fesetround(FE_TONEAREST); }

#endif
}

#ifdef WIN32
#pragma warning( disable : 4756)
#endif

#ifdef USE_HYBRID_KERNEL

#ifdef WIN32
#pragma warning( disable : 4056)
#pragma warning( disable : 4146)
#pragma warning( disable : 4800)
#endif

#include <mpir.h>

#include <stdint.h>
namespace T_MESH
{
	// Use the following to filter out sign and mantissa
#define EXPONENT_MASK	0x7FF0000000000000L

	// Use the following to filter out sign and exponent
#define MANTISSA_MASK	0x000FFFFFFFFFFFFFL

	//This is 2^(-52) << 52
#define MACH_DIFFEPS	0x0340000000000000L

	// Interpret 64 bits as either a double or a uint64_t
	// Useful to quickly play with bits in IEEE 754 standard representation
	typedef union error_approx_type_t
	{
		double d;
		uint64_t u;

		inline error_approx_type_t() {}
		inline error_approx_type_t(double a) : d(a) {}
		inline uint64_t is_negative() const { return u >> 63; }
		inline bool is_finite() const { return ((u & EXPONENT_MASK) != EXPONENT_MASK); }
		inline bool is_nan() const { return (u & MANTISSA_MASK); }
	} casted_double;

	// tmesh_fraction is essentially the same as mpq_class, but it correctly implements
	// the semantics of NaNs and Infinities as specified in IEEE 754 standard.
	// A more efficient solution is possible by tweaking mpir directly,
	// but this solution would require a fork to mpir which becomes tricky to update.
	// Hopefully, future versions of mpir will correctly handle these cases.
	// If and when this will happen, tmesh_fraction should be easily replaceable with
	// mpq_class.

	class tmesh_fraction
	{
	private:
		mpq_t mp;

		// Quick zero-check for mpz numbers
		static inline int is_mpz_zero(const __mpz_struct& z) { return (!z._mp_size); }

		// TRUE if mpq is NAN or INFINITY
		static inline bool is_nan_or_infinity(const mpq_t& op) { return is_mpz_zero(op->_mp_den); }

		// TRUE if mpq is ZERO
		static inline bool is_zero(const mpq_t& op) { return is_mpz_zero(op->_mp_num); }

		// TRUE if either arg is NAN or INFINITY
		static inline bool nan_or_infinity_ops(const mpq_t& op1, const mpq_t& op2) { return (is_nan_or_infinity(op1) || is_nan_or_infinity(op2)); }

		// TRUE if op is NAN
		static inline bool is_nan(const mpq_t& op) { return (is_mpz_zero(op->_mp_num) && is_mpz_zero(op->_mp_den)); }

		// TRUE if op is INFINITY
		static inline bool is_infinity(const mpq_t& op) { return (!is_mpz_zero(op->_mp_num) && is_mpz_zero(op->_mp_den)); }

		// TRUE if op is POSITIVE INFINITY
		static inline bool is_plus_infinity(const mpq_t& op) { return ((op->_mp_num._mp_size > 0) && is_mpz_zero(op->_mp_den)); }

		// TRUE if op is NEGATIVE INFINITY
		static inline bool is_minus_infinity(const mpq_t& op) { return ((op->_mp_num._mp_size < 0) && is_mpz_zero(op->_mp_den)); }

		// TRUE if op is POSITIVE (including infinity)
		static inline bool is_positive(const mpq_t& op) { return (op->_mp_num._mp_size > 0); }

		// TRUE if op is NEGATIVE (including infinity)
		static inline bool is_negative(const mpq_t& op) { return (op->_mp_num._mp_size < 0); }


	public:
		inline tmesh_fraction() { mpq_init(mp); }
		inline tmesh_fraction(const tmesh_fraction& a) { mpq_init(mp); mpq_set(mp, a.mp); }
		inline tmesh_fraction(mpq_srcptr q) { mpq_init(mp); mpq_set(mp, q); }
		inline tmesh_fraction(const mpz_srcptr &d) { mpq_init(mp); mpz_set(mpq_numref(mp), d); mpz_set_ui(mpq_denref(mp), 1); }
		inline tmesh_fraction(const mpz_srcptr &num, const mpz_srcptr &den) { mpq_init(mp); mpz_set(mpq_numref(mp), num); mpz_set(mpq_denref(mp), den); }
		inline tmesh_fraction(double a)
		{
			mpq_init(mp);
			casted_double d(a);
			if (d.is_finite()) mpq_set_d(mp, a);
			else if (d.is_nan()) setAsNAN();
			else if (d.is_negative()) setAsMinusInfinity();
			else setAsPlusInfinity();
		}

		inline void setAsNAN() { mpz_set_ui(mpq_numref(mp), 0); mpz_set_ui(mpq_denref(mp), 0); }
		inline void setAsPlusInfinity() { mpz_set_ui(mpq_numref(mp), 1); mpz_set_ui(mpq_denref(mp), 0); }
		inline void setAsMinusInfinity() { mpz_set_si(mpq_numref(mp), -1); mpz_set_ui(mpq_denref(mp), 0); }

		const mpq_t& get_mpq_t() const { return mp; }

		inline double get_d() const { return mpq_get_d(mp); }

		inline const mpz_srcptr get_den() const { return &mp->_mp_den; }
		inline const mpz_srcptr get_num() const { return &mp->_mp_num; }

		inline tmesh_fraction operator+(const tmesh_fraction& b) const
		{
			tmesh_fraction r;
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp) || (is_plus_infinity(mp) && is_minus_infinity(b.mp)) || (is_plus_infinity(b.mp) && is_minus_infinity(mp))) { r.setAsNAN(); return r; }
				if (is_plus_infinity(mp) || is_plus_infinity(b.mp)) { r.setAsPlusInfinity(); return r; }
				if (is_minus_infinity(mp) || is_minus_infinity(b.mp)) { r.setAsMinusInfinity(); return r; }
			}
			mpq_add(r.mp, mp, b.mp);
			return r;
		}

		inline tmesh_fraction operator-(const tmesh_fraction& b) const
		{
			tmesh_fraction r;
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp) || (is_plus_infinity(mp) && is_plus_infinity(b.mp)) || (is_minus_infinity(b.mp) && is_minus_infinity(mp))) { r.setAsNAN(); return r; }
				if (is_plus_infinity(mp) || is_minus_infinity(b.mp)) { r.setAsPlusInfinity(); return r; }
				if (is_minus_infinity(mp) || is_plus_infinity(b.mp)) { r.setAsMinusInfinity(); return r; }
			}
			mpq_sub(r.mp, mp, b.mp);
			return r;
		}

		inline tmesh_fraction operator*(const tmesh_fraction& b) const
		{
			tmesh_fraction r;
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp) || (is_infinity(mp) && is_zero(b.mp)) || (is_infinity(b.mp) && is_zero(mp))) { r.setAsNAN(); return r; }
				if ((is_negative(mp) && is_negative(b.mp)) || (is_positive(mp) && is_positive(b.mp))) { r.setAsPlusInfinity(); return r; }
				if ((is_negative(mp) && is_positive(b.mp)) || (is_positive(mp) && is_negative(b.mp))) { r.setAsMinusInfinity(); return r; }
			}
			mpq_mul(r.mp, mp, b.mp);
			return r;
		}

		inline tmesh_fraction operator/(const tmesh_fraction& b) const
		{
			tmesh_fraction r;
			if (nan_or_infinity_ops(mp, b.mp) || is_zero(b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp) || (is_zero(mp) && is_zero(b.mp)) || (is_infinity(b.mp) && is_infinity(mp))) { r.setAsNAN(); return r; }
				if (is_infinity(b.mp)) { return r; }
				if (is_positive(mp)) { r.setAsPlusInfinity(); return r; }
				if (is_negative(mp)) { r.setAsMinusInfinity(); return r; }
			}
			mpq_div(r.mp, mp, b.mp);
			return r;
		}

		inline tmesh_fraction& operator+=(const tmesh_fraction& b)
		{
			tmesh_fraction r = (*this) + (b);
			mpq_set(mp, r.mp);
			return *this;
		}

		inline tmesh_fraction& operator-=(const tmesh_fraction& b)
		{
			tmesh_fraction r = (*this) - (b);
			mpq_set(mp, r.mp);
			return *this;
		}

		inline tmesh_fraction& operator*=(const tmesh_fraction& b)
		{
			tmesh_fraction r = (*this) * (b);
			mpq_set(mp, r.mp);
			return *this;
		}

		inline tmesh_fraction& operator/=(const tmesh_fraction& b)
		{
			tmesh_fraction r = (*this) / (b);
			mpq_set(mp, r.mp);
			return *this;
		}

		inline bool operator==(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				return ((is_plus_infinity(mp) && is_plus_infinity(b.mp)) || (is_minus_infinity(mp) && is_minus_infinity(b.mp)));
			}
			return mpq_equal(mp, b.mp);
		}
		inline bool operator!=(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				return (!((is_plus_infinity(mp) && is_plus_infinity(b.mp)) || (is_minus_infinity(mp) && is_minus_infinity(b.mp))));
			}
			return !mpq_equal(mp, b.mp);
		}
		inline bool operator<(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp)) return false;
				return ((is_plus_infinity(b.mp) && !is_plus_infinity(mp)) || ((is_minus_infinity(mp) && !is_minus_infinity(b.mp))));
			}
			return (mpq_cmp(mp, b.mp)<0);
		}
		inline bool operator>(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp)) return false;
				return ((is_plus_infinity(mp) && !is_plus_infinity(b.mp)) || ((is_minus_infinity(b.mp) && !is_minus_infinity(mp))));
			}
			return (mpq_cmp(mp, b.mp)>0);
		}
		inline bool operator<=(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp)) return false;
				return (is_plus_infinity(b.mp) || is_minus_infinity(mp));
			}
			return (mpq_cmp(mp, b.mp) <= 0);
		}
		inline bool operator>=(const tmesh_fraction& b) const {
			if (nan_or_infinity_ops(mp, b.mp))
			{
				if (is_nan(mp) || is_nan(b.mp)) return false;
				return (is_plus_infinity(mp) || is_minus_infinity(b.mp));
			}
			return (mpq_cmp(mp, b.mp) >= 0);
		}

		inline bool operator==(int b) const { return (!is_nan_or_infinity(mp) && (mpq_cmp_si(mp, b, 1) == 0)); }

		inline bool operator!=(int b) const { return (is_nan_or_infinity(mp) || (mpq_cmp_si(mp, b, 1) != 0)); }

		inline bool operator<(int b) const { 
			if (is_nan_or_infinity(mp)) return is_minus_infinity(mp);
			return (mpq_cmp_si(mp, b, 1) < 0); 
		}
		inline bool operator>(int b) const { 
			if (is_nan_or_infinity(mp)) return is_plus_infinity(mp);
			return (mpq_cmp_si(mp, b, 1) > 0);
		}
		inline bool operator<=(int b) const { 
			if (is_nan_or_infinity(mp)) return is_minus_infinity(mp);
			return (mpq_cmp_si(mp, b, 1) <= 0);
		}
		inline bool operator>=(int b) const {
			if (is_nan_or_infinity(mp)) return is_plus_infinity(mp);
			return (mpq_cmp_si(mp, b, 1) >= 0);
		}

		inline int sign() const { return mpq_sgn(mp); }
	};

	inline int sgn(const tmesh_fraction& f) { return f.sign(); }
}

#ifdef USE_LAZY_KERNEL
#include "unprecise_numbers.h"

typedef T_MESH::lazy_num EXACT_NT;
#define EXACT_NT_FRACTION(x)	((x)->exact())
#define EXACT_NT_SIGN(x) ((x).sign())
#define EXACT_NT_TO_DOUBLE(x) ((x).approx())
#define EXACT_NT_DENOMINATOR(x) ((x)->exact().get_den())
#define EXACT_NT_NUMERATOR(x) ((x)->exact().get_num())
#define EXACT_NT_NAN	EXACT_NT(TMESH_NAN)
#define EXACT_NT_FPUT(x,f)	((x)->fput((f)))
#define EXACT_NT_FGET(x,f)	((x)->fget((f)))

#else

namespace T_MESH
{
typedef tmesh_fraction EXACT_NT;
}

#define EXACT_NT_FRACTION(x)	(*(x))
#define EXACT_NT_SIGN(x) (sgn(x))
#define EXACT_NT_TO_DOUBLE(x) ((x).get_d())
#define EXACT_NT_DENOMINATOR(x) ((x)->get_den())
#define EXACT_NT_NUMERATOR(x) ((x)->get_num())
#define EXACT_NT_NAN	EXACT_NT(0, 0)
#define EXACT_NT_FPUT(x,f)	gmp_fprintf((f), "%Qd", (x)->get_mpq_t())
#define EXACT_NT_FGET(x,f)	gmp_fscanf((f), "%Qd", (x)->get_mpq_t())

#endif // USE_LAZY_KERNEL

#endif // USE_HYBRID_KERNEL


namespace T_MESH
{

#ifdef USE_HYBRID_KERNEL

double to_upper_double(const tmesh_fraction& a);
double to_lower_double(const tmesh_fraction& a);

class PM_Rational
	{
		// This represents the current kernel precision (TRUE=exact, FALSE=approximated)
	public:
#ifdef WIN32
		static __declspec(thread) char use_rationals;
#else
		static thread_local char use_rationals;
#endif

		inline static bool isUsingRationals() { return use_rationals; }
		inline static void useRationals(bool v) { use_rationals = v; }

	protected:
		double _val_d;		// Double version. Always present.
		EXACT_NT *_val_p;	// Exact version. If NULL, this number is meant to represent a double

		inline EXACT_NT& getVal() { return *_val_p; }	// Exact value
		inline const EXACT_NT& getVal() const { return *_val_p; }

		inline double& getDVal() { return _val_d; }		// Double value
		inline const double& getDVal() const { return _val_d; }
		
		// Switch to rational
		inline void S2RF() { _val_p = new EXACT_NT(_val_d); } // If known to be not already rational
		inline void S2R() { if (!_val_p) S2RF(); } // With the check
		

	public:
		inline PM_Rational() : _val_p(NULL) {} // Undetermined double
		inline PM_Rational(const EXACT_NT& a) : _val_d(EXACT_NT_TO_DOUBLE(a)), _val_p(new EXACT_NT(a)) { }
		inline PM_Rational(float a) : _val_d(a), _val_p(NULL) { }
		inline PM_Rational(double a) : _val_d(a), _val_p(NULL) { }
		inline PM_Rational(int a) : _val_d(a), _val_p(NULL) { }
		inline PM_Rational(const PM_Rational& a) : _val_d(a._val_d) { if (a._val_p) _val_p = new EXACT_NT(*a._val_p); else _val_p = NULL; }

		inline ~PM_Rational() { if (_val_p) delete (_val_p); }

		inline EXACT_NT toRational() const { return (_val_p) ? (*_val_p) : (EXACT_NT(_val_d)); }

		inline bool isOfRationalType() const { return _val_p; }
		inline bool isOfDoubleType() const { return !_val_p; }

		inline int sign() const { return (_val_p) ? (EXACT_NT_SIGN((*_val_p))) : ((_val_d > 0) - (_val_d < 0)); }

		inline double toDouble() const { return _val_d; }
		inline int toInt() const { return int(toDouble()); }
		inline float toFloat() const { return float(toDouble()); }

		inline double toNearestDouble() const { return (_val_p) ? ((EXACT_NT_FRACTION(_val_p)).get_d()) : (_val_d); }
		inline double toUpperDouble() const { return (_val_p) ? (to_upper_double(EXACT_NT_FRACTION(_val_p))) : (_val_d); }
		inline double toLowerDouble() const { return (_val_p) ? (to_lower_double(EXACT_NT_FRACTION(_val_p))) : (_val_d); }

		inline void operator+=(const PM_Rational& a)
		{
			if (use_rationals) { S2R(); (*(_val_p)) += a.toRational(); } else { _val_d += a.toDouble(); }
		}

		inline void operator-=(const PM_Rational& a)
		{
			if (use_rationals) { S2R(); (*(_val_p)) -= a.toRational(); } else { _val_d -= a.toDouble(); }
		}

		inline void operator*=(const PM_Rational& a)
		{
			if (use_rationals) { S2R(); (*(_val_p)) *= a.toRational(); } else { _val_d *= a.toDouble(); }
		}

		inline void operator/=(const PM_Rational& a)
		{
			if (use_rationals) { S2R(); (*(_val_p)) /= a.toRational(); } else { _val_d /= a.toDouble(); }
		}

		inline PM_Rational operator+(const PM_Rational& a) const
		{
			if (use_rationals) return PM_Rational(toRational() + a.toRational());
			else return PM_Rational(toDouble() + a.toDouble());
		}

		inline PM_Rational operator-(const PM_Rational& a) const
		{
			if (use_rationals) return PM_Rational(toRational() - a.toRational());
			else return PM_Rational(toDouble() - a.toDouble());
		}

		inline PM_Rational operator*(const PM_Rational& a) const
		{
			if (use_rationals) return PM_Rational(toRational() * a.toRational());
			else return PM_Rational(toDouble() * a.toDouble());
		}

		inline PM_Rational operator/(const PM_Rational& a) const
		{
			if (use_rationals) return PM_Rational(toRational() / a.toRational());
			else return PM_Rational(toDouble() / a.toDouble());
		}

		inline bool operator==(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() == a.toRational());
			else return (_val_d == a._val_d);
		}

		inline bool operator!=(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() != a.toRational());
			else return (_val_d != a._val_d);
		}

		inline PM_Rational& operator=(const PM_Rational& a)
		{
			if (_val_p || a._val_p)
			{
				if (this == (&a)) return *this;
				if (_val_p) delete _val_p; 
				if (a._val_p) _val_p = new EXACT_NT(*a._val_p);
				else _val_p = NULL;
			}
			_val_d = a._val_d;
			return *this;
		}

		inline void setFromRational(const EXACT_NT& a)
		{
			if (_val_p) delete _val_p;
			_val_p = new EXACT_NT(a);
			_val_d = EXACT_NT_TO_DOUBLE(a);
		}

		inline void setFromDouble(double a)
		{
			if (_val_p) { delete _val_p; _val_p = NULL; }
			_val_d = a;
		}

		inline bool operator<(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() < a.toRational());
			else return (_val_d < a._val_d);
		}

		inline bool operator>(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() > a.toRational());
			else return (_val_d > a._val_d);
		}

		inline bool operator<=(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() <= a.toRational());
			else return (_val_d <= a._val_d);
		}

		inline bool operator>=(const PM_Rational& a) const
		{
			if (_val_p || a._val_p) return (toRational() >= a.toRational());
			else return (_val_d >= a._val_d);
		}

		friend char orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry);
		friend char orient3D(const class Point3c *t, const class Point3c *a, const class Point3c *b, const class Point3c *c);
		friend char inSphere3D(const class Point3c *pa, const class Point3c *pb, const class Point3c *pc, const class Point3c *pd, const class Point3c *pe);
		friend char inCircle3D(const class Point3c *pa, const class Point3c *pb, const class Point3c *pc, const class Point3c *pd);
		friend PM_Rational operator-(const PM_Rational& a);
		friend PM_Rational ceil(const PM_Rational& a);
		friend PM_Rational floor(const PM_Rational& a);
		friend PM_Rational round(const PM_Rational& a);

		friend class unprecise_number;

		inline void fput(FILE *fp) const { const EXACT_NT& n = toRational(); EXACT_NT_FPUT(&n, fp); }

		int fget(FILE *fp);
};

inline PM_Rational operator-(const PM_Rational& a)
{
	if (a.isOfRationalType()) return PM_Rational(-(*(a._val_p)));
	else return PM_Rational(-a._val_d);
}

PM_Rational ceil(const PM_Rational& a);
PM_Rational floor(const PM_Rational& a);
PM_Rational round(const PM_Rational& a);
inline PM_Rational fabs(const PM_Rational& a) {	if (a < 0) return -a; else return a; }

/**************** I/O operators ****************/

#define TMESH_TO_DOUBLE(x) ((x).toDouble())
#define TMESH_TO_FLOAT(x) ((x).toFloat())
#define TMESH_TO_INT(x) ((x).toInt())
#define TMESH_TO_NEAREST_DOUBLE(x) ((x).toNearestDouble())
#define TMESH_TO_LOWER_DOUBLE(x) ((x).toLowerDouble())
#define TMESH_TO_UPPER_DOUBLE(x) ((x).toUpperDouble())
#define TMESH_IS_ZERO(x) (!(x).sign())

#define TMESH_OMP_CLAUSES firstprivate(tmesh_thread_initializer) copyin(PM_Rational::use_rationals)

#else

typedef double PM_Rational;

#define TMESH_TO_DOUBLE(x) (x)
#define TMESH_TO_FLOAT(x) ((float)(x))
#define TMESH_TO_INT(x) ((int)(x))
#define TMESH_TO_NEAREST_DOUBLE(x) (x)
#define TMESH_TO_LOWER_DOUBLE(x) (x)
#define TMESH_TO_UPPER_DOUBLE(x) (x)
#define TMESH_IS_ZERO(x) ((x)==0)

#define TMESH_OMP_CLAUSES

#endif

typedef PM_Rational coord;

#define TMESH_DETERMINANT3X3(a11, a12, a13, a21, a22, a23, a31, a32, a33) ((a11)*((a22)*(a33) - (a23)*(a32)) - (a12)*((a21)*(a33) - (a23)*(a31)) + (a13)*((a21)*(a32) - (a22)*(a31)))

} //namespace T_MESH

#endif //_COORDINATES_H

