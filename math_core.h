#ifndef MATH_CORE_H
#define MATH_CORE_H
#include <stdint.h>
typedef double Real;
// Custom IEEE 754 Primitives
Real mc_frexp(Real value, int *exp);
Real mc_ldexp(Real value, int exp);
Real mc_fabs(Real x);
int mc_isnan(Real x);
int mc_isinf(Real x);
Real mc_modf(Real x, Real *iptr);
Real mc_fmod(Real x, Real y);
Real mc_round(Real x);
Real mc_floor(Real x); 
Real mc_ceil(Real x);

// Constants exposed for testing
extern const Real PI_C;
extern const Real PI_4_C;
extern const Real E_F_C;
extern const Real PI_2_C;

// Helpers
Real poly_eval(Real x, const Real* coeffs, int n);
Real rational_eval(Real x, const Real* num, int num_len, const Real* den, int den_len);

// Core reductions
int rem_pio2(Real x, Real *y);

// Functions
Real mc_sqrt(Real x);
Real mc_cbrt(Real x);
Real mc_exp(Real x);
Real mc_ln(Real x);
Real mc_ln1p(Real u);

Real mc_sin(Real x);
Real mc_cos(Real x);
Real mc_tan(Real x);
Real mc_cot(Real x);

Real mc_sinh(Real x);
Real mc_cosh(Real x);

Real mc_atan(Real x);
Real mc_atan2(Real y, Real x);
Real mc_asin(Real x);
Real mc_acos(Real x);

Real mc_sinc(Real x);

#endif
