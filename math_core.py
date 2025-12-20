from coeffs import *
from math import frexp, ldexp
from fractions import Fraction

# Helper functions
def poly_eval(x, coeffs):
    """Evaluate sum_{k=0}^n coeffs[k] * x^k using Horner's method."""
    y = 0
    i = len(coeffs) - 1
    while i >= 0:
        y = y * x + coeffs[i]
        i -= 1
    return y

def rational_eval(x, num_coeffs, den_coeffs):
    """Evaluate P(x)/Q(x) where P, Q are given by coefficient lists."""
    p = poly_eval(x, num_coeffs)
    q = poly_eval(x, den_coeffs)
    return p / q

def _as_fraction(x):
    """Coerce inputs to Fraction for exact remainder arithmetic."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, float):
        return Fraction.from_float(x)
    return Fraction(x)

def _reduce_angle(x):
    """Reduce angle to [-pi, pi] exactly using Fraction."""
    xf = _as_fraction(x)
    # Modulo with Fraction preserves precision.
    # Result of % TAU is in [0, TAU)
    r = (xf + PI) % TAU
    return r - PI

def _reduce_pi(x):
    """Reduce angle to [-pi/2, pi/2) by mod pi (exact Fraction)."""
    # x can be float/Fraction
    r = _reduce_angle(x)  # [-pi, pi] Fraction
    if r >= PI_2:
        r -= PI
    elif r < -PI_2:
        r += PI
    return r

def _pade_about_1(x, num_coeffs, den_coeffs):
    """Evaluate Padé approximant about x=1 using t = x-1."""
    t = _as_fraction(x) - 1
    return rational_eval(t, num_coeffs, den_coeffs)

# Implementations
def sqrt(x):
    """Square root via frexp/ldexp range reduction + Padé seed + 2 Newton iterations."""
    xf = float(x)
    if xf < 0.0: # Hanlding for complex set is not implemented yet.
        raise ValueError("sqrt domain")
    if xf == 0.0:
        return 0.0

    m, e = frexp(xf) # making floated number x in format m * 2^e. m is in [0.5, 1) and is a float and e is an integer.
    if e  &1: # if e is odd, then m is in [0.25, 0.5)
        m *= 2.0
        e -= 1
    # if e is even, then m is in [0.5, 1)

    y0 = float(_pade_about_1(m, numerator_sqrt_13_13, denominator_sqrt_13_13))
    y = 0.5 * (y0 + m / y0)
    y = 0.5 * (y + m / y)

    return ldexp(y, e // 2)

def cbrt(x):
    """Cube root via frexp/ldexp range reduction + Padé seed + 2 Newton iterations."""
    xf = float(x)
    if xf == 0.0:
        return 0.0

    sign = 1.0
    if xf < 0.0:
        sign = -1.0
        xf = -xf

    m, e = frexp(xf)
    q = e // 3
    r = e - 3 * q
    a = ldexp(m, r)

    y0 = float(_pade_about_1(a, numerator_cbrt_13_13, denominator_cbrt_13_13))
    y = (2.0 * y0 + a / (y0 * y0)) / 3.0
    y = (2.0 * y + a / (y * y)) / 3.0

    return sign * ldexp(y, q)

def EXP(x):
    xf = float(x)
    k = int(round(xf / LN2_F))
    r = _as_fraction(xf) - k * LN2
    e_power_r = rational_eval(r, numerator_exp, denumerator_exp)
    return float(e_power_r) * ldexp(1.0, k)

def _sin_small(x):
    x2 = x * x
    f = poly_eval(x2, coeffs_sin)
    return x * f

def _cos_small(x):
    x2 = x * x
    return poly_eval(x2, coeffs_cos)

def sin(x):
    xr = _reduce_angle(x)
    if xr < 0:
        return float(-sin(-xr))
    if xr > PI_2:
        xr = PI - xr
    if xr <= PI_4:
        return float(_sin_small(xr))
    return float(_cos_small(PI_2 - xr))

def cos(x):
    xr = _reduce_angle(x)
    if xr < 0:
        xr = -xr
    if xr > PI_2:
        return float(-cos(PI - xr))
    if xr <= PI_4:
        return float(_cos_small(xr))
    return float(_sin_small(PI_2 - xr))

_POLE_EPS = 1e-16 # used for pole detection for undefined values of trig functions.

def tan(x):
    r2 = _reduce_angle(x)
    if abs(r2 - PI_2) < _POLE_EPS:
        return float("inf")
    if abs(r2 + PI_2) < _POLE_EPS:
        return float("-inf")
    r = _reduce_pi(x)
    if abs(abs(r) - PI_2) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")
    s = sin(r)
    c = cos(r)
    if abs(c) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")
    return s / c

def cot(x):
    r = _reduce_pi(x)
    if abs(r) < _POLE_EPS:
        return float("inf")
    s = sin(r)
    c = cos(r)
    if abs(s) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")
    return c / s

def sinh(x):
    return (0.5*(EXP(x) - EXP(-x)))

def cosh(x):
    return (0.5*(EXP(x) + EXP(-x)))

def _atan_series(z):
    z2 = z * z
    return z * poly_eval(z2, coeffs_atan)

def atan(x):
    xf = float(x)
    if xf == 0.0:
        return 0.0
    sign = 1.0
    if xf < 0.0:
        sign = -1.0
        xf = -xf
    if xf > 1.0:
        return sign * (PI_2_F - atan(1.0 / xf))
    n = 0
    t = xf
    while t > 0.3:
        t = t / (1.0 + sqrt(1.0 + t * t))
        n += 1
    a = _atan_series(t)
    return sign * (a * (2.0 ** n))

def atan2(y, x):
    yf = float(y)
    xf = float(x)
    if xf > 0.0:
        return atan(yf / xf)
    if xf < 0.0:
        if yf >= 0.0:
            return atan(yf / xf) + PI_F
        return atan(yf / xf) - PI_F
    if yf > 0.0:
        return PI_2_F
    if yf < 0.0:
        return -PI_2_F
    return 0.0

def asin(x):
    xf = float(x)
    if xf > 1.0 or xf < -1.0:
        raise ValueError("asin domain")
    if xf == 1.0:
        return PI_2_F
    if xf == -1.0:
        return -PI_2_F
    denom = sqrt(1.0 - xf * xf)
    return atan(xf / denom)

def acos(x):
    return PI_2_F - asin(x)

def _atanh_series(z):
    z2 = z * z
    return z * poly_eval(z2, coeffs_atanh)

def ln1p(u):
    uf = float(u)
    if uf <= -1.0:
        raise ValueError("ln1p domain")
    if uf == 0.0:
        return 0.0
    m = 1.0 + uf
    z = (m - 1.0) / (m + 1.0)
    return 2.0 * _atanh_series(z)

def ln(x):
    xf = float(x)
    if xf <= 0.0:
        raise ValueError("Can't evaluate in complex set")
    m, e = frexp(xf)
    z = (m - 1.0) / (m + 1.0)
    return 2.0 * _atanh_series(z) + e * LN2_F

def sinc(x):
    if x==0:
        return 1
    else:
        return rational_eval(x, numerator_sinc, denumerator_sinc)
