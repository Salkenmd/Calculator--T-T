#README.md. Avoid built in modules and functions as much as possible. Some direct hardware communication is allowed but this is still a project where I try to learn numerical analysis.
# Indev version 1.1
from fractions import Fraction
from math import frexp, ldexp

F = Fraction

# Float versions for fast range-reduction / UI-facing ops.
LN2_F = float(F("0.693147180559945309417232121458176568"))
PI_F = float(F("3.1415926535897932384626433832795028841971693993751058"))
TAU_F = 2.0 * PI_F
PI_2_F = 0.5 * PI_F
PI_4_F = 0.25 * PI_F
E_F = 2.71828182845904523536028747135266249775724709369995

### Pade approximants [13/13]
# Store rational coefficients exactly so you can isolate approximation error from FP error.
numerator_LN1P = [
    F(1),
    F(6),
    F(4741, 300),
    F(1441, 60),
    F(107091, 4600),
    F(8638, 575),
    F(263111, 40250),
    F(153081, 80500),
    F(395243, 1101240),
    F(28549, 688275),
    F(605453, 228813200),
    F(785633, 10296594000),
    F(1145993, 1873980108000),
    F(0),
]

denominator_LN1P = [
    F(1),
    F(13, 2),
    F(468, 25),
    F(1573, 50),
    F(1573, 46),
    F(11583, 4600),
    F(10296, 805),
    F(2574, 575),
    F(11583, 10925),
    F(143, 874),
    F(572, 37145),
    F(117, 148580),
    F(13, 742900),
    F(1, 10400600),
]

numerator_arcsin = [
    F(1),
    F(0),
    -F(48612811065193217720717677178392447, 17521502337917614029717740270538300),
    F(0),
    F(1542880498363289248597622552896377147, 537326071696140163578010701629841200),
    F(0),
    -F(25629388878535757822653318119825725571, 18806412509364905725230374557044442000),
    F(0),
    F(59997700449506246047746748284324063301, 205817378502489528256921219152294373248),
    F(0),
    -F(3929590045745452902082508393999567898951, 171057110133180185706863413251017990210560),
    F(0),
    F(17667852791636166403806066845986733616137, 54939518901597871409380837432385778032332800),
    F(0),
]

denominator_arcsin = [
    F(1),
    F(0),
    -F(17177687151615384464112433518938499, 5840500779305871343239246756846100),
    F(0),
    F(70638889839100585257826495902581175, 21493042867845606543120428065193648),
    F(0),
    -F(260976433493184241308750702390067677, 150451300074919245801842996456355536),
    F(0),
    F(7106212846910135271370912693952576463, 16334712579562660972771525329547172480),
    F(0),
    -F(100234105411492207380201907281393735693, 2221520910820521892296927444818415457280),
    F(0),
    F(55355990129067361992384874772991295329, 44430418216410437845938548896368309145600),
    F(1),
]


### 2^n coefficients not used anymore?
pow2_positive=[
    2.0,                       # 2^{2^0} = 2^1
    4.0,                       # 2^2
    16.0,                      # 2^4
    256.0,                     # 2^8
    65536.0,                   # 2^{16}
    4294967296.0,              # 2^{32}
    1.8446744073709552e+19,   # 2^{64}
    3.402823669209385e+38,    # 2^{128}
    1.157920892373162e+77,    # 2^{256}
    1.3407807929942597e+154   # 2^{512}
]

pow2_negative=[
    0.5,                       # 2^{-1}
    0.25,                      # 2^{-2}
    0.0625,                    # 2^{-4}
    0.00390625,                # 2^{-8}
    1.52587890625e-05,         # 2^{-16}
    2.3283064365386963e-10,   # 2^{-32}
    5.421010862427522e-20,    # 2^{-64}
    2.938735877055719e-39,    # 2^{-128}
    8.636168555094445e-78,    # 2^{-256}
    7.458340731200207e-155    # 2^{-512}
]


### Pade for 2^x [13/13] not used anymore?
numerator_pow2=[
    1,
    0.34657359027997264,
    0.057654361670184166, 
    0.006105451953130373, 
    0.0004599974790604923, 
    2.6087396373446953e-05, 
    1.148089221866952e-06, 
    3.9789740358416966e-08, 
    1.0886899872887758e-09, 
    2.3290814666452394e-11, 
    3.7985794162342824e-13, 
    4.488024907663356e-15, 
    3.456513123366294e-17, 
    1.3164133659504377e-19
]

denumerator_pow2=[
    1, 
    -0.34657359027997264, 
    0.057654361670184166, 
    -0.006105451953130373, 
    0.0004599974790604923, 
    -2.6087396373446953e-05, 
    1.148089221866952e-06, 
    -3.9789740358416966e-08, 
    1.0886899872887758e-09, 
    -2.3290814666452394e-11, 
    3.7985794162342824e-13, 
    -4.488024907663356e-15, 
    3.456513123366294e-17, 
    -1.3164133659504377e-19
]

numerator_exp = [
    F(1),
    F(1, 2),
    F(3, 25),
    F(11, 600),
    F(11, 5520),
    F(3, 18400),
    F(1, 96600),
    F(1, 1932000),
    F(1, 48944000),
    F(1, 1585785600),
    F(1, 67395888000),
    F(1, 3953892096000),
    F(1, 355850288640000),
    F(1, 64764752532480000),
]

denumerator_exp = [
    F(1),
    -F(1, 2),
    F(3, 25),
    -F(11, 600),
    F(11, 5520),
    -F(3, 18400),
    F(1, 96600),
    -F(1, 1932000),
    F(1, 48944000),
    -F(1, 1585785600),
    F(1, 67395888000),
    -F(1, 3953892096000),
    F(1, 355850288640000),
    -F(1, 64764752532480000),
]

### Some other constants. ###
# Keep these as decimal -> exact rational. (Denominators will be huge; that's OK for indev accuracy checks.)
LN2 = F("0.693147180559945309417232121458176568")
PI = F("3.1415926535897932384626433832795028841971693993751058")
PI_2 = PI / 2
PI_4 = PI / 4
TAU = 2 * PI


# Padé approximants around x=1.
# Wolfram Alpha typically returns these in variable t = (x - 1) for a series/Padé about x=1.
# So evaluate them with t = x - 1.

# sqrt(x) Padé [13/13] about x=1
numerator_sqrt_13_13 = [
    F(1, 1),
    F(27, 4),
    F(81, 4),
    F(2277, 64),
    F(10395, 256),
    F(32319, 1024),
    F(8721, 512),
    F(26163, 4096),
    F(53703, 32768),
    F(36465, 131072),
    F(3861, 131072),
    F(7371, 4194304),
    F(819, 16777216),
    F(27, 67108864),
]

denominator_sqrt_13_13 = [
    F(1, 1),
    F(25, 4),
    F(69, 4),
    F(1771, 64),
    F(7315, 256),
    F(20349, 1024),
    F(4845, 512),
    F(12597, 4096),
    F(21879, 32768),
    F(12155, 131072),
    F(1001, 131072),
    F(1365, 4194304),
    F(91, 16777216),
    F(1, 67108864),
]

# cbrt(x) Padé [13/13] about x=1
numerator_cbrt_13_13 = [
    F(1, 1),
    F(20, 3),
    F(296, 15),
    F(13838, 405),
    F(214489, 5589),
    F(272986, 9315),
    F(779960, 50301),
    F(857956, 150903),
    F(214489, 150903),
    F(8579560, 36669429),
    F(2624336, 110008287),
    F(149110, 110008287),
    F(104377, 2970223749),
    F(2294, 8910671247),
]

denominator_cbrt_13_13 = [
    F(1, 1),
    F(19, 3),
    F(266, 15),
    F(11704, 405),
    F(169708, 5589),
    F(200564, 9315),
    F(114608, 10935),
    F(114608, 32805),
    F(25636, 32805),
    F(179452, 1594323),
    F(232232, 23914845),
    F(10556, 23914845),
    F(5278, 645700815),
    F(58, 1937102445),
]


def _pade_about_1(x, num_coeffs, den_coeffs):
    """Evaluate Padé approximant about x=1 using t = x-1."""
    t = _as_fraction(x) - 1
    return rational_eval(t, num_coeffs, den_coeffs)


def sqrt(x):
    """Square root via frexp/ldexp range reduction + Padé seed + 2 Newton iterations.

    Returns float.
    """
    xf = float(x)
    if xf < 0.0:
        raise ValueError("sqrt domain")
    if xf == 0.0:
        return 0.0

    m, e = frexp(xf)  # xf = m * 2**e, with m in [0.5, 1)

    # Make exponent even: if e is odd, pull one factor of 2 into the mantissa.
    # After this:
    # - e is even
    # - m is in [0.5, 1) if e was even, or in [1, 2) if e was odd
    if e & 1:
        m *= 2.0
        e -= 1

    # sqrt(xf) = sqrt(m) * 2**(e/2)
    y0 = float(_pade_about_1(m, numerator_sqrt_13_13, denominator_sqrt_13_13))

    # Newton-Raphson on y^2 - m = 0
    y = 0.5 * (y0 + m / y0)
    y = 0.5 * (y + m / y)

    return ldexp(y, e // 2)


def cbrt(x):
    """Cube root via frexp/ldexp range reduction + Padé seed + 2 Newton iterations.

    Note: reduction keeps mantissa in a limited range, but not as tight as sqrt.
    Returns float.
    """
    xf = float(x)
    if xf == 0.0:
        return 0.0

    sign = 1.0
    if xf < 0.0:
        sign = -1.0
        xf = -xf

    m, e = frexp(xf)  # m in [0.5, 1)

    # Choose q,r with e = 3q + r, r in {0,1,2}
    q = e // 3
    r = e - 3 * q

    # Fold remainder into mantissa: a = m * 2**r (a in [0.5, 4))
    a = ldexp(m, r)

    y0 = float(_pade_about_1(a, numerator_cbrt_13_13, denominator_cbrt_13_13))

    # Newton-Raphson on y^3 - a = 0: y_{n+1} = (2y + a/y^2)/3
    y = (2.0 * y0 + a / (y0 * y0)) / 3.0
    y = (2.0 * y + a / (y * y)) / 3.0

    return sign * ldexp(y, q)


numerator_sinc = [
     1.0,                                          
     0.0,                                          
    -0.15329450413674683,                          
     0.0,                                          
     0.006193549628917415,                         
     0.0,                                          
    -0.00010152460010152460,                       
     0.0,                                          
     7.840236846805112e-7,                         
     0.0,                                          
    -2.884389601750644e-9,                         
     0.0,                                          
     4.160167953514099e-12                         
]

denumerator_sinc = [
     1.0,                                          
     0.0,                                          
    13.361546804547310,                            
     0.0,                                          
    89.07819178082192,                             
     0.0,                                          
   385.0704225352113,                              
     0.0,                                          
   1160.228571428571,                              
     0.0,                                          
   2361.538461538461,                              
     0.0,                                          
   2609.538461538461                               
]


### polynomial and computation O(n) :(( ###
# Horner's method is already optimal O(n) for general polynomials.
# This implementation is type-generic: it works for float or Fraction coefficients.
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
    """Coerce inputs to Fraction for exact remainder arithmetic.

    Note: if x is a float, Fraction.from_float keeps the *binary64* value exactly.
    If you later want decimal-exact parsing, feed in strings and use Fraction("...").
    """
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, float):
        return Fraction.from_float(x)
    # Fallback: try Fraction constructor (works for str like "0.1")
    return Fraction(x)

def fpow2(n): ###Need to code how to compute blazingly fast. Let 2^n = 2^r * 2^k 
    global pow2_positive
    global pow2_negative
    global numerator_pow2
    global denumerator_pow2
    if n>0:
        r=0.01 #initial guess
        k= n-r
        while (r< 25 and r > -25) and (k <= 512 and k>=0):
            r+=1 ### wtf am I doing
            k=n-r #Need to somehow build a method to find closest value to k in set pow2_positive powers.
            #????
    

def EXP(x):
    global numerator_exp
    global denumerator_exp

    # Compute k in float-space (cheap), then compute remainder r exactly as a Fraction.
    # This keeps the Padé evaluation exact w.r.t. the chosen remainder.
    xf = float(x)
    k = int(round(xf / LN2_F))

    r = _as_fraction(xf) - k * LN2
    e_power_r = rational_eval(r, numerator_exp, denumerator_exp)  # Fraction

    # 2^k scaling using ldexp (later you can replace with your own exponent logic).
    return float(e_power_r) * ldexp(1.0, k)

def sinc(x):
    global numerator_sinc
    global denumerator_sinc #???
    if x==0:
        return 1
    else:
        return rational_eval(x, numerator_sinc, denumerator_sinc)

def _reduce_angle(x):
    """Reduce angle to [-pi, pi] exactly using Fraction."""
    xf = _as_fraction(x)
    # Modulo with Fraction preserves precision.
    # Result of % TAU is in [0, TAU)
    r = (xf + PI) % TAU
    return r - PI


def _sin_small(x):
    """sin(x) on small |x| using Taylor series with exact coefficients."""
    # sin(x) = x * (1 - x^2/3! + x^4/5! - ...)
    x2 = x * x
    
    # Coefficients:
    # 1/3! = 1/6
    # 1/5! = 1/120
    # 1/7! = 1/5040
    # 1/9! = 1/362880
    # 1/11! = 1/39916800
    # 1/13! = 1/6227020800
    # 1/15! = 1/1307674368000
    # 1/17! = 1/355687428096000
    
    f = F(1)
    f += x2 * -F(1, 6)
    f += x2 * x2 * F(1, 120)
    f += x2 * x2 * x2 * -F(1, 5040)
    f += x2 * x2 * x2 * x2 * F(1, 362880)
    f += x2 * x2 * x2 * x2 * x2 * -F(1, 39916800)
    f += x2 * x2 * x2 * x2 * x2 * x2 * F(1, 6227020800)
    f += x2 * x2 * x2 * x2 * x2 * x2 * x2 * -F(1, 1307674368000)
    f += x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * F(1, 355687428096000)
    
    return x * f


def _cos_small(x):
    """cos(x) on small |x| using Taylor series with exact coefficients."""
    # cos(x) = 1 - x^2/2! + x^4/4! - ...
    x2 = x * x
    
    # Coefficients:
    # 1/2! = 1/2
    # 1/4! = 1/24
    # 1/6! = 1/720
    # 1/8! = 1/40320
    # 1/10! = 1/3628800
    # 1/12! = 1/479001600
    # 1/14! = 1/87178291200
    # 1/16! = 1/20922789888000
    
    f = F(1)
    f += x2 * -F(1, 2)
    f += x2 * x2 * F(1, 24)
    f += x2 * x2 * x2 * -F(1, 720)
    f += x2 * x2 * x2 * x2 * F(1, 40320)
    f += x2 * x2 * x2 * x2 * x2 * -F(1, 3628800)
    f += x2 * x2 * x2 * x2 * x2 * x2 * F(1, 479001600)
    f += x2 * x2 * x2 * x2 * x2 * x2 * x2 * -F(1, 87178291200)
    f += x2 * x2 * x2 * x2 * x2 * x2 * x2 * x2 * F(1, 20922789888000)
    
    return f


def sin(x):
    # Perform reduction exactly
    xr = _reduce_angle(x)

    # Use symmetries to map to [0, pi]
    if xr < 0:
        return float(-sin(-xr)) # Recurse with positive reduced angle

    # Now xr in [0, pi]
    if xr > PI_2:
        # sin(x) = sin(pi - x)
        xr = PI - xr

    # Now xr in [0, pi/2]
    if xr <= PI_4:
        return float(_sin_small(xr))

    # sin(x) = cos(pi/2 - x)
    return float(_cos_small(PI_2 - xr))


def cos(x):
    # Perform reduction exactly
    xr = _reduce_angle(x)

    # cos is even
    if xr < 0:
        xr = -xr

    # Now xr in [0, pi]
    if xr > PI_2:
        # cos(x) = -cos(pi - x)
        return float(-cos(PI - xr)) # Recurse with simpler angle

    # Now xr in [0, pi/2]
    if xr <= PI_4:
        return float(_cos_small(xr))

    # cos(x) = sin(pi/2 - x)
    return float(_sin_small(PI_2 - xr))


_POLE_EPS = 1e-16


def _reduce_pi(x):
    """Reduce angle to [-pi/2, pi/2) by mod pi (exact Fraction)."""
    # x can be float/Fraction
    r = _reduce_angle(x)  # [-pi, pi] Fraction
    if r >= PI_2:
        r -= PI
    elif r < -PI_2:
        r += PI
    return r


def tan(x):
    # First, detect exact poles using 2pi reduction so tan(pi/2) gives +inf and tan(-pi/2) gives -inf.
    r2 = _reduce_angle(x)  # [-pi, pi] Fraction
    if abs(r2 - PI_2) < _POLE_EPS:
        return float("inf")
    if abs(r2 + PI_2) < _POLE_EPS:
        return float("-inf")

    # Use pi-reduced angle for regular evaluation (keeps values bounded, fixes sign flip across poles).
    r = _reduce_pi(x) # Fraction

    # Poles at r = +/- pi/2
    if abs(abs(r) - PI_2) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")

    s = sin(r)
    c = cos(r)
    if abs(c) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")
    return s / c


def cot(x):
    # Reduce by pi so pole at 0 represents k*pi.
    r = _reduce_pi(x) # Fraction

    # Poles at r = 0
    if abs(r) < _POLE_EPS:
        # match usual limit: cot(0+) = +inf, cot(0-) = -inf; choose +inf at exactly 0
        return float("inf")

    s = sin(r)
    c = cos(r)
    if abs(s) < _POLE_EPS:
        return float("inf") if r > 0 else float("-inf")
    return c / s


def atan2(y, x):
    """Quadrant-correct arctangent.

    Returns float in [-pi, pi].
    """
    yf = float(y)
    xf = float(x)

    if xf > 0.0:
        return atan(yf / xf)
    if xf < 0.0:
        if yf >= 0.0:
            return atan(yf / xf) + PI_F
        return atan(yf / xf) - PI_F

    # xf == 0
    if yf > 0.0:
        return PI_2_F
    if yf < 0.0:
        return -PI_2_F
    return 0.0

def _atan_series(z):
    """atan(z) for small |z| using odd-power Taylor series.

    atan(z) = z - z^3/3 + z^5/5 - ...
    Good when |z| is small (we enforce this via argument reduction).
    """
    z2 = z * z
    term = z
    s = term
    k = 3.0
    sign = -1.0
    # Up to z^25/25 is plenty for |z| <= ~0.3.
    for _ in range(1, 13):
        term *= z2
        s += sign * (term / k)
        sign = -sign
        k += 2.0
    return s


def atan(x):
    """atan(x) using argument reduction + Taylor series.

    Strategy:
      - handle sign
      - if |x| > 1: atan(x) = sign*(pi/2 - atan(1/|x|))
      - repeatedly apply half-angle reduction:
            atan(t) = 2*atan( t / (1 + sqrt(1+t^2)) )
        until |t| <= 0.3, then use series.

    Returns float.
    """
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
    # Reduce until t is small enough for the Taylor series.
    while t > 0.3:
        t = t / (1.0 + sqrt(1.0 + t * t))
        n += 1

    a = _atan_series(t)
    # Undo reductions: atan(x) = 2^n * atan(t)
    return sign * (a * (2.0 ** n))


def asin(x):
    """asin(x) via atan and sqrt.

    asin(x) = atan( x / sqrt(1 - x^2) ) for |x|<1.
    For |x|==1: +/-pi/2.

    Returns float.
    """
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

def sinh(x):
    return (0.5*(EXP(x) - EXP(-x)))

def cosh(x):
    return (0.5*(EXP(x) + EXP(-x)))

# atan is implemented above (argument reduction + series)

def _atanh_series(z):
    """atanh(z) for |z| <= ~1/3 using odd-power series.

    atanh(z) = z + z^3/3 + z^5/5 + ...
    """
    z2 = z * z
    term = z
    s = term
    k = 3.0
    # Up to z^23/23 is plenty for |z|<=1/3 in double.
    for _ in range(1, 12):
        term *= z2
        s += term / k
        k += 2.0
    return s


def ln1p(u):
    """Natural log of (1+u).

    Uses transform:
      ln(m) = 2 * atanh((m-1)/(m+1))
    which converges quickly for m near 1.

    Returns float.
    """
    uf = float(u)
    if uf <= -1.0:
        raise ValueError("ln1p domain")
    if uf == 0.0:
        return 0.0
    m = 1.0 + uf
    z = (m - 1.0) / (m + 1.0)
    return 2.0 * _atanh_series(z)


def ln(x):
    """Natural logarithm using frexp range reduction.

    x = m * 2**e with m in [0.5, 1).
    ln(x) = ln(m) + e*ln2, and ln(m) = 2*atanh((m-1)/(m+1)).

    Returns float.
    """
    xf = float(x)
    if xf <= 0.0:
        raise ValueError("Can't evaluate in complex set")

    m, e = frexp(xf)  # xf = m * 2**e, m in [0.5, 1)
    z = (m - 1.0) / (m + 1.0)
    return 2.0 * _atanh_series(z) + e * LN2_F


###Manual freexp inferior and crude.

#def frexp_manual(x):
    #if x <= 0:
    #    raise ValueError("x must be positive")
    #n = 0
    #m = x
    #while m >= 2.0:
    #    m /= 2.0
    #    n += 1
    #while m < 1.0:
    #    m *= 2.0
    #    n -= 1
#return m, n   


# --- Expression engine (tokenize -> shunting-yard -> eval) ---

class CalcError(Exception):
    pass


def _is_ident_start(ch):
    o = ord(ch)
    return (65 <= o <= 90) or (97 <= o <= 122) or ch == "_"


def _is_ident_continue(ch):
    o = ord(ch)
    return _is_ident_start(ch) or (48 <= o <= 57)


def tokenize(expr):
    """Convert expression string into tokens.

    Tokens are tuples: (type, value)
      - ("num", float)
      - ("id",  str)
      - ("op",  one of + - * / ^)
      - ("lpar", "(") / ("rpar", ")")
      - ("comma", ",")
    """
    tokens = []
    i = 0
    n = len(expr)
    while i < n:
        ch = expr[i]
        if ch.isspace():
            i += 1
            continue

        # number: 123, 1.23, .5, 1e-3
        if ch.isdigit() or ch == ".":
            j = i
            saw_dot = (ch == ".")
            i += 1
            while i < n:
                c = expr[i]
                if c.isdigit():
                    i += 1
                    continue
                if c == "." and not saw_dot:
                    saw_dot = True
                    i += 1
                    continue
                break
            # exponent part
            if i < n and expr[i] in ("e", "E"):
                k = i + 1
                if k < n and expr[k] in ("+", "-"):
                    k += 1
                k0 = k
                while k < n and expr[k].isdigit():
                    k += 1
                if k == k0:
                    raise CalcError(f"Invalid exponent at position {i}")
                i = k
            num_str = expr[j:i]
            try:
                tokens.append(("num", float(num_str)))
            except ValueError:
                raise CalcError(f"Invalid number '{num_str}'")
            continue

        if _is_ident_start(ch):
            j = i
            i += 1
            while i < n and _is_ident_continue(expr[i]):
                i += 1
            ident = expr[j:i]
            tokens.append(("id", ident))
            continue

        if ch in "+-*/^":
            tokens.append(("op", ch))
            i += 1
            continue

        if ch == "(":
            tokens.append(("lpar", ch))
            i += 1
            continue
        if ch == ")":
            tokens.append(("rpar", ch))
            i += 1
            continue
        if ch == ",":
            tokens.append(("comma", ch))
            i += 1
            continue

        raise CalcError(f"Unexpected character '{ch}' at position {i}")

    return tokens


# Operator precedence/associativity
_OP_INFO = {
    "+": (1, "L"),
    "-": (1, "L"),
    "*": (2, "L"),
    "/": (2, "L"),
    "^": (3, "R"),
    "u-": (4, "R"),  # unary minus
}


def to_rpn(tokens):
    """Shunting-yard algorithm with function calls.

    Produces tokens:
      - ("num", float)
      - ("var", name)
      - ("op", op)
      - ("func", name, argc)
    """
    output = []
    stack = []

    prev = None  # None | "value" | "op" | "lpar" | "comma"

    i = 0
    while i < len(tokens):
        ttype, val = tokens[i]

        if ttype == "num":
            output.append(("num", val))
            prev = "value"
            i += 1
            continue

        if ttype == "id":
            # If next token is '(', treat as function.
            if i + 1 < len(tokens) and tokens[i + 1][0] == "lpar":
                stack.append(("func", val, 0))  # argc counter (commas seen)
                prev = "op"
            else:
                output.append(("var", val))
                prev = "value"
            i += 1
            continue

        if ttype == "comma":
            # Pop until lpar
            while stack and stack[-1][0] != "lpar":
                output.append(stack.pop())
            if not stack:
                raise CalcError("Misplaced comma or missing '('")
            # increment nearest func argc
            j = len(stack) - 1
            while j >= 0 and stack[j][0] != "func":
                j -= 1
            if j < 0:
                raise CalcError("Comma used outside function call")
            ftype, fname, argc = stack[j]
            stack[j] = (ftype, fname, argc + 1)
            prev = "comma"
            i += 1
            continue

        if ttype == "op":
            op = val
            if op == "-" and (prev is None or prev in ("op", "lpar", "comma")):
                op = "u-"
            p1, a1 = _OP_INFO[op]
            while stack and stack[-1][0] == "op":
                top = stack[-1][1]
                p2, _ = _OP_INFO[top]
                if (a1 == "L" and p1 <= p2) or (a1 == "R" and p1 < p2):
                    output.append(stack.pop())
                else:
                    break
            stack.append(("op", op))
            prev = "op"
            i += 1
            continue

        if ttype == "lpar":
            stack.append(("lpar", val))
            prev = "lpar"
            i += 1
            continue

        if ttype == "rpar":
            while stack and stack[-1][0] != "lpar":
                output.append(stack.pop())
            if not stack:
                raise CalcError("Mismatched parentheses")
            stack.pop()  # discard lpar

            # If this closes a function call, emit it.
            if stack and stack[-1][0] == "func":
                _, fname, argc = stack.pop()
                # If previous token was lpar, this was an empty-arg call: f()
                if prev == "lpar":
                    argc_final = 0
                else:
                    argc_final = argc + 1
                output.append(("func", fname, argc_final))

            prev = "value"
            i += 1
            continue

        raise CalcError(f"Unhandled token type {ttype}")

    while stack:
        if stack[-1][0] == "lpar":
            raise CalcError("Mismatched parentheses")
        output.append(stack.pop())

    return output


_FUNCS = {
    "sin": (1, sin),
    "cos": (1, cos),
    "tan": (1, tan),
    "cot": (1, cot),
    "asin": (1, asin),
    "acos": (1, acos),
    "atan": (1, atan),
    "atan2": (2, atan2),
    "sqrt": (1, sqrt),
    "cbrt": (1, cbrt),
    "exp": (1, EXP),
    "ln": (1, ln),
    "log": (1, ln),
    "abs": (1, abs),
    "pow": (2, lambda a, b: a ** b),
}

_CONSTS = {
    "pi": PI_F,
    "tau": TAU_F,
    "e": E_F,
}


def eval_rpn(rpn, variables=None):
    if variables is None:
        variables = {}

    stack = []

    for tok in rpn:
        ttype = tok[0]

        if ttype == "num":
            stack.append(tok[1])
            continue

        if ttype == "var":
            name = tok[1]
            low = name.lower()
            if name in variables:
                stack.append(float(variables[name]))
            elif low in _CONSTS:
                stack.append(_CONSTS[low])
            else:
                raise CalcError(f"Unknown identifier '{name}'")
            continue

        if ttype == "op":
            op = tok[1]
            if op == "u-":
                if not stack:
                    raise CalcError("Missing operand for unary -")
                stack.append(-stack.pop())
                continue

            if len(stack) < 2:
                raise CalcError(f"Missing operand for operator '{op}'")
            b = stack.pop()
            a = stack.pop()
            if op == "+":
                stack.append(a + b)
            elif op == "-":
                stack.append(a - b)
            elif op == "*":
                stack.append(a * b)
            elif op == "/":
                stack.append(a / b)
            elif op == "^":
                stack.append(a ** b)
            else:
                raise CalcError(f"Unknown operator '{op}'")
            continue

        if ttype == "func":
            _, name, argc = tok
            low = name.lower()
            if low not in _FUNCS:
                raise CalcError(f"Unknown function '{name}'")
            expected, fn = _FUNCS[low]
            if argc != expected:
                raise CalcError(f"Function '{name}' expects {expected} arg(s), got {argc}")
            if len(stack) < argc:
                raise CalcError(f"Not enough values for function '{name}'")
            if argc == 0:
                res = fn()
            elif argc == 1:
                a = stack.pop()
                res = fn(a)
            elif argc == 2:
                b = stack.pop()
                a = stack.pop()
                res = fn(a, b)
            else:
                raise CalcError("Too many args (not supported)")
            stack.append(float(res))
            continue

        raise CalcError(f"Unhandled RPN token {tok}")

    if len(stack) != 1:
        raise CalcError("Invalid expression")
    return float(stack[0])


def _insert_implicit_multiplication(tokens):
    """Insert '*' tokens where multiplication is implied.

    Examples:
      2pi -> 2*pi
      (1+2)(3+4) -> (1+2)*(3+4)
      2sin(0.3) -> 2*sin(0.3)

    Rule (roughly):
      if prev is a value token and next starts a value/expression, insert '*',
      except for known function calls like sin( ... ).
    """
    if not tokens:
        return tokens

    out = [tokens[0]]

    def is_value(tok):
        return tok[0] in ("num", "id", "rpar")

    def starts_value(tok):
        return tok[0] in ("num", "id", "lpar")

    for i in range(1, len(tokens)):
        prev = out[-1]
        cur = tokens[i]

        if is_value(prev) and starts_value(cur):
            # Don't insert between function name and '(' for known functions.
            if prev[0] == "id" and cur[0] == "lpar":
                name = prev[1].lower()
                if name in _FUNCS:
                    out.append(cur)
                    continue
            out.append(("op", "*"))

        out.append(cur)

    return out


def evaluate(expr, variables=None):
    """High-level API for the GUI: expression string -> float."""
    toks = tokenize(expr)
    toks = _insert_implicit_multiplication(toks)
    rpn = to_rpn(toks)
    return eval_rpn(rpn, variables=variables)
