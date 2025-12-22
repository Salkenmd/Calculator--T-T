#include "math_core.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// --- Coefficients Embedded ---

const Real E_F_C = 2.7182818284590450907955982984276488L;
static const Real LN2_C = 0.69314718055994528622676398299518041L;
static const Real LN2_F_C = 0.69314718055994528622676398299518041L;
const Real PI_C = 3.1415926535897931159979634685441852L;
const Real PI_2_C = 1.5707963267948965579989817342720926L;
static const Real PI_2_F_C = 1.5707963267948965579989817342720926L;
const Real PI_4_C = 0.78539816339744827899949086713604629L;
static const Real PI_4_F_C = 0.78539816339744827899949086713604629L;
static const Real PI_F_C = 3.1415926535897931159979634685441852L;
static const Real TAU_C = 6.2831853071795862319959269370883703L;
static const Real TAU_F_C = 6.2831853071795862319959269370883703L;

static const Real coeffs_atan[] = {
    1L, -0.33333333333333331482961625624739099L, 0.2000000000000000111022302462515654L,
    -0.14285714285714284921269268124888185L, 0.11111111111111110494320541874913033L,
    -0.090909090909090911614143237784446683L, 0.076923076923076927347011633173679002L,
    -0.06666666666666666574148081281236955L, 0.05882352941176470506601248189326725L,
    -0.052631578947368418130992040460114367L, 0.047619047619047616404230893749627285L,
    -0.043478260869565216184540190624829847L, 0.040000000000000000832667268468867405L,
};
static const int coeffs_atan_len = 13;

static const Real coeffs_atanh[] = {
    1L, 0.33333333333333331482961625624739099L, 0.2000000000000000111022302462515654L,
    0.14285714285714284921269268124888185L, 0.11111111111111110494320541874913033L,
    0.090909090909090911614143237784446683L, 0.076923076923076927347011633173679002L,
    0.06666666666666666574148081281236955L, 0.05882352941176470506601248189326725L,
    0.052631578947368418130992040460114367L, 0.047619047619047616404230893749627285L,
    0.043478260869565216184540190624829847L,
};
static const int coeffs_atanh_len = 12;

static const Real coeffs_cos[] = {
    1L, -0.5L, 0.041666666666666664353702032030923874L,
    -0.0013888888888888889418943284326246612L, 2.480158730158730156578963943481142e-05L,
    -2.7557319223985888275785856998956191e-07L, 2.0876756987868100186555149434461387e-09L,
    -1.1470745597729724507296570450374645e-11L, 4.7794773323873852534461526327767315e-14L,
};
static const int coeffs_cos_len = 9;

static const Real coeffs_sin[] = {
    1L, -0.1666666666666666574148081281236955L, 0.0083333333333333332176851016015461937L,
    -0.00019841269841269841252631711547849136L, 2.7557319223985892510950593270457887e-06L,
    -2.5052108385441720223866179321353664e-08L, 1.6059043836821613340862918294945196e-10L,
    -7.6471637318198164055138442124427704e-13L, 2.8114572543455205981105182743959941e-15L,
};
static const int coeffs_sin_len = 9;

static const Real denominator_LN1P[] = {
    1L, 6.5L, 18.719999999999998863131622783839703L, 31.460000000000000852651282912120223L,
    34.195652173913046567577112000435591L, 2.5180434782608696053785024560056627L,
    12.790062111801242750175333640072495L, 4.4765217391304350513792087440378964L,
    1.0602288329519451437477073341142386L, 0.16361556064073226601784938338823849L,
    0.015399111589715977488346254631323973L, 0.0007874545699286579101758287713153095L,
    1.7498990442859066251702140903923066e-05L, 9.6148299136588279162458860107459691e-08L,
};
static const int denominator_LN1P_len = 14;

static const Real denominator_arcsin[] = {
    1L, 0L, -2.9411325844659690353921632777201012L, 0L, 3.2865932605931287291411990736378357L,
    0L, -1.7346239837291370644578591964091174L, 0L, 0.43503751977865495748076796189707238L,
    0L, -0.045119586731448141359024361918272916L, 0L, 0.0012459029725860547820059931467540082L, 1L,
};
static const int denominator_arcsin_len = 14;

static const Real denominator_cbrt_13_13[] = {
    1L, 6.3333333333333330372738600999582559L, 17.733333333333334280723647680133581L,
    28.898765432098766581248128204606473L, 30.364644838074788424364669481292367L,
    21.531293612453033858855633297935128L, 10.480841335162322636165299627464265L,
    3.4936137783874410267515031591756269L, 0.7814662399024538963487884757341817L,
    0.11255686582957154240336450357062859L, 0.0097107884245120554278507540857390268L,
    0.00044139947384145704518954023143351151L, 8.1740643303973526886951894709909539e-06L,
    2.9941627583873086383941652249829346e-08L,
};
static const int denominator_cbrt_13_13_len = 14;

static const Real denominator_sqrt_13_13[] = {
    1L, 6.25L, 17.25L, 27.671875L, 28.57421875L, 19.8720703125L, 9.462890625L,
    3.075439453125L, 0.667694091796875L, 0.09273529052734375L, 0.00763702392578125L,
    0.0003254413604736328125L, 5.424022674560546875e-06L, 1.490116119384765625e-08L,
};
static const int denominator_sqrt_13_13_len = 14;

static const Real denumerator_exp[] = {
    1L, -0.5L, 0.11999999999999999555910790149937384L, -0.018333333333333333425851918718763045L,
    0.0019927536231884057128793674706912498L, -0.00016304347826086957803926047461118287L,
    1.0351966873706003395257460586442733e-05L, -5.1759834368530021211452039203715358e-07L,
    2.04315135665250075022251601561743e-08L, -6.3060227057175949591547168724629119e-10L,
    1.4837700484041399485741487398290288e-11L, -2.5291534915979658214332080792540263e-13L,
    2.8101705462199623130636742246532468e-15L, -1.5440497506703087532742710210031064e-17L,
};
static const int denumerator_exp_len = 14;

static const Real denumerator_sinc[] = {
    1L, 0L, 13.361546804547309363897511502727866L, 0L, 89.078191780821924794508959166705608L,
    0L, 385.07042253521132124660653062164783L, 0L, 1160.2285714285710582771571353077888L,
    0L, 2361.5384615384609787724912166595459L, 0L, 2609.5384615384609787724912166595459L,
};
static const int denumerator_sinc_len = 13;

static const Real numerator_LN1P[] = {
    1L, 6L, 15.803333333333332788583902583923191L, 24.016666666666665719276352319866419L,
    23.280652173913043867514716112054884L, 15.022608695652174048973392928019166L,
    6.5369192546583851921582208888139576L, 1.9016273291925465116491977823898196L,
    0.35890723184773526410751287585298996L, 0.041479059968762485577364884647977306L,
    0.0026460580071429443994868169198753094L, 7.630027949048006036085467496121737e-05L,
    6.1152890316592413428979386341710622e-07L, 0L,
};
static const int numerator_LN1P_len = 14;

static const Real numerator_arcsin[] = {
    1L, 0L, -2.7744659177993025167552332277409732L, 0L, 2.8714044965154674748930574423866346L,
    0L, -1.3628005269890395290133255912223831L, 0L, 0.29150939967287808896045930850959849L,
    0L, -0.022972386489436107115391649813318509L, 0L, 0.00032158732265713947358690139033399191L, 0L,
};
static const int numerator_arcsin_len = 14;

static const Real numerator_cbrt_13_13[] = {
    1L, 6.6666666666666669627261399000417441L, 19.733333333333334280723647680133581L,
    34.167901234567899848570959875360131L, 38.37699051708713682273810263723135L,
    29.306065485775629753106841235421598L, 15.505854754378640691925284045282751L,
    5.6854800766055015870392708166036755L, 1.4213700191513753967598177041509189L,
    0.23397037352286015332580859649169724L, 0.023855802790566132259852949459855154L,
    0.0013554433403730756571570204016552452L, 3.514112363930196343943593295655603e-05L,
    2.5744412922565540418910127844209779e-07L,
};
static const int numerator_cbrt_13_13_len = 14;

static const Real numerator_exp[] = {
    1L, 0.5L, 0.11999999999999999555910790149937384L, 0.018333333333333333425851918718763045L,
    0.0019927536231884057128793674706912498L, 0.00016304347826086957803926047461118287L,
    1.0351966873706003395257460586442733e-05L, 5.1759834368530021211452039203715358e-07L,
    2.04315135665250075022251601561743e-08L, 6.3060227057175949591547168724629119e-10L,
    1.4837700484041399485741487398290288e-11L, 2.5291534915979658214332080792540263e-13L,
    2.8101705462199623130636742246532468e-15L, 1.5440497506703087532742710210031064e-17L,
};
static const int numerator_exp_len = 14;

static const Real numerator_sinc[] = {
    1L, 0L, -0.15329450413674683417752930836286396L, 0L, 0.0061935496289174151243739352423745004L,
    0L, -0.00010152460010152459688945181248698191L, 0L, 7.8402368468051121814507481086464757e-07L,
    0L, -2.8843896017506438908679855732172109e-09L, 0L, 4.1601679535140991798873542454890572e-12L,
};
static const int numerator_sinc_len = 13;

static const Real numerator_sqrt_13_13[] = {
    1L, 6.75L, 20.25L, 35.578125L, 40.60546875L, 31.5615234375L, 17.033203125L, 6.387451171875L,
    1.638885498046875L, 0.27820587158203125L, 0.02945709228515625L, 0.0017573833465576171875L,
    4.8816204071044921875e-05L, 4.0233135223388671875e-07L,
};
static const int numerator_sqrt_13_13_len = 14;

static const Real two_over_pi[] = {
    10680707.0L, 7228996.0L, 1387004.0L, 2578385.0L, 16069853.0L, 12639074.0L, 9804092.0L,
    4427841.0L, 16666979.0L, 11263675.0L, 12935607.0L, 2387514.0L, 4345298.0L, 14681673.0L,
    3074569.0L, 13734428.0L, 16653803.0L, 1880361.0L, 10960616.0L, 8533493.0L, 3062596.0L,
    8710556.0L, 7349940.0L, 6258241.0L, 3772886.0L, 3769171.0L, 3798172.0L, 8675211.0L,
    12450088.0L, 3874808.0L, 9961438.0L, 366607.0L, 15675153.0L, 9132554.0L, 7151469.0L,
    3571407.0L, 2607881.0L, 12013382.0L, 4155038.0L, 6285869.0L, 7677882.0L, 13102053.0L,
    15825725.0L, 473591.0L, 9065106.0L, 15363067.0L, 6271263.0L, 9264392.0L, 5636912.0L,
    4652155.0L, 7056368.0L, 13614112.0L, 10155062.0L, 1944035.0L, 9527646.0L, 15080200.0L,
    6658437.0L, 6231200.0L, 6832269.0L, 16767104.0L, 5075751.0L, 3212806.0L, 1398474.0L,
    7579849.0L, 6349435.0L, 12618859.0L, 4703257.0L, 12806093.0L, 14477321.0L, 2786137.0L,
    12875403.0L, 9837734.0L, 14528324.0L, 13719321.0L, 343717.0L, 16713477.0L, 4161075.0L,
    15217346.0L, 14569368.0L, 3308987.0L, 12795174.0L, 15690526.0L, 6224031.0L, 3809077.0L,
    13300351.0L, 1935345.0L, 2199676.0L, 8135786.0L, 16412373.0L, 7810352.0L, 4406037.0L,
    12981429.0L, 10295747.0L, 12764333.0L, 4279596.0L, 6094860.0L, 4619654.0L, 2978275.0L,
    10143387.0L, 25139.0L, 8180404.0L, 9938868.0L, 13980983.0L, 16137943.0L, 1577123.0L,
    16545357.0L, 2792804.0L, 11261808.0L, 16284771.0L, 5746810.0L, 15144215.0L, 5654976.0L,
    14276155.0L, 3703975.0L, 13312804.0L, 7834326.0L, 2315354.0L, 12132096.0L, 1772273.0L,
    14667289.0L, 16724383.0L, 6954598.0L, 6379417.0L, 4717484.0L, 14188414.0L, 12018978.0L,
    9037874.0L, 6340582.0L, 13485295.0L, 603756.0L, 13909853.0L, 14147094.0L, 14564184.0L,
    9608158.0L, 2630354.0L, 15238696.0L, 5069026.0L, 3328710.0L, 1499912.0L, 13336032.0L,
    5292055.0L, 10952179.0L, 6021144.0L, 3412782.0L, 6427267.0L, 84099.0L, 6000373.0L,
    8368301.0L, 15919390.0L, 4409928.0L,
};

// --- Core Logic ---

Real poly_eval(Real x, const Real* coeffs, int n) {
    Real y = 0;
    int i = n - 1;
    while (i >= 0) {
        y = y * x + coeffs[i];
        i--;
    }
    return y;
}

Real rational_eval(Real x, const Real* num, int num_len, const Real* den, int den_len) {
    Real p = poly_eval(x, num, num_len);
    Real q = poly_eval(x, den, den_len);
    return p / q;
}

static Real _pade_about_1(Real x, const Real* num, int num_len, const Real* den, int den_len) {
    Real t = x - 1.0L;
    return rational_eval(t, num, num_len, den, den_len);
}

int rem_pio2l(Real x, Real *y) {
    if (isnan(x) || isinf(x)) {
        *y = NAN;
        return 0;
    }
    
    Real ax = fabsl(x);
    if (ax <= PI_4_C) {
        *y = x;
        return 0;
    }
    
    if (ax < 3 * PI_4_C) {
        if (x > 0) {
            *y = x - PI_2_C;
            return 1;
        } else {
            *y = x + PI_2_C;
            return -1;
        }
    }
    
    int e_int;
    Real m = frexpl(ax, &e_int);
    int jk = e_int / 24;
    int k_start = jk - 2;
    if (k_start < 0) k_start = 0;
    
    Real accum_frac = 0.0L;
    int accum_quad = 0;
    
    for (int k = k_start; k < k_start + 5; k++) {
        Real term = m * two_over_pi[k]; 
        int shift = e_int - 24 * (k + 1);
        Real val = ldexpl(term, shift);
        Real ipart;
        Real fpart = modfl(val, &ipart);
        long long i_ll = (long long)fmodl(ipart, 16.0L); 
        accum_quad += i_ll;
        accum_frac += fpart;
    }
    
    Real f_ipart;
    Real f_fpart = modfl(accum_frac, &f_ipart);
    accum_quad += (long long)f_ipart;
    accum_frac = f_fpart;
    
    int n = accum_quad & 3;
    
    if (accum_frac >= 0.5L) {
        n += 1;
        accum_frac -= 1.0L;
    } else if (accum_frac < -0.5L) {
        n -= 1;
        accum_frac += 1.0L;
    }
    
    n &= 3;
    *y = accum_frac * PI_2_C;
    
    if (x < 0) {
        n = -n;
        *y = -(*y);
    }
    
    return n & 3; 
}

Real mc_sqrt(Real x) {
    if (x < 0) return NAN;
    if (x == 0) return 0;
    
    int e;
    Real m = frexpl(x, &e);
    if (e & 1) {
        m *= 2.0L;
        e--;
    }
    Real y0 = _pade_about_1(m, numerator_sqrt_13_13, numerator_sqrt_13_13_len, 
                            denominator_sqrt_13_13, denominator_sqrt_13_13_len);
    y0 = 0.5L * (y0 + m / y0);
    y0 = 0.5L * (y0 + m / y0);
    
    return ldexpl(y0, e / 2);
}

Real mc_cbrt(Real x) {
    if (x == 0) return 0;
    int sign = 1;
    if (x < 0) { sign = -1; x = -x; }
    
    int e;
    Real m = frexpl(x, &e);
    int q = e / 3;
    int r = e - 3 * q;
    if (r < 0) { 
        q--;
        r += 3;
    }
    
    Real a = ldexpl(m, r);
    Real y0 = _pade_about_1(a, numerator_cbrt_13_13, numerator_cbrt_13_13_len,
                             denominator_cbrt_13_13, denominator_cbrt_13_13_len);
                             
    y0 = (2.0L * y0 + a / (y0 * y0)) / 3.0L;
    y0 = (2.0L * y0 + a / (y0 * y0)) / 3.0L;
    
    return sign * ldexpl(y0, q);
}

Real mc_exp(Real x) {
    Real k_float = roundl(x / LN2_F_C);
    int k = (int)k_float;
    Real r = x - k * LN2_C; 
    
    Real e_r = rational_eval(r, numerator_exp, numerator_exp_len, denumerator_exp, denumerator_exp_len);
    return e_r * ldexpl(1.0L, k);
}

static Real _sin_small(Real x) {
    Real x2 = x*x;
    Real f = poly_eval(x2, coeffs_sin, coeffs_sin_len);
    return x * f;
}
static Real _cos_small(Real x) {
    Real x2 = x*x;
    return poly_eval(x2, coeffs_cos, coeffs_cos_len);
}

Real mc_sin(Real x) {
    Real y;
    int n = rem_pio2l(x, &y);
    int quadrant = n & 3;
    
    switch(quadrant) {
        case 0: return _sin_small(y);
        case 1: return _cos_small(y);
        case 2: return -_sin_small(y);
        case 3: return -_cos_small(y);
    }
    return 0;
}

Real mc_cos(Real x) {
    Real y;
    int n = rem_pio2l(x, &y);
    int quadrant = n & 3;
    
    switch(quadrant) {
        case 0: return _cos_small(y);
        case 1: return -_sin_small(y);
        case 2: return -_cos_small(y);
        case 3: return _sin_small(y);
    }
    return 0;
}

Real mc_tan(Real x) {
    Real y;
    int n = rem_pio2l(x, &y);
    Real s = _sin_small(y);
    Real c = _cos_small(y);
    int odd = n & 1;
    
    if (odd) {
        return -c/s;
    } else {
        return s/c;
    }
}

Real mc_cot(Real x) {
    Real y;
    int n = rem_pio2l(x, &y);
    Real s = _sin_small(y);
    Real c = _cos_small(y);
    int odd = n & 1;
    
    if (odd) {
        return -s/c;
    } else {
        return c/s;
    }
}

static Real _atan_series(Real z) {
    Real z2 = z*z;
    return z * poly_eval(z2, coeffs_atan, coeffs_atan_len);
}

Real mc_atan(Real x) {
    Real xf = x;
    if (xf == 0) return 0;
    int sign = 1;
    if (xf < 0) { sign = -1; xf = -xf; }
    
    if (xf > 1.0L) {
        return sign * (PI_2_C - mc_atan(1.0L/xf));
    }
    
    int n = 0;
    while (xf > 0.3L) {
        xf = xf / (1.0L + mc_sqrt(1.0L + xf*xf));
        n++;
    }
    Real a = _atan_series(xf);
    return sign * (a * ldexpl(1.0L, n));
}

Real mc_atan2(Real y, Real x) {
    if (x > 0) return mc_atan(y/x);
    if (x < 0) {
        if (y >= 0) return mc_atan(y/x) + PI_C;
        return mc_atan(y/x) - PI_C;
    }
    if (y > 0) return PI_2_C;
    if (y < 0) return -PI_2_C;
    return 0;
}

Real mc_asin(Real x) {
    if (fabsl(x) > 1) return NAN;
    if (x == 1) return PI_2_C;
    if (x == -1) return -PI_2_C;
    return mc_atan(x / mc_sqrt(1.0L - x*x));
}

Real mc_acos(Real x) {
    return PI_2_C - mc_asin(x);
}

static Real _atanh_series(Real z) {
    Real z2 = z * z;
    return z * poly_eval(z2, coeffs_atanh, coeffs_atanh_len);
}

Real mc_ln1p(Real u) {
    if (u <= -1) return NAN;
    if (u == 0) return 0;
    Real m = 1.0L + u;
    Real z = (m - 1.0L) / (m + 1.0L);
    return 2.0L * _atanh_series(z);
}

Real mc_ln(Real x) {
    if (x <= 0) return NAN;
    int e;
    Real m = frexpl(x, &e);
    Real z = (m - 1.0L) / (m + 1.0L);
    return 2.0L * _atanh_series(z) + e * LN2_C;
}

Real mc_sinh(Real x) {
    return 0.5L * (mc_exp(x) - mc_exp(-x));
}

Real mc_cosh(Real x) {
    return 0.5L * (mc_exp(x) + mc_exp(-x));
}

Real mc_sinc(Real x) {
    if (x == 0) return 1.0L;
    return rational_eval(x, numerator_sinc, numerator_sinc_len, denumerator_sinc, denumerator_sinc_len);
}
