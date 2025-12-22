#include "math_core.h"

#include <stdint.h>

#ifndef NAN
#define NAN (0.0/0.0)
#endif
#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

typedef union {
    double f;
    uint64_t u;
} DoubleU;

#define DBL_EXP_MASK    0x7FF0000000000000ULL
#define DBL_MANT_MASK   0x000FFFFFFFFFFFFFULL
#define DBL_EXP_SHIFT   52
#define DBL_EXP_BIAS    1023
#define DBL_SIGN_MASK   0x8000000000000000ULL

Real mc_fabs(Real x) {
    if (x < 0) return -x;
    return x;
}

int mc_isnan(Real x) {
    DoubleU val;
    val.f = x;
    uint64_t e = val.u & DBL_EXP_MASK;
    uint64_t m = val.u & DBL_MANT_MASK;
    return (e == DBL_EXP_MASK) && (m != 0);
}

int mc_isinf(Real x) {
    DoubleU val;
    val.f = x;
    uint64_t e = val.u & DBL_EXP_MASK;
    uint64_t m = val.u & DBL_MANT_MASK;
    return (e == DBL_EXP_MASK) && (m == 0);
}

Real mc_frexp(Real value, int *exp) {
    DoubleU val;
    val.f = value;
    
    if (value == 0.0) {
        *exp = 0;
        return 0.0;
    }
    
    if (mc_isnan(value) || mc_isinf(value)) {
        *exp = 0;
        return value;
    }

    uint64_t raw_exp = (val.u & DBL_EXP_MASK) >> DBL_EXP_SHIFT;
    
    if (raw_exp == 0) {
        int sub_exp;
        Real norm_val = value * 18014398509481984.0;
        Real ret = mc_frexp(norm_val, &sub_exp);
        *exp = sub_exp - 54;
        return ret;
    }

    int true_exp = (int)raw_exp - DBL_EXP_BIAS;
    *exp = true_exp + 1;
    
    uint64_t new_exp = 0x3FEULL;
    val.u = (val.u & ~DBL_EXP_MASK) | (new_exp << DBL_EXP_SHIFT);
    
    return val.f;
}

Real mc_ldexp(Real value, int exp) {
    if (value == 0.0) return 0.0;
    if (mc_isnan(value) || mc_isinf(value)) return value;

    int current_exp;
    Real m = mc_frexp(value, &current_exp);
    
    long long new_exp_unbiased = (long long)current_exp + exp;
    
    long long biased_exp = new_exp_unbiased + 1022;
    
    DoubleU u;
    
    if (biased_exp >= 2047) {
        u.f = value;
        u.u = (u.u & DBL_SIGN_MASK) | DBL_EXP_MASK;
        return u.f;
    }
    
    if (biased_exp <= 0) {
        long long shift = 1 - biased_exp;
        if (shift > 53) {
            u.f = value;
            u.u &= DBL_SIGN_MASK;
            return u.f;
        }
        
        u.f = m;
        uint64_t mant = (u.u & DBL_MANT_MASK) | (1ULL << 52);
        mant >>= shift;
        
        DoubleU sign_u; 
        sign_u.f = value;
        u.u = (sign_u.u & DBL_SIGN_MASK) | mant;
        return u.f;
    }
    
    u.f = m;
    uint64_t mant = u.u & DBL_MANT_MASK;
    
    DoubleU res;
    DoubleU sign_u;
    sign_u.f = value;
    
    res.u = (sign_u.u & DBL_SIGN_MASK) | ((uint64_t)biased_exp << DBL_EXP_SHIFT) | mant;
    return res.f;
}

Real mc_modf(Real x, Real *iptr) {
    if (mc_isnan(x)) { *iptr = x; return x; }
    if (mc_isinf(x)) { *iptr = x; return 0.0; }
    if (x == 0.0) { *iptr = x; return x; }

    DoubleU u;
    u.f = x;
    uint64_t e_bits = (u.u & DBL_EXP_MASK) >> DBL_EXP_SHIFT;
    int e = (int)e_bits - DBL_EXP_BIAS;

    if (e < 0) {
        u.u &= DBL_SIGN_MASK;
        *iptr = u.f;
        return x;
    }
    if (e >= 52) {
        *iptr = x;
        u.u &= DBL_SIGN_MASK;
        return u.f;
    }

    uint64_t mask = -1ULL;
    mask <<= (52 - e);

    DoubleU i_val;
    i_val.u = u.u & mask;
    *iptr = i_val.f;

    return x - *iptr;
}

Real mc_fmod(Real x, Real y) {
    if (y == 0.0) return NAN; 
    Real q = x / y;
    Real iptr;
    mc_modf(q, &iptr);
    return x - iptr * y;
}

Real mc_floor(Real x) {
    if (mc_isnan(x) || mc_isinf(x)) return x;
    Real iptr;
    Real frac = mc_modf(x, &iptr);
    if (frac < 0.0) return iptr - 1.0;
    return iptr;
}

Real mc_ceil(Real x) {
    if (mc_isnan(x) || mc_isinf(x)) return x;
    Real iptr;
    Real frac = mc_modf(x, &iptr);
    if (frac > 0.0) return iptr + 1.0;
    return iptr;
}

Real mc_round(Real x) {
    if (mc_isnan(x) || mc_isinf(x)) return x;
    Real iptr;
    Real frac = mc_modf(x, &iptr);
    Real abs_frac = mc_fabs(frac);

    if (abs_frac > 0.5) {
        if (x > 0) return iptr + 1.0;
        return iptr - 1.0;
    }
    if (abs_frac == 0.5) {
        if (x > 0) return iptr + 1.0;
        return iptr - 1.0;
    }
    return iptr;
}

// --- Coefficients Embedded ---

const Real E_F_C = 2.7182818284590450907955982984276488;
static const Real LN2_C = 0.69314718055994528622676398299518041;
static const Real LN2_F_C = 0.69314718055994528622676398299518041;
const Real PI_C = 3.1415926535897931159979634685441852;
const Real PI_2_C = 1.5707963267948965579989817342720926;
static const Real PI_2_F_C = 1.5707963267948965579989817342720926;
const Real PI_4_C = 0.78539816339744827899949086713604629;
static const Real PI_4_F_C = 0.78539816339744827899949086713604629;
static const Real PI_F_C = 3.1415926535897931159979634685441852;
static const Real TAU_C = 6.2831853071795862319959269370883703;
static const Real TAU_F_C = 6.2831853071795862319959269370883703;

static const Real coeffs_atan[] = {
    1, -0.33333333333333331482961625624739099, 0.2000000000000000111022302462515654,
    -0.14285714285714284921269268124888185, 0.11111111111111110494320541874913033,
    -0.090909090909090911614143237784446683, 0.076923076923076927347011633173679002,
    -0.06666666666666666574148081281236955, 0.05882352941176470506601248189326725,
    -0.052631578947368418130992040460114367, 0.047619047619047616404230893749627285,
    -0.043478260869565216184540190624829847, 0.040000000000000000832667268468867405,
};
static const int coeffs_atan_len = 13;

static const Real coeffs_atanh[] = {
    1, 0.33333333333333331482961625624739099, 0.2000000000000000111022302462515654,
    0.14285714285714284921269268124888185, 0.11111111111111110494320541874913033,
    0.090909090909090911614143237784446683, 0.076923076923076927347011633173679002,
    0.06666666666666666574148081281236955, 0.05882352941176470506601248189326725,
    0.052631578947368418130992040460114367, 0.047619047619047616404230893749627285,
    0.043478260869565216184540190624829847,
};
static const int coeffs_atanh_len = 12;

static const Real coeffs_cos[] = {
    1, -0.5, 0.041666666666666664353702032030923874,
    -0.0013888888888888889418943284326246612, 2.480158730158730156578963943481142e-05,
    -2.7557319223985888275785856998956191e-07, 2.0876756987868100186555149434461387e-09,
    -1.1470745597729724507296570450374645e-11, 4.7794773323873852534461526327767315e-14,
};
static const int coeffs_cos_len = 9;

static const Real coeffs_sin[] = {
    1, -0.1666666666666666574148081281236955, 0.0083333333333333332176851016015461937,
    -0.00019841269841269841252631711547849136, 2.7557319223985892510950593270457887e-06,
    -2.5052108385441720223866179321353664e-08, 1.6059043836821613340862918294945196e-10,
    -7.6471637318198164055138442124427704e-13, 2.8114572543455205981105182743959941e-15,
};
static const int coeffs_sin_len = 9;

static const Real denominator_LN1P[] = {
    1, 6.5, 18.719999999999998863131622783839703, 31.460000000000000852651282912120223,
    34.195652173913046567577112000435591, 2.5180434782608696053785024560056627,
    12.790062111801242750175333640072495, 4.4765217391304350513792087440378964,
    1.0602288329519451437477073341142386, 0.16361556064073226601784938338823849,
    0.015399111589715977488346254631323973, 0.0007874545699286579101758287713153095,
    1.7498990442859066251702140903923066e-05, 9.6148299136588279162458860107459691e-08,
};
static const int denominator_LN1P_len = 14;

static const Real denominator_arcsin[] = {
    1, 0, -2.9411325844659690353921632777201012, 0, 3.2865932605931287291411990736378357,
    0, -1.7346239837291370644578591964091174, 0, 0.43503751977865495748076796189707238,
    0, -0.045119586731448141359024361918272916, 0, 0.0012459029725860547820059931467540082, 1,
};
static const int denominator_arcsin_len = 14;

static const Real denominator_cbrt_13_13[] = {
    1, 6.3333333333333330372738600999582559, 17.733333333333334280723647680133581,
    28.898765432098766581248128204606473, 30.364644838074788424364669481292367,
    21.531293612453033858855633297935128, 10.480841335162322636165299627464265,
    3.4936137783874410267515031591756269, 0.7814662399024538963487884757341817,
    0.11255686582957154240336450357062859, 0.0097107884245120554278507540857390268,
    0.00044139947384145704518954023143351151, 8.1740643303973526886951894709909539e-06,
    2.9941627583873086383941652249829346e-08,
};
static const int denominator_cbrt_13_13_len = 14;

static const Real denominator_sqrt_13_13[] = {
    1, 6.25, 17.25, 27.671875, 28.57421875, 19.8720703125, 9.462890625,
    3.075439453125, 0.667694091796875, 0.09273529052734375, 0.00763702392578125,
    0.0003254413604736328125, 5.424022674560546875e-06, 1.490116119384765625e-08,
};
static const int denominator_sqrt_13_13_len = 14;

static const Real denumerator_exp[] = {
    1, -0.5, 0.11999999999999999555910790149937384, -0.018333333333333333425851918718763045,
    0.0019927536231884057128793674706912498, -0.00016304347826086957803926047461118287,
    1.0351966873706003395257460586442733e-05, -5.1759834368530021211452039203715358e-07,
    2.04315135665250075022251601561743e-08, -6.3060227057175949591547168724629119e-10,
    1.4837700484041399485741487398290288e-11, -2.5291534915979658214332080792540263e-13,
    2.8101705462199623130636742246532468e-15, -1.5440497506703087532742710210031064e-17,
};
static const int denumerator_exp_len = 14;

static const Real denumerator_sinc[] = {
    1, 0, 13.361546804547309363897511502727866, 0, 89.078191780821924794508959166705608,
    0, 385.07042253521132124660653062164783, 0, 1160.2285714285710582771571353077888,
    0, 2361.5384615384609787724912166595459, 0, 2609.5384615384609787724912166595459,
};
static const int denumerator_sinc_len = 13;

static const Real numerator_LN1P[] = {
    1, 6, 15.803333333333332788583902583923191, 24.016666666666665719276352319866419,
    23.280652173913043867514716112054884, 15.022608695652174048973392928019166,
    6.5369192546583851921582208888139576, 1.9016273291925465116491977823898196,
    0.35890723184773526410751287585298996, 0.041479059968762485577364884647977306,
    0.0026460580071429443994868169198753094, 7.630027949048006036085467496121737e-05,
    6.1152890316592413428979386341710622e-07, 0,
};
static const int numerator_LN1P_len = 14;

static const Real numerator_arcsin[] = {
    1, 0, -2.7744659177993025167552332277409732, 0, 2.8714044965154674748930574423866346,
    0, -1.3628005269890395290133255912223831, 0, 0.29150939967287808896045930850959849,
    0, -0.022972386489436107115391649813318509, 0, 0.00032158732265713947358690139033399191, 0,
};
static const int numerator_arcsin_len = 14;

static const Real numerator_cbrt_13_13[] = {
    1, 6.6666666666666669627261399000417441, 19.733333333333334280723647680133581,
    34.167901234567899848570959875360131, 38.37699051708713682273810263723135,
    29.306065485775629753106841235421598, 15.505854754378640691925284045282751,
    5.6854800766055015870392708166036755, 1.4213700191513753967598177041509189,
    0.23397037352286015332580859649169724, 0.023855802790566132259852949459855154,
    0.0013554433403730756571570204016552452, 3.514112363930196343943593295655603e-05,
    2.5744412922565540418910127844209779e-07,
};
static const int numerator_cbrt_13_13_len = 14;

static const Real numerator_exp[] = {
    1, 0.5, 0.11999999999999999555910790149937384, 0.018333333333333333425851918718763045,
    0.0019927536231884057128793674706912498, 0.00016304347826086957803926047461118287,
    1.0351966873706003395257460586442733e-05, 5.1759834368530021211452039203715358e-07,
    2.04315135665250075022251601561743e-08, 6.3060227057175949591547168724629119e-10,
    1.4837700484041399485741487398290288e-11, 2.5291534915979658214332080792540263e-13,
    2.8101705462199623130636742246532468e-15, 1.5440497506703087532742710210031064e-17,
};
static const int numerator_exp_len = 14;

static const Real numerator_sinc[] = {
    1, 0, -0.15329450413674683417752930836286396, 0, 0.0061935496289174151243739352423745004,
    0, -0.00010152460010152459688945181248698191, 0, 7.8402368468051121814507481086464757e-07,
    0, -2.8843896017506438908679855732172109e-09, 0, 4.1601679535140991798873542454890572e-12,
};
static const int numerator_sinc_len = 13;

static const Real numerator_sqrt_13_13[] = {
    1, 6.75, 20.25, 35.578125, 40.60546875, 31.5615234375, 17.033203125, 6.387451171875,
    1.638885498046875, 0.27820587158203125, 0.02945709228515625, 0.0017573833465576171875,
    4.8816204071044921875e-05, 4.0233135223388671875e-07,
};
static const int numerator_sqrt_13_13_len = 14;

static const Real two_over_pi[] = {
    10680707.0, 7228996.0, 1387004.0, 2578385.0, 16069853.0, 12639074.0, 9804092.0,
    4427841.0, 16666979.0, 11263675.0, 12935607.0, 2387514.0, 4345298.0, 14681673.0,
    3074569.0, 13734428.0, 16653803.0, 1880361.0, 10960616.0, 8533493.0, 3062596.0,
    8710556.0, 7349940.0, 6258241.0, 3772886.0, 3769171.0, 3798172.0, 8675211.0,
    12450088.0, 3874808.0, 9961438.0, 366607.0, 15675153.0, 9132554.0, 7151469.0,
    3571407.0, 2607881.0, 12013382.0, 4155038.0, 6285869.0, 7677882.0, 13102053.0,
    15825725.0, 473591.0, 9065106.0, 15363067.0, 6271263.0, 9264392.0, 5636912.0,
    4652155.0, 7056368.0, 13614112.0, 10155062.0, 1944035.0, 9527646.0, 15080200.0,
    6658437.0, 6231200.0, 6832269.0, 16767104.0, 5075751.0, 3212806.0, 1398474.0,
    7579849.0, 6349435.0, 12618859.0, 4703257.0, 12806093.0, 14477321.0, 2786137.0,
    12875403.0, 9837734.0, 14528324.0, 13719321.0, 343717.0, 16713477.0, 4161075.0,
    15217346.0, 14569368.0, 3308987.0, 12795174.0, 15690526.0, 6224031.0, 3809077.0,
    13300351.0, 1935345.0, 2199676.0, 8135786.0, 16412373.0, 7810352.0, 4406037.0,
    12981429.0, 10295747.0, 12764333.0, 4279596.0, 6094860.0, 4619654.0, 2978275.0,
    10143387.0, 25139.0, 8180404.0, 9938868.0, 13980983.0, 16137943.0, 1577123.0,
    16545357.0, 2792804.0, 11261808.0, 16284771.0, 5746810.0, 15144215.0, 5654976.0,
    14276155.0, 3703975.0, 13312804.0, 7834326.0, 2315354.0, 12132096.0, 1772273.0,
    14667289.0, 16724383.0, 6954598.0, 6379417.0, 4717484.0, 14188414.0, 12018978.0,
    9037874.0, 6340582.0, 13485295.0, 603756.0, 13909853.0, 14147094.0, 14564184.0,
    9608158.0, 2630354.0, 15238696.0, 5069026.0, 3328710.0, 1499912.0, 13336032.0,
    5292055.0, 10952179.0, 6021144.0, 3412782.0, 6427267.0, 84099.0, 6000373.0,
    8368301.0, 15919390.0, 4409928.0,
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
    Real t = x - 1.0;
    return rational_eval(t, num, num_len, den, den_len);
}

int rem_pio2(Real x, Real *y) {
    if (mc_isnan(x) || mc_isinf(x)) {
        *y = NAN;
        return 0;
    }
    
    Real ax = mc_fabs(x);
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
    Real m = mc_frexp(ax, &e_int);
    int jk = e_int / 24;
    int k_start = jk - 2;
    if (k_start < 0) k_start = 0;
    
    Real accum_frac = 0.0;
    int accum_quad = 0;
    
    for (int k = k_start; k < k_start + 5; k++) {
        Real term = m * two_over_pi[k]; 
        int shift = e_int - 24 * (k + 1);
        Real val = mc_ldexp(term, shift);
        Real ipart;
        Real fpart = mc_modf(val, &ipart);
        long long i_ll = (long long)mc_fmod(ipart, 16.0); 
        accum_quad += i_ll;
        accum_frac += fpart;
    }
    
    Real f_ipart;
    Real f_fpart = mc_modf(accum_frac, &f_ipart);
    accum_quad += (long long)f_ipart;
    accum_frac = f_fpart;
    
    int n = accum_quad & 3;
    
    if (accum_frac >= 0.5) {
        n += 1;
        accum_frac -= 1.0;
    } else if (accum_frac < -0.5) {
        n -= 1;
        accum_frac += 1.0;
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
    Real m = mc_frexp(x, &e);
    if (e & 1) {
        m *= 2.0;
        e--;
    }
    Real y0 = _pade_about_1(m, numerator_sqrt_13_13, numerator_sqrt_13_13_len, 
                            denominator_sqrt_13_13, denominator_sqrt_13_13_len);
    y0 = 0.5 * (y0 + m / y0);
    y0 = 0.5 * (y0 + m / y0);
    
    return mc_ldexp(y0, e / 2);
}

Real mc_cbrt(Real x) {
    if (x == 0) return 0;
    int sign = 1;
    if (x < 0) { sign = -1; x = -x; }
    
    int e;
    Real m = mc_frexp(x, &e);
    int q = e / 3;
    int r = e - 3 * q;
    if (r < 0) { 
        q--;
        r += 3;
    }
    
    Real a = mc_ldexp(m, r);
    Real y0 = _pade_about_1(a, numerator_cbrt_13_13, numerator_cbrt_13_13_len,
                             denominator_cbrt_13_13, denominator_cbrt_13_13_len);
                             
    y0 = (2.0 * y0 + a / (y0 * y0)) / 3.0;
    y0 = (2.0 * y0 + a / (y0 * y0)) / 3.0;
    
    return sign * mc_ldexp(y0, q);
}

Real mc_exp(Real x) {
    Real k_float = mc_round(x / LN2_F_C);
    int k = (int)k_float;
    Real r = x - k * LN2_C; 
    
    Real e_r = rational_eval(r, numerator_exp, numerator_exp_len, denumerator_exp, denumerator_exp_len);
    return e_r * mc_ldexp(1.0, k);
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
    int n = rem_pio2(x, &y);
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
    int n = rem_pio2(x, &y);
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
    int n = rem_pio2(x, &y);
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
    int n = rem_pio2(x, &y);
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
    
    if (xf > 1.0) {
        return sign * (PI_2_C - mc_atan(1.0/xf));
    }
    
    int n = 0;
    while (xf > 0.3) {
        xf = xf / (1.0 + mc_sqrt(1.0 + xf*xf));
        n++;
    }
    Real a = _atan_series(xf);
    return sign * (a * mc_ldexp(1.0, n));
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
    if (mc_fabs(x) > 1) return NAN;
    if (x == 1) return PI_2_C;
    if (x == -1) return -PI_2_C;
    return mc_atan(x / mc_sqrt(1.0 - x*x));
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
    Real m = 1.0 + u;
    Real z = (m - 1.0) / (m + 1.0);
    return 2.0 * _atanh_series(z);
}

Real mc_ln(Real x) {
    if (x <= 0) return NAN;
    int e;
    Real m = mc_frexp(x, &e);
    Real z = (m - 1.0) / (m + 1.0);
    return 2.0 * _atanh_series(z) + e * LN2_C;
}

Real mc_sinh(Real x) {
    return 0.5 * (mc_exp(x) - mc_exp(-x));
}

Real mc_cosh(Real x) {
    return 0.5 * (mc_exp(x) + mc_exp(-x));
}

Real mc_sinc(Real x) {
    if (x == 0) return 1.0;
    return rational_eval(x, numerator_sinc, numerator_sinc_len, denumerator_sinc, denumerator_sinc_len);
}
