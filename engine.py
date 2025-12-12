#README.md. Addition, substraction, multiplication, division and taking modulus are allowed only (with raising to integer and simple rational powers too). No/Some direct hardware communication is allowed :((. Build some built in functions manually.

### Pade approximants [13/13]
numerator_LN1P=[
    1,
    6,
    4741/300,
    1441/60,
    107091/4600,
    8638/575,
    263111/40250,
    153081/80500,
    395243/1101240,
    28549/688275,
    605453/228813200,
    785633/10296594000,
    1145993/1873980108000,
    0
]

denominator_LN1P = [
    1,
    13/2,
    468/25,
    1573/50,
    1573/46,
    11583/4600,
    10296/805,
    2574/575,
    11583/10925,
    143/874,
    572/37145,
    117/148580,
    13/742900,
    1/10400600
]

numerator_arcsin = [
    1,
    0,
    -48612811065193217720717677178392447 / 17521502337917614029717740270538300,
    0,
    1542880498363289248597622552896377147 / 537326071696140163578010701629841200,
    0,
    -25629388878535757822653318119825725571 / 18806412509364905725230374557044442000,
    0,
    59997700449506246047746748284324063301 / 205817378502489528256921219152294373248,
    0,
    -3929590045745452902082508393999567898951 / 171057110133180185706863413251017990210560,
    0,
    17667852791636166403806066845986733616137 / 54939518901597871409380837432385778032332800,
    0
]

denominator_arcsin = [
    1,
    0,
    -17177687151615384464112433518938499 / 5840500779305871343239246756846100,
    0,
    70638889839100585257826495902581175 / 21493042867845606543120428065193648,
    0,
    -260976433493184241308750702390067677 / 150451300074919245801842996456355536,
    0,
    7106212846910135271370912693952576463 / 16334712579562660972771525329547172480,
    0,
    -100234105411492207380201907281393735693 / 2221520910820521892296927444818415457280,
    0,
    55355990129067361992384874772991295329 / 44430418216410437845938548896368309145600,
    1
]


### 2^n coefficients not used anymore
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


### Pade for 2^x [13/13] not used anymore
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
    1,
    1/2,
    3/25,
    11/600,
    11/5520,
    3/18400,
    1/96600,
    1/1932000,
    1/48944000,
    1/1585785600,
    1/67395888000,
    1/3953892096000,
    1/355850288640000,
    1/64764752532480000
]

denumerator_exp = [
    1,
   -1/2,
    3/25,
   -11/600,
    11/5520,
   -3/18400,
    1/96600,
   -1/1932000,
    1/48944000,
   -1/1585785600,
    1/67395888000,
   -1/3953892096000,
    1/355850288640000,
   -1/64764752532480000
]

### Some other constants. ###
LN2 = 0.693147180559945309417232121458176568
PI=3.1415926535897932384626433832795028841971693993751058
PI_2=PI/2
PI_4=PI/4



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
def poly_eval(x, coeffs):
    """Evaluate sum_{k=0}^n coeffs[k] * x^k using Horner's method."""
    y = 0.0
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

###def fpow2(n): ###Need to code how to compute 2^k fast :( . Still stuck trying to find my way out against bit shifting. Thinking of taking square root of 2^k for range reduction.
 #   global pow2_positive
 #   global pow2_negative
  #  global numerator_pow2
  #  global denumerator_pow2
  #  pass
        
    

def EXP(x):
    global numerator_exp
    global denumerator_exp
    k=round(x/LN2)
    r=x-k*LN2
    e_power_r= rational_eval(r,numerator_exp, denumerator_exp)
    pow2_k = 2**k ### need to replace with bit shifting thingy or mathematical equivalent. Building from scratch##
    return e_power_r * pow2_k

def sinc(x):
    global numerator_sinc
    global denumerator_sinc #???
    if x==0:
        return 1
    else:
        return rational_eval(x, numerator_sinc, denumerator_sinc)

def sin(x):
    return sinc(x)*x #???

def cos(x):
    return sin(PI_2 - x)

def tan(x):
    return ((sin(x))/(cos(x)))

def cot(x):
    return 1/tan(x)

def asin(x): #replace methods with newton-raphson
    if x>0.9:
        pass
    if x< -0.9:
        pass    
    return rational_eval(x, numerator_arcsin, denominator_arcsin)

def acos(x):
    return PI_2 - asin(x)

def sinh(x):
    return (0.5*(EXP(x) - EXP(-x)))

def cosh(x):
    return (0.5*(EXP(x) + EXP(-x)))

def atan(x):
    return asin(x/(1+x*x))

def ln(x):
    if x<=0:
        raise ValueError("Can't evaluate in complex set")
    #### need to finish later...


###Manual freexp inferior and crude.

def frexp_manual(x):
    if x <= 0:
        raise ValueError("x must be positive")
    n = 0
    m = x
    while m >= 2.0:
        m /= 2.0
        n += 1
    while m < 1.0:
        m *= 2.0
        n -= 1
    return m, n   

