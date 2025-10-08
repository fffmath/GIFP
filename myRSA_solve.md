### myRSA_solve

myRSA writeup from [link](https://mp.weixin.qq.com/s/0sBfu94em2sR82OYDZF6zQ)

First, we notice that LFSR has a cycle. By calculation or random seed simulation, the cycle of LFSR is observed to be 2697. Given that we have obtained 4096 random bits, the subsequent 4096-2697 = 1399 bits will certainly be repetitious. This denotes that the initial 1399 bits of 'a' and the subsequent 1399 bits of 'b' overlap. Proceeding with this analysis, we can conclude that 

```python
bin(q2)[2:][256:256+877] == bin(p2)[2:][905:905+877]
# This leads us to reference the paper Generalized Implicit Factorization Problem.
```

The paper even provides the execution code fffmath/gifp.

We merely need to comment out some debugging information in the original code and then make minor adjustments to run it.

However, running it directly does not yield the answer. 

According to the guidance in readme, we should check the Groebner basis at this step.

Take note that the Groebner basis does not necessarily need to have four polynomials.

Noticing that the first two polynomials are comparably simple, we can directly ascertain the proportionality constants of x, y, and w. Since w itself is a prime number (as the paper suggests that w is p2), w corresponds to the denominator of the two equations, while x and y correspond to the numerator of the equations.

After solving the equation, we can directly decompose N1 and N2, and then follow the RSA process to obtain the flag.

```python
from Crypto.Util.number import long_to_bytes

def eliminate_N2(f, modular):
    tmp_poly = 0
    for mono in f.monomials():
        if f.monomial_coefficient(mono) % modular == 0:
            tmp_poly += mono * modular
        else:
            tmp_poly += mono * (f.monomial_coefficient(mono) % modular)
    return tmp_poly


## for simplicity, we exchange N1 and N2, and assume p > q
N1 = 15100254697650550107773880032815145863356657719287915742500114525591753087962467826081728465512892164117836132237310655696249972190691781679185814089899954980129273157108546566607320409558512492474972517904901612694329245705071789171594962844907667870548108438624866788136327638175865250706483350097727472981522495023856155253124778291684107340441685908190131143526592231859940556416271923298043631447630144435617140894108480182678930181019645093766210388896642127572162172851596331016756329494450522133805279328640942549500999876562756779916153474958590607156569686953857510763692124165713467629066731049974996526071
N2 = 11195108435418195710792588075406654238662413452040893604269481198631380853864777816171135346615239847585274781942826320773907414919521767450698159152141823148113043170072260905000812966959448737906045653134710039763987977024660093279241536270954380974093998238962759705207561900626656220185252467266349413165950122829268464816965028949121220409917771866453266522778041491886000765870296070557269360794230165147639201703312790976341766891628037850902489808393224528144341687117276366107884626925409318998153959791998809250576701129098030933612584038842347204032289231076557168670724255156010233010888918002630018693299
e = 65537
c = 4814924495615599863001719377787452659409530421473568305028025012566126400664362465878829645297418472853978736123334486954531551369698267539790007454131291197238666548347098462574649698959650399163371262093116196849478966536838813625531493756454672767925688817717023571267320336512019086040845336203876733170680765788230657626290346741730982737645243576591521512217549028162039336681342312618225110504010746912604698079363106352791499951364510694837846064947033912634697178711135807010770985698854383359436879061864935030256963597840531276583023488437671584864430423908673190679945703404235633491522955548722332086120

modulus_bit_length = 2048
alpha_bit_length = 256
share_bit_length = 877
beta1_bit_length = 10      
beta2_bit_length = 659

gamma_bit_length = share_bit_length
m = 6
R.<x, y, z, w> = PolynomialRing(ZZ, 4, order='lex')
f = x * z + 2 ** (beta2_bit_length + gamma_bit_length) * y * z + N2
X = Integer(2**beta2_bit_length)
Y = Integer(2**(modulus_bit_length - alpha_bit_length - gamma_bit_length - beta1_bit_length))
Z = Integer(2**alpha_bit_length)
W = Integer(2 ** (modulus_bit_length - alpha_bit_length))
M = Integer(2 ** (beta2_bit_length - beta1_bit_length))
alpha = alpha_bit_length / modulus_bit_length
t = round((1 - sqrt(alpha)) * m) 
s = round(sqrt(alpha) * m)

f = f.change_ring(ZZ)

qr = R.quotient(z*w - N2)

shifts = Sequence([], f.parent())

modular = M^m*N1^t

N2_inverse = inverse_mod(N2, modular)

for ii in range(m + 1):
    for jj in range(m - ii + 1):
        g = (y*z) ** jj * w ** s * f ** (ii) * M ** (m - ii) * N1 ** max(t - ii, 0) * N2_inverse ** min(ii + jj, s)
        g = qr(g).lift()
        g = eliminate_N2(g, modular)
        shifts.append(g)

L, monomials = shifts.coefficient_matrix()
monomials = vector(monomials)
bounds = [X, Y, Z, W]
factors = [monomial(*bounds) for monomial in monomials]
for i, factor in enumerate(factors):
    L.rescale_col(i, factor)

print(L.dimensions())
print(L.rank())
L = L.dense_matrix().LLL()
print('LLL done')
L = L.change_ring(QQ)
for i, factor in enumerate(factors):
    L.rescale_col(i, 1/factor)
H = Sequence([], f.parent().change_ring(QQ))
for h in filter(None, L*monomials):
    H.append(h)
H.insert(0, w * z - N2)
def check(denominator):
    if denominator == 1:
        return 0
    d = gcd(denominator, N2)
    if d == 1 or d == N2:
        return 0
    return d

G = ideal(H[:20]).groebner_basis()
print('groebner done')
for polynomial in G:
    for coefficient, monomial in polynomial:
        if check(coefficient.denominator()):
            p2 = check(coefficient.denominator())
            print(p2)
            break
    else:
        continue
    break

q2 = N2 // p2
assert p2 * q2 == N2
x0 = int(G[0](w=p2).univariate_polynomial().roots()[0][0])
y0 = int(G[1](w=p2).univariate_polynomial().roots()[0][0])

p1 = gcd(N1, p2+x0+y0*(2**int(gamma_bit_length+beta2_bit_length)))
q1 = N1 // p1

phi = (p1-1)*(q1-1)

d = inverse_mod(e, phi)

print(long_to_bytes(pow(c, d, N1)))
```
