# Universal minimax optimality of antipodal randomness extraction under heterogeneous independent biases

2026-11-27. AI-assisted research note for independent review. Novelty is unconfirmed.

## Theorem

Let n>=2 and let independent binary inputs X_i have absolute biases |E[(-1)^X_i]|<=rho_i, with 0<=rho_i<=1. For every deterministic f:{0,1}^n -> {0,1}^{n-1}, let its risk be the supremum over these input laws of TV(Law(f(X)), U_{n-1}). The minimax risk is attained by F(x)=(x_1 XOR x_n,...,x_{n-1} XOR x_n).

Put p(x)=2^{-n} product_i(1+rho_i(-1)^x_i), and let xbar denote the bitwise complement. The exact minimax value is

epsilon* = (1/4) sum_x |p(x)+p(xbar)-2^{1-n}|.

## Proof

Let N=2^n, M=N/2, c=1/M. By separate convexity in each bias coordinate (the output law is affine in each coordinate), the maximum risk is attained at a bias-box vertex. Write its distribution p_s(x)=p(x XOR s), s in {0,1}^n. Let H(A)=E_s min(sum_{x in A}p_s(x), c) for uniform s. The overlap 1-TV for an extractor is the sum of H over its fibers when averaged over s.

All singleton H values equal h1. Set h2=H({0,1^n}). For every distinct x,y, the variables p_s(x),p_s(y) have identical marginal distributions as s varies. Among all couplings of a finite equiprobable list with itself, the expectation of min(a+b,c) is maximized by oppositely ordered pairing. Indeed, for a<=a' and b<=b', concavity of phi(t)=min(t,c) implies phi(a+b)+phi(a'+b') <= phi(a+b')+phi(a'+b); successive exchanges sort the pairing in opposite order. For all rho_i<1, p(z)p(zbar)=4^{-n} product_i(1-rho_i^2)>0, so complement pairs achieve that opposite ordering (ties cause no difficulty). Boundary cases rho_i=1 follow by continuity. Hence every two-element set A has H(A)<=h2.

H is monotone submodular, because each summand is a nondecreasing concave function of a nonnegative additive set weight. For a set A of size k>=2, choose distinct a,b in A. Starting with {a}, each subsequent element z contributes at most H({a,z})-H({a}) <= h2-h1. Thus H(A)<=h1+(k-1)(h2-h1)=h2+(k-2)(h2-h1). This also holds for k=1 with equality, and for k=0 because submodularity gives h2<=2h1. Summing over M fibers A_j of any deterministic extractor, with sum_j |A_j|=N=2M, yields sum_j H(A_j)<=M h2.

The antipodal partition attains equality in this average bound. Furthermore, translation by s permutes its fibers, so the output TV is identical for all vertex laws. Its worst-case risk therefore equals its average risk 1-M h2, and every competitor has worst-case risk at least its average risk >=1-M h2. The formula follows by counting each antipodal fiber twice. QED.

## Computation

At rho_i=0.1 for all i, the raw TV and extracted TV are respectively:
n=4: 0.07848125, 0.0075125
n=8: 0.11367506113281253, 0.019523657890624938
n=16: 0.16108964736627623, 0.03610817503688844

Independent exhaustive enumeration in this submission audit verified the minimum for all 16 maps for n=2 and all 65,536 maps for n=3, for four bias vectors each. Reproduction (Python 3 + NumPy):

```python
from itertools import product
import numpy as np
for n in [2,3]:
    N=2**n; M=N//2
    X=np.array(list(product([-1,1],repeat=n)))
    f=(np.arange(M**N)[:,None]//(M**np.arange(N)))%M
    for rho in [np.full(n,.1),np.linspace(.1,.8,n),
                np.linspace(0,1,n),np.full(n,.9)]:
        p=np.prod((1+X*rho)/2,axis=1)
        Q=p[np.bitwise_xor(np.arange(N)[:,None],np.arange(N)[None,:])]
        risks=np.zeros(len(f))
        for q in Q:
            v=np.stack([(f==j)@q for j in range(M)],axis=1)
            risks=np.maximum(risks,.5*np.abs(v-1/M).sum(axis=1))
        expected=.25*np.abs(p+p[::-1]-1/M).sum()
        assert abs(risks.min()-expected)<1e-10
    print("passed",n,len(f))
```

## Context and limitations

The XOR construction is a standard linear corrector (its row space is the even-weight code). The candidate contribution is an exact minimax total-variation proof over all deterministic maps, including unbalanced and nonlinear partitions, for this particular independent heterogeneous bias-box model and output length. No originality claim is established. The proof uses classical rearrangement and submodularity arguments. This does not establish low TV for arbitrary n or rho, security with side information, correlated-source performance, or hardware suitability.

Related: Skorski, Soto, and Günlü, Exact Bias of Linear TRNG Correctors — Spectral Approach, arXiv:2509.26393, https://arxiv.org/abs/2509.26393 (linear corrector analysis and references to earlier work); Dichtl, Bad and Good Ways of Post-Processing Biased Physical Random Numbers (2007); Grujić, Optimizing Linear Correctors: A Tight Output Min-Entropy Bound and Selection Technique (2024). Independent checks for prior art and proof correctness are invited.
