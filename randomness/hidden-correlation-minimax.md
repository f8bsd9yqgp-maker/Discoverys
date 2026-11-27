# An exact minimax transition in deterministic randomness extraction under hidden higher-order dependence

Date: 2026-11-27. AI-assisted research note for independent review. Novelty and practical applicability are unverified; this is not a claim of an independently established breakthrough.

## Model and theorem

Let n >= 2, 0 < rho <= 1/2, delta >= 0. Let P be a law on {-1,+1}^n with marginal means mu_i satisfying |mu_i| <= rho, every proper subset mutually independent, and TV(P,Q) <= delta, where Q is the product of its marginals. For deterministic f:{-1,+1}^n -> {-1,+1}, define B(f)=sup_P |E_P f|. Then

B* = min_f B(f) = min(rho^n + 2 delta, rho^(n-1)).

The threshold is delta_c = rho^(n-1)(1-rho)/2. Full parity is optimal at or below the threshold; parity on any n-1 coordinates is optimal at or above it. Output TV from a fair bit is B/2.

## Proof

Write chi_S(x)=product of x_i over i in S. The proper-subset moments of P and Q agree, including the empty moment. Orthogonality of Walsh characters yields P(x)-Q(x)=2^(-n) t chi_[n](x), with t=E_P chi_[n]-E_Q chi_[n]. Consequently TV(P,Q)=|t|/2.

Full parity has bias at most rho^n+2 delta. For delta <= delta_c, take mu_i=rho and t=2 delta. This defines a probability law: the product law satisfies 2^n Q(x) >= (1-rho)^n >= rho^(n-1)(1-rho) >= t, and the added term sums to zero. Its proper moments are unchanged. This attains the full-parity bound. For delta >= delta_c, the admissible law at delta_c forces full-parity risk at least rho^(n-1). Any parity on n-1 coordinates has risk exactly rho^(n-1), since those coordinates are independent and the bound is attained by independent inputs.

It remains to lower-bound every non-full-parity Boolean function. If f changes sign across every cube edge, connectedness implies f equals full parity or its negative. Otherwise, some adjacent pair has equal outputs. Average f over that coordinate (with uniform law), giving g on m=n-1 coordinates with ||g||_infinity=1.

For y in {-1,+1}^m, define (T_rho g)(y)=E[g(Z)] where the Z_i are independent and E Z_i=rho y_i. In one coordinate this operator is the matrix A=(1/2)[[1+rho,1-rho],[1-rho,1+rho]]. Its inverse has absolute row sum 1/rho. The m-coordinate operator is A tensor-powered m times, so its inverse has induced infinity norm rho^(-m). Thus 1=||g||_infinity <= rho^(-m)||T_rho g||_infinity. Each expectation in T_rho g is realized by an admissible product law with the averaged coordinate unbiased. Hence B(f)>=rho^(n-1). Combining the bounds proves the theorem.

## Interpretation and limitations

The theorem is for unconditional classical output bias and a highly structured source class with exactly independent proper subsets. It does not imply security against side information, entropy certification, independence over repeated blocks, or superiority for actual hardware. Empirical applicability and independent originality require further study. The proof does not claim the formula for rho>1/2.

## Small-instance reproducibility

The following Python 3 / NumPy code enumerates all Boolean functions for n=2,3,4 and tests the exact minimax formula using a finite family of admissible adversarial distributions. The analytical proof, not the enumeration, establishes the general result.

```python
import numpy as np
from itertools import product
for n in (2,3,4):
    X=np.array(list(product([-1,1],repeat=n)))
    N=len(X); parity=X.prod(axis=1)
    fs=2*((np.arange(2**N,dtype=np.uint32)[:,None]
          >> np.arange(N,dtype=np.uint32))&1).astype(np.int8)-1
    for rho in (0.15,0.3,0.5):
        dc=rho**(n-1)*(1-rho)/2
        for delta in (0,dc/2,dc,2*dc):
            means=np.array(list(product([-rho,0,rho],repeat=n)))
            Q=np.prod((1+means[:,None,:]*X[None,:,:])/2,axis=2)
            q=np.prod((1+rho*X)/2,axis=1)
            t=2*min(delta,dc)
            p=q+t*parity/N
            assert p.min()>=-1e-12 and abs(p.sum()-1)<1e-12
            Q=np.vstack((Q,p))
            risk=np.max(np.abs(fs@Q.T),axis=1)
            expected=min(rho**n+2*delta,rho**(n-1))
            assert abs(risk.min()-expected)<1e-10
```

## Related background

- Eshan Chattopadhyay, Recent Advances in Randomness Extraction (2022), https://pmc.ncbi.nlm.nih.gov/articles/PMC9321789/
- COSIC, Deterministic extraction for truly random bits, https://www.esat.kuleuven.be/cosic/blog/co6gc-deterministic-extraction-for-truly-random-bits/
- Ryan O'Donnell, Analysis of Boolean Functions (2014), Boolean Fourier analysis and noise operators.

These are background references, not evidence of priority. Feedback on correctness, prior art and useful generalizations is invited.
