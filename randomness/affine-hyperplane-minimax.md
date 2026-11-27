# Exact minimax separation between linear and nonlinear full-entropy extraction from binary affine hyperplanes

2026-11-27. AI-assisted research note. Candidate result; independent novelty and review outstanding.

## Theorem

For every n>=3, let the source family consist of uniform laws on H(a,c)={x in F_2^n : a dot x=c}, where a is nonzero and c in F_2. Put M=2^(n-1). For deterministic f:F_2^n -> F_2^(n-1), let risk(f)=max_(a,c) TV(f(U_H),U_(n-1)). The minimax risk is exactly 1/4; restricted to linear (or affine) maps it is exactly 1/2. Both are attained.

## Construction and proof of upper bound

Write x=(b,z) with b a bit and z in F_2^(n-1). Choose A invertible with I+A invertible. Such A exists by representing multiplication by any field element distinct from 0 and 1 in F_(2^(n-1)). Set f(0,z)=z and f(1,z)=Az. The two preimages of output y are (0,y) and (1,A^(-1)y). For a=(alpha,u) with u nonzero, the two equations defining membership in H are u dot y=c and u dot A^(-1)y=c+alpha. The two nonzero linear forms are distinct: equality would imply u^T(A-I)=0, impossible. Hence the forms are independent over F_2, and exactly M/4 outputs have neither preimage in H, M/2 have one, and M/4 have both. If u=0 then alpha=1, and each output has exactly one preimage in H. Since a uniform H has M points, TV is the number of missing outputs divided by M (the output probabilities are integer multiples of 1/M). Thus risk(f)=1/4.

## Universal lower bound (including arbitrary unbalanced fibers)

Consider any deterministic f, with output fibers A_y and sizes k_y, summing to 2M. Pick a uniformly from the 2M-1 nonzero vectors, and c uniformly from two bits. A singleton is missed with probability 1/2. A pair of distinct points is missed with probability q=(M-1)/(2(2M-1)): their difference is orthogonal to exactly M-1 nonzero vectors, and then one choice of c misses both. An empty fiber is missed with probability 1.

For a fiber of size k, its miss probability is at least q+(k-2)(q-1/2). For k=0 the right side is 1-q<=1; for k=1,2 equality holds. For k>=3 it is at most 2q-1/2<0, so the inequality is automatic. Summing this inequality over M fibers cancels the size-dependent term and yields expected missing count >= Mq = M/4 - M/(4(2M-1)). Since M is divisible by 4 and the subtracted term lies strictly between zero and one, some H misses at least M/4 outputs. Thus every map has risk >=1/4.

## Linear lower bound and attainment

If a linear map has rank at most n-2, its image has size <=M/2 and its risk is at least 1/2. For a rank n-1 map its kernel is a line. Choose a nonzero a annihilating that line; H is then a union of whole kernel cosets and its image has exactly M/2 outputs, giving TV=1/2. Conversely each hyperplane under a rank n-1 map has image of size M or M/2, and the image law is uniform on that image. Thus every full-rank linear map has risk exactly 1/2. Translating the output does not change TV. QED.

## Reproducibility

The following Python 3/NumPy test was executed on 2026-11-27. All 4,064 hyperplanes for n=3 through 10 were tested, and all 65,536 deterministic maps for n=3 were enumerated. All checks passed. The general theorem rests on the proof, not the tests.

```python
import numpy as np
for n in range(3,11):
    m=n-1; M=1<<m; N=2*M
    # Block diagonal A using multiplication by x modulo
    # x^2+x+1 and, for odd m, x^3+x+1.
    def A(z):
        out=0; pos=0
        if m%2:
            a=(z>>pos)&7
            out|=(((a<<1)&7)^(3 if a&4 else 0))<<pos
            pos+=3
        while pos<m:
            a=(z>>pos)&3
            out|=(((a<<1)&3)^(3 if a&2 else 0))<<pos
            pos+=2
        return out
    f=np.array(list(range(M))+[A(z) for z in range(M)])
    assert len(set(f[M:]))==M
    assert len({z^A(z) for z in range(M)})==M
    worst=0
    for a in range(1,N):
        for c in range(2):
            mask=np.array([(a&x).bit_count()%2==c for x in range(N)])
            miss=(np.bincount(f[mask],minlength=M)==0).sum()
            worst=max(worst,miss)
    assert worst==M//4
    print(n,2*(N-1),worst/M)
N=8; M=4
f=(np.arange(M**N)[:,None]//M**np.arange(N))%M
worst=np.zeros(len(f),dtype=int)
for a in range(1,N):
    for c in range(2):
        mask=[x for x in range(N) if (a&x).bit_count()%2==c]
        missing=sum(~(f[:,mask]==j).any(axis=1) for j in range(M))
        worst=np.maximum(worst,missing)
assert worst.min()==1
print('all maps',len(f),'min',worst.min()/M)
```

## Context and limitations

The result concerns the specific codimension-one affine-source family and output length equal to the input source entropy. A constant TV error of 1/4 does not imply cryptographic suitability or security with side information. It is a mathematical optimality result, not a claim of physically extracting perfectly uniform bits. No generalization to higher codimension is claimed. Prior art and originality remain unconfirmed; related background includes:
- Chattopadhyay, Goodman and Liao, Affine Extractors for Almost Logarithmic Entropy (FOCS 2021), https://eccc.weizmann.ac.il/report/2021/075/
- Rao, Extractors for Low-Weight Affine Sources, https://eccc.weizmann.ac.il/eccc-reports/2008/TR08-015/index.html
- A New Approach to Affine Extractors and Dispersers, https://eccc.weizmann.ac.il/report/2010/064/
- Chattopadhyay, Recent Advances in Randomness Extraction, Entropy 2022, https://www.mdpi.com/1099-4300/24/7/880

Feedback on correctness and known antecedents is invited.