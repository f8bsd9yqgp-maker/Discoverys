"""Exact symbolic certificate for non-elementary four-qubit probe optimization.
Requires Python 3 and SymPy. No numerical optimization is needed.
Run: python quantum_sensor_monodromy_certificate.py
"""
import sympy as s
from math import comb
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ

p,q,a,z=s.symbols('p q a z')
R=s.Rational

def block(n,probs,mul):
    K=s.Matrix(n+1,n+1,lambda i,j:sum(comb(i,k)*comb(n-i,j-k)*a**(i+j-2*k) for k in range(max(0,i+j-n),min(i,j)+1)))
    rho=mul*K*s.diag(*probs)
    G=s.diag(*[R(n,2)-i for i in range(n+1)])
    U=s.Matrix.hstack(*[s.eye(n+1)[:,i]+s.eye(n+1)[:,n-i] if i!=n-i else s.eye(n+1)[:,i] for i in range(n//2+1)],*[s.eye(n+1)[:,i]-s.eye(n+1)[:,n-i] for i in range(n//2)])
    rho=U.inv()*rho*U; G=U.inv()*G*U; d=n//2+1
    A=rho[:d,:d]; B=rho[d:,d:]; C=G[:d,d:]; E=G[d:,:d]
    D=A*C-C*B; DD=B*E-E*A
    M=s.kronecker_product(A,s.eye(n//2))+s.kronecker_product(s.eye(d),B.T)
    dm=DomainMatrix.from_Matrix(M); dv=DomainMatrix.from_Matrix(s.Matrix(list(D)))
    dm,dv=dm.unify(dv); v,den=dm.solve_den(dv); L=v.to_Matrix().reshape(d,n//2)
    return s.factor(s.cancel(-4*s.trace(DD*L)/dm.domain.to_sympy(den)))

probs=[p/2,q/8,(1-p-q)/6,q/8,p/2]
F=s.factor(block(4,probs,1)+3*block(2,probs[1:4],1-a*a)).subs(a**2,z)
n,d=s.fraction(F); n=s.expand(n); d=s.expand(d)
def stationarity(v):
    factors=s.factor_list(s.diff(n,v)*d-n*s.diff(d,v))[1]
    nonconstant=[(x,e) for x,e in factors if x.has(p,q)]
    assert len(nonconstant)==1 and nonconstant[0][1]==1
    return nonconstant[0][0]
P,Q=stationarity(p),stationarity(q)
assert all(s.degree(x,q)==10 and s.degree(x,z)==7 for x in (P,Q))

B=4*(z-1)**4*(9*z+5)*p**2-4*(z-1)**2*(9*z**3-19*z**2-13*z-1)*p+(z+1)**3*(9*z*z+2*z+1)
for v,real_count in [(R(1,4),2),(R(1,100),4)]:
    r=s.Poly(s.resultant(P.subs(z,v),Q.subs(z,v),q),p)
    assert r.degree()==58
    h,rem=r.div(s.Poly(p**20*(p-1)**16*B.subs(z,v)**4,p)); assert rem.is_zero
    assert h.degree()==14
    h=h.clear_denoms()[1].primitive()[1]
    assert s.gcd(h,h.diff()).degree()==0
    assert h.count_roots(-s.oo,s.oo)==real_count
    if v==R(1,4):
        for prime,target in [(19,[1,13]),(131,[1,1,2,5,5]),(149,[14])]:
            fs=s.factor_list(h,modulus=prime)[1]
            assert h.LC()%prime!=0 and all(e==1 for _,e in fs)
            assert sorted(f.degree() for f,e in fs)==target
    print('Specialization',v,'passed; real roots:',real_count,flush=True)

f=s.cancel(F.subs(z,R(1,4)))
assert f.subs({p:R(1,8),q:R(1,2)})==1
for ff,x,t in [(s.cancel(f.subs(p,0)),q,R(13,20)),(s.cancel(f.subs(q,0)),p,R(9,20)),(s.cancel(f.subs(q,1-p)),p,R(29,100))]:
    df=s.diff(ff,x).subs(x,t); val=ff.subs(x,t)
    assert max(val-t*df,val+(1-t)*df)<1
assert all(c>0 for c in s.Poly(B.subs(z,R(1,4)),p).all_coeffs())
print('Interior optimality at z=1/4 certified',flush=True)

# R(z,p)=Res_q(P,Q) has z-degree <= 140 (20x20 Sylvester determinant).
# Degree in p <= 58 and divisibility by the monic p^20(p-1)^16 follow
# from 441 evaluations. The quotient U has p-degree <=22 and z-degree <=140.
# The pseudo-remainder prem_p(U,B^4) has z-degree <=140+15*20=440.
# Therefore 441 distinct zero evaluations prove the full polynomial identity.
ring2,x,y=ring('q,p',ZZ); ring1,x1=ring('p',ZZ)
termsP=s.Poly(P,z,p,q).terms(); termsQ=s.Poly(Q,z,p,q).terms()
for v in range(2,443):
    powers=[v**i for i in range(8)]
    def ev(terms):
        dd={}
        for (i,j,k),c in terms: dd[k,j]=dd.get((k,j),0)+int(c)*powers[i]
        return ring2.from_dict(dd)
    pp,qq=ev(termsP),ev(termsQ)
    assert pp.degree(x)==10 and qq.degree(x)==10
    r=ring1.from_dict(dict(pp.resultant(qq)))
    assert r.degree()<=58
    u,rem=divmod(r,x1**20*(x1-1)**16); assert not rem
    b=4*(v-1)**4*(9*v+5)*x1**2-4*(v-1)**2*(9*v**3-19*v**2-13*v-1)*x1+(v+1)**3*(9*v*v+2*v+1)
    assert not u.prem(b**4)
    if v%50==0: print('Checked interpolation point',v,flush=True)
print('All 441 exact resultant identity checks passed.',flush=True)
