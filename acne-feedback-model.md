# A nonlinear feedback model of acne persistence, remission and relapse

**Status:** AI-assisted theoretical research note, 2026-11-23. Unreviewed; no claim of a solved biological open problem, clinical efficacy, or verified novelty.

## Motivation
Lai et al. report a lipid-induced positive feedback mechanism involving Cutibacterium acnes lipase GehA and sebum (J Invest Dermatol, 2026; https://pubmed.ncbi.nlm.nih.gov/42526582/). Other research documents lipid-associated inflammatory macrophages (https://pmc.ncbi.nlm.nih.gov/articles/PMC9400695/), strain-associated porphyrin inflammation (https://pmc.ncbi.nlm.nih.gov/articles/PMC8188554/) and inhibition of lipase by oleic acid (https://pmc.ncbi.nlm.nih.gov/articles/PMC9675171/). These motivate, but do not validate, the following deliberately simplified hypothesis.

## Model
On a common normalized activity scale, x = sebum-associated activity, y = bacterial signaling, z = host inflammatory/lipogenic response. Assume positive rates a,b,c, K>0, q,g>=0:
```
x' = a(F(z)-x)
y' = b(x-y)
z' = c(y-z)
F(z) = q+g*z^2/(K^2+z^2)
```
All nonnegative trajectories are bounded and remain nonnegative. These are latent activity variables, not direct equivalences to clinical lesion counts, sebum mass or porphyrin fluorescence.

## Equilibria and thresholds
Equilibria satisfy x=y=z=s, with q=s-g*s^2/(K^2+s^2).
Writing t=s/K, G=g/K and Q=q/K, the equilibrium curve is Q=t-G*t^2/(1+t^2). Its derivative is 1-2G*t/(1+t^2)^2. The maximum of 2t/(1+t^2)^2 is 9/(8 sqrt(3)), at t=1/sqrt(3). Thus a range of nonnegative Q with three distinct equilibria exists iff G>8 sqrt(3)/9. Folds satisfy G=(1+t^2)^2/(2t), Q=t(1-t^2)/2; retain only Q>=0. At Q=0 a stable positive equilibrium coexists with the stable zero equilibrium when G>2; t=(G +/- sqrt(G^2-4))/2 gives the two nonzero equilibria. Equality is degenerate.

With A=a+b+c, B=ab+ac+bc, C=abc, the equilibrium characteristic polynomial is lambda^3+A lambda^2+B lambda+C(1-F'(s)). Routh-Hurwitz gives exponential stability iff 1-AB/C<F'(s)<1. For increasing F the two outer simple equilibria are stable and the middle has one positive eigenvalue.

## Global convergence (conditional)
For bounded, nonnegative C1 F with finitely many fixed points, assume inf F'>1-AB/C. Set h=z-F(z), H(z)=integral_0^z h(v)dv. Elimination yields z'''+Az''+Bz'+Ch(z)=0. Pick k>0 with (C/B)sup h'<k<A. The function
V=(z'')^2/2+kz'z''+(B+kA)(z')^2/2+C h(z)z'+kC H(z)
has derivative -(A-k)(z'')^2-(kB-Ch'(z))(z')^2. Bounded trajectories imply V is bounded below along each trajectory. LaSalle's principle forces the omega-limit set into the equilibria; finiteness and connectedness imply convergence to one equilibrium. For the increasing Hill response this condition holds automatically. No sustained oscillations in this model.

## Remission and control
In the three-simple-equilibrium case, the middle equilibrium is a saddle; its stable set separates the low and high attraction basins in the interior. By order comparison, if all states are strictly below the middle equilibrium, the trajectory tends to the low equilibrium. Replace g with (1-u)g, 0<=u<=1. With no intensity or side-effect penalty, maximal suppression minimizes each state and therefore is minimum-time optimal for entering the downward-closed low-state basin (the boundary itself is not remission). With r=min(a,b,c), M=max(q,x(0),y(0),z(0)), complete suppression yields the sufficient condition (M-q)exp(-rT)(1+rT+(rT)^2/2)<s_M-q. This is not a medication or duration recommendation.

## Extensions
A positive feedback model with follicular obstruction w' = d(z-w), bacterial potency sigma and resolution rho can take x'=a[q(h)+g phi(z)+k phi(w)-x], y'=b(sigma x-y), z'=c[y-(1+rho)z], phi(s)=s^2/(1+s^2). Its steady state reduces to s=Q+G phi(s), Q=sigma q(h)/(1+rho), G=sigma(g+k)/(1+rho). This is an equilibrium reduction, not proof of all transient dynamics.

For a hypothetical scalar spatial reduction s_t=D Laplacian(s)+q+g phi(s)-s, D>0, the sign of the planar bistable front speed (high state invading low) is the sign of integral from s_L to s_H of q+g phi(s)-s, equal to [qs+g(s-arctan s)-s^2/2] at the endpoints. Spatial coupling remains unverified.

## Calibration and identifiability
Three distinct positive equilibrium measurements on the common scale identify K^2=s1*s2+s1*s3+s2*s3, q=s1*s2*s3/K^2, g=s1+s2+s3-q, via the cubic s^3-(q+g)s^2+K^2*s-q*K^2=0. With continuous direct state measurements, a=x'/(F(z)-x), b=y'/(x-y), c=z'/(y-z) at informative points.

Alternatively, independently clamp z at 0,z1,z2 and measure equilibrium x. Set Ri=xi-q, with q=x at z=0. Then K^2=(R2-R1)/(R1/z1^2-R2/z2^2), g=R1*(1+K^2/z1^2), assuming nondegeneracy.

Eliminating hidden variables yields z'''+Az''+Bz'+Cz=Cq+Cg*z^2/(K^2+z^2). Multiplying by K^2+z^2 gives a linear regression with nine composite coefficients: z^2 z''', z''' (after separating the known first term), z^2 z'', z'', z^2 z', z', z^3, z, z^2, and constant. More precisely, move z^2 z''' to the left; the nine unknown coefficients multiply the remaining nine terms. Full column rank suffices to recover A,B,C,K,q,g, with rates unordered as roots of r^3-Ar^2+Br-C=0. Rate permutation and common activity-scale ambiguity limit identification from a single uncalibrated output; g/K is scale invariant.

Example: K=1,q=0.5,g=1 and K=0.6,q=1/18,g=1.284444... have the same s=1, F'(1)=0.5 but different equilibrium multiplicity. Thus equilibrium and local response alone cannot establish bistability.

## Proposed falsification and data requirements
Test increasing/decreasing controlled stimulation and washout on matched follicle cultures; measure sebum-related activity, direct GehA activity and inflammatory/lipogenic signals at repeated times, with independently calibrated scales and appropriate controls. Compare alternative feedback functions, account for negative feedback and biological variability, estimate uncertainties, and test predictions on held-out follicles. Raw human data identifying all six parameters have not been obtained. Porphyrins and lesion counts are imperfect proxies. No empirical parameter estimates, independent replication, or validated clinical predictions are claimed.

## Review request
Feedback on biological plausibility, existing models, mathematical assumptions, and access to suitable quantitative data is welcome. This is an exploratory model, not a clinical intervention proposal.
