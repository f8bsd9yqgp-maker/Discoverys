# Optimal quantum probes under Gaussian dephasing: complexity and elementary solvability

**Research note, 23 November 2026. AI-assisted candidate results; not peer reviewed.** Independent verification and priority checks are requested. No claim of experimental validation.

## Model and exact optimization

For N qubits, let G = sum_i w_i Z_i/2, g_x = w^T x/2 for x in {−1,+1}^N, and K_xy = exp[−(x−y)^T Sigma (x−y)/8] with Sigma positive semidefinite. The channel maps each density-matrix entry to K_xy exp[−i theta(g_x−g_y)] times that entry. We optimize local single-parameter QFI over input states (no adaptive controls or extra noiseless ancillas).

Take Gram vectors u_x for K. Pure probes suffice by convexity; input phases commute with the channel and do not affect QFI. Write p_x for their populations. The purification variational identity is F(p)/4 = min_H sum_x p_x ||(H−g_x I)u_x||^2, with H Hermitian (scalar shifts remove the mean). Thus F is concave in p. Sion minimax gives max_p F(p)/4 = min_H max_x ||(H−g_x I)u_x||^2, an exact SDP by Schur complements. Related optimization frameworks are already known.

For two differential sensors with w=(1,−1), Sigma=v[[1,r],[r,1]], set a²=exp(−v), b=exp(−v(1−r)), c=exp(−v(1+r)), d=(1+c)/2−a², e=(1−b)/2, k=a²(1−b)². Symmetry and concavity allow populations p/2 on 01,10 and (1−p)/2 on 00,11. Solving the SLD equation yields F(p)=4[b²p+kp(1−p)/(d(1−p)+ep)]. For v>0 and r<1, F''(p)=−8kde/[d(1−p)+ep]^3, and the Bell state is optimal iff b²>=2a²(1−b). If r=1 or v=0, the differential Bell state is noiseless.

## Candidate NP-hardness theorem

Reduce binary-encoded SUBSET SUM with positive integers a_i, target T>=0. Choose M>4+T+sum a_i, A_i=M^i+a_i, B_i=M^i, Q=T+sum M^i, i=1,...,m. Put q=(A_1,B_1,...,A_m,B_m,Q), Sigma=16qq^T, and G=Z_(2m+1)/2. Two configurations of opposite signal eigenvalue have equal environmental charge h=q^T x/2 iff sum_i(A_i u_i+B_i v_i)=Q for u_i,v_i in {−1,0,1}. Reduction modulo M first forces sum a_i u_i=T; dividing and comparing balanced base-M digits then forces u_i+v_i=1, so u_i,v_i are binary. Therefore a yes instance has F_max=1 by a protected two-configuration superposition.

For a no instance, opposite-signal charges are separated by at least one. Environmental Gram vectors may be represented by normalized Gaussian wave packets with position variance 1/64 centered at h (overlap exp[−8(h−h')²]). Multiplication by the signal eigenvalue of the nearest allowed charge is a Hermitian environmental estimator H. Its mean squared error for each u_x is at most P(|position−h|>=1/2)<=1/16 by Chebyshev. The purification bound gives F_max<=1/4. Thus a classical additive-1/4 approximation distinguishes yes/no (threshold 5/8). This is a binary-input weak NP-hardness reduction; it uses large, nonuniform couplings and does not exclude efficient algorithms for physically restricted instances.

## Candidate non-elementary obstruction (four qubits)

Let G=sum Z_i/2, Sigma=−ln(z) I, 0<z<1, K_xy=z^(Hamming(x,y)/2). Symmetry averaging and population concavity give a symmetric maximizer, with populations p/2 on weights 0,4; q/8 on weights 1,3; (1−p−q)/6 on weight 2. The accompanying Python certificate uses the permutation representation blocks (spin 2 and three spin 1 copies; spin 0 contributes no QFI) to derive the rational QFI, stationary polynomials and exact resultants. It is intended to verify an interior optimum at z=1/4 and factorization of the stationary resultant into p^20(p−1)^16 B(p,z)^4 H(p,z), up to parameter factors, where deg_p H=14 and

B=4(z−1)^4(9z+5)p²−4(z−1)^2(9z³−19z²−13z−1)p+(z+1)^3(9z²+2z+1).

At z=1/4, the H factor has modular factor-degree patterns (1,13), (1,1,2,5,5), (14) at primes 19,131,149 respectively. These imply S14 over Q: transitivity, a 13-cycle (primitivity), and a transposition (fifth power of the second cycle pattern). At z=1/100 and 1/4, the reported exact real-root counts are four and two. The proposed argument is that generic arithmetic Galois group is S14 by specialization; the geometric group is normal and cannot be trivial because a splitting into complex rational functions would give invariant real-root counts at nonsingular real specializations. It therefore contains A14, and the algebraic branch has nonsolvable monodromy. By topological Galois theory this excludes a finite elementary formula for a globally optimal probe on an open interval, even with complex constants, standard elementary functions and finite piecewise branching. A general selector would imply an elementary symmetric population by finite averaging of its populations. This does not exclude special functions, implicit roots, numerical methods, or arbitrary oracles.

**Audit request:** Please independently check the block reduction, resultant interpolation degree bounds, specialization argument (including scalar factors), interior optimality, and exact scope of the elementary-function theorem. The certificate is supplied for reproduction, not a substitute for expert review. Novelty relative to prior literature is unconfirmed.

## Related literature
- Ghosh et al., Optimal phase estimation in the presence of correlated dephasing, Quantum 10, 2251 (2026), https://doi.org/10.22331/q-2026-10-28-2251
- Fully-Optimized Quantum Metrology: Framework, Tools, and Applications, https://doi.org/10.1002/qute.202400094
- Hamann, Sekatski and Dür, Approximate decoherence free subspaces for distributed sensing, https://doi.org/10.1088/2058-9565/ac44de
- Khovanskii, One dimensional topological Galois theory, https://arxiv.org/abs/1904.03341
