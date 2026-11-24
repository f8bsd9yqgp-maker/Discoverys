# Certifiable HDR Imaging Under Adversarial Motion, Occlusion, and Sensor Noise

Date: 2026-11-23 (US Central time)
Status: Preliminary AI-assisted research; candidate contributions, not independently peer-reviewed or confirmed novel. Not a claim to solve unrestricted HDR deghosting.

## Abstract
We formulate fixed-exposure HDR recovery as a robust distinguishability problem, derive an exact exposure count in a clipped bounded-noise model, and extend the reconstruction to calibrated stochastic sensors and finite motion-trajectory graphs. The goal is certified radiance uncertainty, including abstention when observations do not constrain radiance sufficiently.

## 1. Candidate sharp bounded-error exposure count
Let x be log-radiance in [0,L], W>0 the clipping range, a_i fixed exposure offsets, and uncorrupted observations y_i=clip_[0,W](x+a_i+e_i), |e_i|<=delta. At most f measurements may be arbitrarily changed. Assume 0<=delta<=epsilon<L/2 and f a nonnegative integer. With continuous exposure offsets and unrestricted repetition, the proposed exact minimum is

N_min = (2f+1) ceil((L-2epsilon)/(W-2delta+2epsilon)).

Proof sketch: possible uncorrupted measurements are closed clipped intervals. For x<x', they are disjoint iff x'-x>2delta, x+a_i+delta<W, and x'+a_i-delta>0. A pair is distinguishable under f arbitrary replacements iff at least 2f+1 coordinates have disjoint intervals: otherwise split those coordinates between the two corruption patterns and choose shared measurements elsewhere. Uniform epsilon-accuracy is equivalent to distinguishability of every pair with separation greater than 2epsilon, using the midpoint of the feasible radiance extrema.

Write b_i=-a_i and H=W-2delta+2epsilon>0. For each x in [0,L-2epsilon), taking x' decreasing to x+2epsilon requires at least 2f+1 intervals [b_i+delta-2epsilon,b_i+W-delta) to cover x. This follows from finiteness of the exposure list; intervals with lower endpoints strictly above x cannot distinguish sufficiently close pairs. Conversely, every interval covering x distinguishes all x'>x+2epsilon. With m=ceil((L-2epsilon)/H), the points jH, j=0,...,m-1, belong to the target domain. No length-H half-open interval contains two of these points. Thus at least (2f+1)m exposures are necessary. Repeating each of m adjacent length-H intervals 2f+1 times attains the bound. This is a proof for a stylized model, not a physical camera guarantee. Independent mathematical review is requested.

## 2. Stochastic sensor extension
Model photons as P_i~Poisson(t_i R), read noise as independent Gaussian G_i with calibrated variance, and observations as quantized, clipped P_i+G_i (with gain/offset incorporated into calibration). Construct confidence sets I_i by inverting calibrated tests at level alpha/N. For a clean potential measurement vector, simultaneous coverage is at least 1-alpha by the union bound; the adversary may replace up to f entries even after observing this vector. Define S={R: R belongs to at least N-f observed confidence sets}. On simultaneous coverage, the true R belongs to S. If S is nonempty and bounded, the midpoint of its extrema has certified radius diam(S)/2 on the coverage event. Empty/unbounded sets require abstention. Calibration error, model misspecification and spatially simultaneous coverage need separate control.

For finite radiance, output and exposure alphabets and positive capture costs, exact nonadaptive minimax scheduling can be posed as a finite linear program per schedule: choose decoder probabilities d(a|z) summing to one; minimize rho subject to sum_y p_s(y|R)u_(R,y)<=rho for every R and u_(R,y)>=sum_(a:|a-R|>epsilon)d(a|z) for all z differing from y in at most f positions. Enumerate budget-feasible schedules to find a global optimum. This is a computational characterization, not a claim of efficient unrestricted optimization.

## 3. Joint motion and radiance certification
Use a finite layered directed acyclic graph with nodes representing candidate scene-point positions at exposure midpoints and edges representing admissible displacement. Require the reference node to be visited. A known speed envelope, discretization tolerance and radiance-variation envelope determine conservative interval expansion for each node. For example, under spatial Lipschitz constant K_s, temporal along-trajectory Lipschitz constant K_t, alignment error rho_i and speed bound v, the average-radiance discrepancy can be bounded by K_s(rho_i+v t_i/4)+K_t(|c_i-t_0|+t_i/4) under the corresponding image-formation assumptions.

For a candidate radiance R, let w_i(p,R)=0 if R is in the expanded interval at node p, and 1 otherwise. Forward dynamic programming computes F_i(p,R)=w_i(p,R)+min_(q->p) F_(i-1)(q,R). The backward recurrence B is symmetric. Enforce the reference anchor by restricting its layer; inaccessible states have infinite cost. Accept radiance when a complete path has cost <=f. An edge p->q is feasible iff F_i(p,R)+B_(i+1)(q,R)<=f. Maximizing the displacement of feasible edges over accepted radiances gives a sharp displacement bound within the finite graph model. It does not bound unobserved instantaneous velocity.

Node costs only change at interval endpoints; process all endpoints and intervening open cells separately. With V vertices and E edges this gives O(V(V+E)) time and O(V+E) working memory if output is streamed. Under simultaneous coverage across all V nodes (e.g., individual error alpha/V), correct motion envelopes and at most f corrupted frames, the true discretized trajectory and radiance remain feasible with probability >=1-alpha. A whole-image guarantee needs an additional simultaneous-coverage budget.

## 4. Identifiability and implementation limits
For finitely many exposure integrations, a smooth time-varying radiance perturbation can be nonzero at the reference instant while having zero integral over each exposure. Two nonnegative scenes can therefore have identical sensor-observation laws and different instantaneous reference radiances. For target separation >2epsilon, at least one must have failure probability >=1/2. This rules out universal guarantees without temporal assumptions; it does not solve practical unrestricted deghosting.

For full-resolution conservative fusion, replace each layer's candidate intervals with their enclosing interval, then use order statistics L_(N-f) and U_(f+1) of lower and upper endpoints. These enclose all radiances consistent with at least N-f intervals. They need not themselves describe the exact feasible set. Monotone calibrated endpoint functions permit neighborhood min/max before interval lookup.

The project explored Python/C++ CPU and FPGA prototypes. Earlier conversational benchmarks used simulated RAW inputs, not physical measurements. FPGA technology mapping and behavioral tests were explored, but routed implementation, timing closure, deployment, and sensor calibration on real hardware remain unverified. This note does not include independent benchmark reproduction or claim hardware performance.

## 5. Related research and novelty questions
- Hasinoff, Durand, Freeman, Noise-Optimal Capture for High Dynamic Range Photography (2010): https://people.csail.mit.edu/hasinoff/hdrnoise/
- Kalantari and Ramamoorthi, Deep High Dynamic Range Imaging of Dynamic Scenes (2017): https://cseweb.ucsd.edu/~viscomp/projects/SIG17HDR/
- Debevec and Malik, Recovering High Dynamic Range Radiance Maps from Photographs (1997): https://www.pauldebevec.com/Research/HDR/
- HDR Imaging with Quanta Image Sensors: Theoretical Limits and Optimal Reconstruction: https://arxiv.org/abs/2011.03614

Potential contributions for review are the precise robust exposure-count formula and joint finite-trajectory radiance/displacement certification. Existing interval covering, robust estimation, statistical testing and dynamic programming methods are not themselves new. No exhaustive prior-art search has established priority.

## 6. Remaining work
Independent proof review; comparison against prior art; efficient statistical schedule optimization beyond enumeration; physical motion-envelope validation and calibrated RAW measurements; comparison against real HDR deghosting benchmarks; device-specific FPGA synthesis, routing and setup/hold signoff; end-to-end hardware experiments. This project is a candidate research contribution, not a validated discovery or a commercially ready camera.
