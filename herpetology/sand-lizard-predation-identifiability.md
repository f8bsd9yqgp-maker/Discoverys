# Sand lizard predation inference: structural identifiability and reanalysis

**Status:** Preliminary AI-assisted technical note, 2026-11-27. Not peer reviewed; novelty and independent cause-specific predictive validation not established.

## Source
Smolinský et al. (2022), High predation risk decimates survival during the reproduction season, Ecology and Evolution, https://doi.org/10.1002/ece3.9407 . Data: https://doi.org/10.5061/dryad.q83bk3jm9 .

## Identifiability
Write S_i=1+E_i-P_i, i in {M,F}, and T=E_M/E_F>0. Equation (8) of the paper depends on parameters only through S_M, S_F, T, and S_M/S_F. At fixed positive S_M,S_F<=1 and T, the entire physically admissible fiber (0<=E_i<=P_i<=1) is E_F=x, E_M=Tx, P_F=1-S_F+x, P_M=1-S_M+Tx, with 0<x<=min(S_F,S_M/T). Thus the likelihood is constant along this fiber wherever the original parameter bounds permit. The posterior can still be proper and informative under priors; convergence does not establish likelihood identifiability. Exact matching of all four target summaries requires the target sex ratio to equal S_M/S_F; noisy observations need not match exactly. This is a statement about the specified likelihood, not a refutation of all biological conclusions.

For a closed cohort with initial and final male fractions q0,q1 and sex-specific survival, S_M/S_F=[q1/(1-q1)]/[q0/(1-q0)]. Cross-sectional age groups are not necessarily a closed cohort, and apparent survival can include emigration.

## Data reanalysis
Dryad archive version 20220818: 507 rows, 344 individual IDs, 180 predation.attempt flags. The README defines this flag as first recorded regeneration or reduced tail length at recapture, not a witnessed attack. Dropping missing tailL leaves 506 records and 343 individuals. Taking the last capture for each ID gives adults M=112,F=135; young M=48,F=48; regenerated adult tails M=68,F=79. Ratio of male fractions adult/young=0.906883; ratio of sex odds=0.829630. Across consecutive recapture intervals, flags occur in 13/78 male and 14/85 female intervals; intervals are not independent and vary in duration. Exploratory first-recapture-only counts: 9/48 male and 8/57 female, Fisher exact p=0.599. These do not identify actual attack probabilities.

Reproduce the counts (Python/pandas; supply path to Dryad ZIP):
```python
import pandas as pd, zipfile, sys
with zipfile.ZipFile(sys.argv[1]) as z:
    d = pd.read_csv(z.open("Lacerta-agilis.txt"), sep="\t")
print(len(d), d.animal.nunique(), d["predation.attempt"].sum())
d=d.dropna(subset=["tailL"]).copy()
d["date"]=pd.to_datetime(d[["year","month","day"]])
d=d.sort_values(["animal","date"])
print(d.groupby("animal").tail(1).groupby(["sex","age"]).regenerate.agg(["count","sum"]))
print(d[d.groupby("animal").cumcount()>0].groupby("sex")["predation.attempt"].agg(["count","sum"]))
```

## Competing risks and validation
For a fixed horizon, define S as true survival, M as probability of non-predatory death and F_p as fraction of survivors with new predator-caused tail loss. If deaths have exclusive, exhaustive cause labels, probability of fatal predation or surviving with predator-caused tail loss is 1-S-M+S F_p. This is an event decomposition, not by itself a predictive model or a novel biological law. Attacks without either outcome remain unobserved. Apparent survival must not be substituted for true survival without addressing emigration. If N monitored animals include D confirmed predator deaths and U unresolved fates, fatal predation fraction is bounded by D/N and (D+U)/N, assuming the remaining outcomes are verified non-predator-death outcomes within the horizon.

External candidate data: Wieczorek et al. (2020), https://doi.org/10.1163/15685381-bja10018 ; Simbula et al. (2025), Lizards unplugged, Herpetological Journal 35(4). Need individual-level dated attack observations, new tail-loss causes, confirmed mortality causes, censoring and detection effort. A prospective test should hold out sites or cohorts, freeze predictions before revealing outcomes, evaluate calibration and Brier/log loss with censoring addressed, and distinguish fatal from nonfatal attacks.

**Limitations:** No measured independent predictive accuracy for actual predation; no claim of a resolved biological mystery or a priority discovery. Researcher feedback requested. AI-assisted analysis should be independently checked before publication.
