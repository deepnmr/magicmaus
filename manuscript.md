---
title: "magicmaus: certainty-bounded, score-committed assignment of methyl NMR resonances"
---

**Structural bioinformatics**

# magicmaus: certainty-bounded, score-committed assignment of methyl NMR resonances

**Donghan Lee** ^1,\*^

^1^ Korea Basic Science Institute, Ochang, Republic of Korea.

^\*^ To whom correspondence should be addressed.

Associate Editor: (to be assigned)

---

## Abstract

**Summary:** Automated assignment of methyl NMR probes in large proteins is
approached in two irreconcilable ways: as a constraint-satisfaction problem
(MAUS), which returns for every peak the set of methyls consistent with all hard
constraints and provably never excludes the correct one, but abstains on any
residual degeneracy; and as a scoring problem (MAGIC), which commits to a single
best assignment but, over the full candidate space, does so on a near-flat
objective and is frequently wrong. We present **magicmaus**, a hybrid that uses
the satisfiability layer to bound the search to a certifiably truth-containing
set of per-peak candidates, then applies an intensity-weighted NOE score *within*
those bounds to commit to a single, globally coherent, injective assignment
carrying a per-peak confidence tier. On maltose-binding protein (192 methyls),
scored on an identical intensity NOESY, magicmaus assigns 72.9% of methyls
correctly (79.7% with ambiguous-NOE evidence) versus 5.7% for the scoring method
alone, while retaining the constraint method's 100% never-exclude guarantee as an
explicit ambiguity envelope.

**Availability and implementation:** magicmaus is implemented in Python 3 (NumPy,
PySAT) and released under the MIT license at
<https://github.com/deepnmr/magicmaus>, with tutorials and the benchmark data.

**Contact:** kbsi.bionmr@gmail.com

---

## 1 Introduction

Methyl groups are the workhorse probes of solution NMR on large proteins and
molecular machines: methyl-TROSY spectroscopy yields sharp signals well beyond
100 kDa (Tugarinov *et al.*, 2003). Exploiting them requires assigning each
observed methyl ^1^H/^13^C HMQC cross peak to a specific methyl in the protein
structure, a laborious step that has motivated several automated,
structure-based methods driven by methyl–methyl NOEs.

Two paradigms dominate. **MAGIC** (Monneau *et al.*, 2017) scores a global
peak→methyl map with a confidence-weighted NOE objective and returns the highest
scorer — a single answer per peak. **MAUS** (Nerli *et al.*, 2021) casts the same
problem as subgraph isomorphism and solves it with a SAT solver, returning for
each peak the *set* of methyls consistent with every hard constraint, and
provably never excluding the correct assignment.

The two have complementary failure modes. MAUS is safe but indecisive: where the
data do not force a unique answer — geminal methyl pairs, shift-degenerate peaks
— it reports the full option set and commits to nothing. MAGIC is decisive but,
searching the entire type-matched candidate space, sits on an objective whose
runner-up assignments lie within ~1–2% of the true optimum on structure-derived
NOE data (Monneau *et al.*, 2017), so its bounded search commits to
near-optimal-but-wrong answers for most peaks. Neither is dominant. We show that
the natural synthesis — bound the space with certainty, then score within the
residual degeneracy — is both simple and decisively better than either parent.

## 2 Approach

magicmaus runs the two layers in sequence (Fig. 1A).

**Layer 1 (bound).** The MAUS SAT encoding is reused verbatim: variables map each
HMQC peak to a candidate methyl of the matching residue type; hard clauses
enforce exactly-one assignment per peak, injectivity, and that every firm NOE
cross peak lands on a structural contact within the distance cutoff. For each
peak the set of methyls appearing in at least one satisfying assignment is
enumerated with the solver's assumption interface (Glucose; Audemard and Simon,
2009). Because the correct global assignment is itself a model, the true methyl
is present in every peak's option set — the *never-exclude* guarantee — and the
option sets prune each peak's candidates from up to ~60 to typically 1–3.

**Layer 2 (commit).** The per-peak option sets are enumerated *independently*, so
their product is not jointly realizable and does not identify a single best map.
This is precisely the question a score answers. Restricting every peak's domain
to its option set — which removes only methyls appearing in no satisfying map,
and so preserves the solution space and the truth — magicmaus obtains one
jointly-consistent assignment from the SAT solver and refines it by
feasibility-preserving coordinate ascent on a MAGIC-style objective: each firm
NOE contributes intensity·(1/r^6^) for the structural contact it is placed on, so
a strong NOE is driven onto a close contact. Every move stays injective and
NOE-consistent, so the output is always a valid bijection — a property a naïve
per-cluster search cannot guarantee on the single 138-peak degeneracy cluster of
the benchmark below. Optionally, the ambiguous NOE cross peaks MAUS discards are
folded back in as diluted, intensity-weighted soft evidence.

Each peak is reported three ways (Fig. 1A): the committed single call; its MAUS
option set, retained as an explicit ambiguity envelope; and a confidence tier
from the local scoring margin — **unique** (forced by hard constraints),
**scored** (resolved by the NOE score), or **ambiguous** (a tied alternative
exists, i.e. a genuine symmetry). Because the scored search runs only over the
tiny pruned domains, the whole pipeline completes in ~1 s for 192 methyls.

## 3 Results

We benchmarked on *E. coli* maltose-binding protein (MBP) using experimental
methyl ^13^C/^1^H shifts from BMRB entry 7114 and structure PDB 1ANF (192
methyls: 44 Ala, 22 Ile, 60 Leu, 6 Met, 20 Thr, 40 Val). As BMRB deposits no
NOESY peak list, a 3D (H)CCH network was simulated from 1ANF geometry with
1/r^6^ intensities; the identical network was supplied to all three engines
(Fig. 1B).

The scoring method (MAGIC) assigned 5.7% of methyls correctly (10.4% at
residue level), consistent with its reported 4–10% on structure-simulated NOESY
and with the near-flat-landscape limitation of full-space scoring. The
constraint method (MAUS) resolved 26.6% of peaks uniquely — all correct — and
abstained on the remaining 73%, with the truth in the option set for 100% of
peaks; its result was unchanged by the intensity column, as its boolean
constraints cannot use it. magicmaus committed a single answer for all 192 peaks
at **72.9%** correct (**79.2%** at residue level), rising to **79.7%** (85.4%)
when the discarded ambiguous NOEs were included, while preserving the **100%**
never-exclude envelope throughout. The residual errors are dominated by geminal
methyl pairs and shift-degenerate peaks — genuine symmetries an achiral NOE
network cannot resolve — which magicmaus flags as `ambiguous` and reports as full
option sets rather than guessing. Running the same scoring inside the MAUS bounds
thus yields ~13× the single-answer accuracy of scoring over the full space, at no
cost to the certainty guarantee.

## 4 Conclusion

magicmaus shows that constraint satisfaction and scoring are not competing
solutions to methyl assignment but complementary stages of one: SAT to bound the
answer with certainty, an intensity-weighted NOE score to commit within the
bound. The result is a single, coherent assignment with calibrated per-peak
confidence and an always-correct ambiguity envelope, obtained in seconds. The
scoring layer is intensity-aware, so the approach improves directly with the
information content of the NOESY — from a boolean network to real intensities to
4D experiments — and accepts tentative anchors that propagate through both
layers. The two clean-room parent implementations and magicmaus are released
together to support reuse and further hybridization.

## Acknowledgements

The author thanks the maintainers of PySAT and the BMRB.

## Funding

None declared.

*Conflict of Interest:* none declared.

## References

Audemard,G. and Simon,L. (2009) Predicting learnt clauses quality in modern SAT
solvers. In *Proc. IJCAI*, pp. 399–404.

Ignatiev,A., Morgado,A. and Marques-Silva,J. (2018) PySAT: a Python toolkit for
prototyping with SAT oracles. In *Proc. SAT*, pp. 428–437.

Monneau,Y.R. *et al.* (2017) Automatic methyl assignment in large proteins by the
MAGIC algorithm. *J. Biomol. NMR*, **69**, 215–227.

Nerli,S., De Paula,V.S., McShan,A.C. and Sgourakis,N.G. (2021)
Backbone-independent NMR resonance assignments of methyl probes in large
proteins. *Nat. Commun.*, **12**, 691.

Tugarinov,V., Hwang,P.M., Ollerenshaw,J.E. and Kay,L.E. (2003) Cross-correlated
relaxation enhanced ^1^H–^13^C NMR spectroscopy of methyl groups in very high
molecular weight proteins and protein complexes. *J. Am. Chem. Soc.*, **125**,
10420–10428.

---

![**Fig. 1.** (**A**) The magicmaus pipeline. A SAT layer (MAUS) bounds each peak
to a per-peak option set that provably contains the truth (100% envelope) and
prunes candidates from up to ~60 to 1–3; an intensity-weighted NOE score
(MAGIC-style), applied only within those bounds via a SAT-feasible seed and
feasibility-preserving coordinate ascent, commits to a single coherent map with a
per-peak confidence tier (unique / scored / ambiguous). (**B**) Methyl-level
accuracy on maltose-binding protein (192 methyls), all engines scored on the same
1/r^6^ intensity NOESY. MAUS is decisive on 26.6% of peaks (rest abstain); MAGIC
and magicmaus commit on all 192. Green markers denote the 100% truth-in-envelope
guarantee, preserved by MAUS and magicmaus.](figure1.png)
