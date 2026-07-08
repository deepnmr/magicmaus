---
title: "magicmaus: certainty-bounded, score-committed assignment of methyl NMR resonances"
---

<!-- Build: pandoc manuscript.md --reference-doc=reference.docx -o magicmaus_AppNote.docx
     (reference.docx sets the Table style to 9 pt) -->

**Structural bioinformatics**

# magicmaus: certainty-bounded, score-committed assignment of methyl NMR resonances

Joonhyeok Choi ^1,2,†^, Sang-Yeop Lee ^3,4,†^, Sooho Koh ^1^, Minjune Yang ^1^, Kyoung-Seok Ryu ^1,4^, Donghan Lee ^1,4,\*^

^1^ Center for Protein Structure and Drug Mechanism Research, Korea Basic Science Institute, Cheongju, 28119, Republic of Korea

^2^ Department of Applied Pharmacy, College of Pharmacy, Chungbuk National University, Cheongju, 28160, Republic of Korea

^3^ Center for Target-to-Therapeutics Research, Korea Basic Science Institute, Cheongju, 28119, Republic of Korea

^4^ Department of Convergent Analytical Science, University of Science and Technology, Daejeon, 34113, Republic of Korea

^\*^ Corresponding author. Donghan Lee, Center for Protein Structure and Drug Mechanism Research, Korea Basic Science Institute, Cheongju, 28119, Republic of Korea; Department of Convergent Analytical Science, University of Science and Technology, Daejeon, 34113, Republic of Korea (dlee04@kbsi.re.kr)

^†^ These authors contributed equally.

Associate Editor: (to be assigned)

---

## Abstract

**Summary:** Automated assignment of methyl NMR probes in large proteins is
approached in two irreconcilable ways: as a constraint-satisfaction problem
(MAUS), which returns for every peak the set of methyls consistent with all hard
constraints and provably never excludes the correct one but abstains on any
residual degeneracy; and as a scoring problem (MAGIC), which commits to a single
best assignment but, searching the full candidate space, does so on a near-flat
objective and is frequently wrong. We present **magicmaus**, a hybrid that uses
the satisfiability layer to bound the search to a certifiably truth-containing set
of per-peak candidates, then maximises an intensity-weighted NOE score *within*
those bounds — by a diverse-seed multistart that reaches the objective's optimum
where a single-seed local search traps 10–20% short — to commit one globally
coherent, injective assignment with a per-peak confidence tier. We score accuracy
residue-wise (correct residue, independent of the prochiral δ1/δ2 or γ1/γ2 methyl,
which an achiral NOE network cannot orient) because the residue is the assignment
problem's hard part. Across seven structure-simulated benchmark targets (43–257
methyls) built from real shifts and structures with realistic measurement noise,
magicmaus commits a single answer for every peak and is the most accurate of the
three engines on every target it can bound tightly — 56.5–100% residue-wise, versus
0–79% for full-space scoring and 39–84% for the constraint layer forced to guess —
while retaining the constraint layer's **100% never-exclude** guarantee as an
explicit ambiguity envelope. On maltose-binding protein (192 methyls) it assigns
60.9% of peaks to the correct residue versus 34.4% for scoring alone. Malate
synthase G (257 methyls) marks the method's scale limit and is reported honestly:
its exact option-set enumeration is intractable, and over the tractable but loose
fallback bounds its sparse (H)CCH network leaves magicmaus (14.8%) below full-space
scoring (30.0%). A real-experimental multimer — the TNF-α homotrimer, measured HMQC,
3D (H)CCH NOESY-HMQC and HMBC-HMQC peak lists scored against an AlphaFold3 model —
demonstrates the method on genuine spectra and honestly bounds it: treating the
three chains as symmetry images is required to explain the inter-subunit NOEs, the
never-exclude envelope is 96.5% (real firm edges are occasionally wrong against a
predicted fold), and magicmaus is still the best engine at 30.6% residue-wise
(versus 20.0% scoring, 10.6% constraint-forced), the gap to its 24.7% methyl-level
figure being entirely unresolved geminal swaps.

**Availability and implementation:** magicmaus is implemented in Python 3 (NumPy,
PySAT) and released under the MIT license at
<https://github.com/deepnmr/magicmaus>, with tutorials and the benchmark data.

**Contact:** dlee04@kbsi.re.kr

---

## 1 Introduction

Methyl groups are the workhorse probes of solution NMR on large proteins and
molecular machines. Selective ^13^CH~3~ labelling of Ile/Leu/Val (and Ala, Met,
Thr) methyls on a perdeuterated background (Tugarinov and Kay, 2004), combined
with methyl-TROSY spectroscopy, yields sharp signals well beyond 100 kDa and up
to the ~1 MDa regime (Tugarinov *et al.*, 2003; Rosenzweig and Kay, 2014).
Exploiting these probes requires assigning each observed methyl ^1^H/^13^C HMQC
cross peak to a specific methyl in the protein structure — a laborious,
error-prone step widely regarded as the principal bottleneck of methyl NMR, which
has motivated a family of automated, structure-based methods that match an
experimental methyl–methyl NOE network onto the contact graph of a known
structure (reviewed in Pritišanac *et al.*, 2020).

These methods differ chiefly in how they treat residual ambiguity, and magicmaus
builds on one representative of each of two paradigms. **MAGIC** (Monneau *et
al.*, 2017) scores a global peak→methyl map with a confidence-weighted NOE
objective and returns the highest scorer — a single answer per peak. **MAUS**
(Nerli *et al.*, 2021) — following the graph-matching formulation of MAGMA
(Pritišanac *et al.*, 2017) but recast as Boolean satisfiability and solved with
a SAT solver — returns for each peak the *set* of methyls consistent with every
hard constraint, provably never excluding the correct assignment. Other tools
occupy intermediate points but demand richer experimental input; MethylFLYA
(Pritišanac *et al.*, 2019), for instance, reaches ~1% error by consensus over
many statistical runs on multi-spectrum data.

The two have complementary failure modes. MAUS is safe but indecisive: where the
data do not force a unique answer — geminal methyl pairs, shift-degenerate peaks —
it reports the full option set and commits to nothing. MAGIC is decisive but,
searching the entire type-matched candidate space, sits on an objective whose
runner-up assignments lie within ~1–2% of the true optimum on structure-derived
NOE data (Monneau *et al.*, 2017), so its bounded search commits to
near-optimal-but-wrong answers for most peaks. Neither is dominant. We show that
the natural synthesis — bound the space with certainty, then score within the
residual degeneracy — is both simple and, where the bounds can be made tight,
decisively better than either parent.

## 2 Approach

magicmaus runs the two layers in sequence (Fig. 1A). Its inputs are a 2D methyl
HMQC peak list, a 3D (H)CCH NOESY-HMQC list and an optional 3D HMBC-HMQC list
(both given as `label C2 C1 H1` triples — the detected methyl at (C1, H1), the
partner methyl by carbon C2 only), and a structure in mmCIF or PDB; a homo-oligomer
is parsed with all chains retained as symmetry images of each methyl.

**Layer 1 (bound).** The MAUS SAT encoding maps each HMQC peak to a candidate
methyl of the matching residue type; hard clauses enforce exactly-one assignment
per peak, injectivity, and that every firm NOE cross peak lands on a structural
contact within the distance cutoff (an HMBC-HMQC geminal link, when supplied,
additionally ties the two linked peaks to the two geminal methyls of one residue).
For each peak the set of methyls appearing in at least one satisfying assignment is
enumerated with the solver's assumption interface (Glucose, Audemard and Simon,
2009, via PySAT, Ignatiev *et al.*, 2018), accelerated by a unit-propagation
pre-filter that discards candidates refuted without a full solve. Because the
correct global assignment is itself a model, the true methyl is present in every
peak's option set — the *never-exclude* guarantee — and the option sets prune each
peak's candidates from up to ~60 to typically 1–3. This is also why the tempting
shortcut of assigning a peak to the last unclaimed methyl of its type is unsound:
the structure routinely carries many more methyls of a type than are observed, so
an unclaimed methyl may simply be unobserved rather than the peak's answer. Because
MAUS models a methyl as free to go unassigned, its enumeration performs every
elimination injectivity licenses — but no more — and thereby keeps the truth in the
set instead of forcing a wrong pin.

**Layer 2 (commit).** The per-peak option sets are enumerated *independently*, so
their product is not jointly realizable and does not identify a single best map.
That is the question a score answers. Restricting every peak's domain to its option
set — which removes only methyls appearing in no satisfying map, preserving the
solution space and the truth — magicmaus maximises a MAGIC-style objective over the
pruned domains: each firm NOE contributes intensity·(1/r^6^) for the structural
contact it is placed on, so a strong NOE is driven onto a close contact, and the
score is maximised over injective, NOE-consistent assignments. On a network
simulated from the true structure this objective's global optimum is the truth (a
truth-seeded ascent stays there), but the truth is a strong yet *narrow* optimum in
a rugged landscape: a single feasible seed's greedy ascent — and, we found,
simulated annealing from one seed — settles 10–20% below it, trapped because the
hard NOE constraints block the moves that reach it. magicmaus instead runs a
**diverse-seed multistart**. Many independent feasible seeds are drawn by
re-solving the option-set SAT with randomised variable phases; each is polished by
a feasibility-preserving greedy ascent (relocate and swap moves that keep the map
injective and every firm NOE on a contact); and the highest-objective assignment
over all restarts is committed. Because Layer 1's pruning collapses each peak's
domain to 1–3 candidates, a few dozen diverse restarts reach the objective's
optimum inside the truth-containing space — the decisive difference from MAGIC,
whose un-pruned domains leave astronomically many near-flat basins the same
multistart cannot escape. The whole pipeline completes in seconds for 192 methyls.

We report accuracy **residue-wise**: a committed call is correct if it names the
right residue, independent of which of the two prochiral methyls (Leu δ1/δ2, Val
γ1/γ2) it assigns. The geminal orientation is a near-symmetric coin flip an achiral
NOE network cannot settle for pairs whose two methyls make similar contacts, and
the residue — not the individual methyl — is the assignment problem's hard part and
the unit that downstream analysis consumes. We additionally report methyl-level
(atom-exact) accuracy to expose the residual geminal swaps (Table 3).

Each peak is reported three ways (Fig. 1A): the committed single call; its MAUS
option set, retained as an explicit ambiguity envelope; and a confidence tier from
the local scoring margin — **unique** (forced by hard constraints), **scored**
(resolved by the NOE score), or **ambiguous** (a tied alternative exists, i.e. a
genuine symmetry). Optionally, the ambiguous NOE cross peaks Layer 1 discards are
folded back in as diluted, intensity-weighted soft evidence (`--soft-ambiguous`).

## 3 Results

**Benchmark construction.** For the seven simulated targets, all inputs are built
from a PDB (or mmCIF) structure and the matching BMRB chemical-shift deposition:
methyl carbon shifts are the deposited carbon values, methyl proton shifts the mean
of the three methyl protons. To model realistic measurement scatter, each methyl's
shift is drawn once from a normal distribution about its deposited value —
σ = 0.02 ppm for ^1^H and 0.10 ppm for ^13^C — and this one measured shift is used
consistently across all experiments, so the frequency degeneracy that drives
assignment ambiguity is realistic. Each observed methyl becomes an anonymous HMQC
peak (P1…Pn) whose structural identity is written only to a separate truth key.
Because BMRB deposits no NOESY, a methyl–methyl (H)CCH NOESY network (Rossi *et
al.*, 2016; Wen *et al.*, 2012) is simulated from the structure (a cross peak for
every methyl pair within 8 Å, both directions) with 1/r^6^ intensities; the
identical network is supplied to all three engines. We benchmark seven targets:
ubiquitin (BMRB 6457, PDB 1UBQ), *E. coli* maltose-binding protein (MBP; BMRB 7114,
PDB 1ANF; 192 methyls), the four de-novo blind targets of the MAUS study — IL-2 and
the HNH, REC2 and REC3 domains of *S. pyogenes* Cas9 — and malate synthase G (MSG;
PDB 1D8C; 257 methyls, the classic large-protein methyl benchmark), spanning 43–257
observed methyls (Table 1). To test the method on genuine spectra we add one
real-experimental target: measured AILTV methyl HMQC, 3D (H)CCH NOESY-HMQC and
HMBC-HMQC peak lists of the tumour-necrosis-factor-α homotrimer (85 assigned
methyls), scored against an AlphaFold3 model of the trimer. This is also our only
multimer; because the three protomers are magnetically equivalent, each residue
gives one HMQC peak whose NOE partners may lie in the same or a neighbouring
subunit, so the structure is parsed with all three chains retained as symmetry
images and every contact scored by the minimum distance over subunit pairs.

**The three engines, residue-wise (Table 1).** MAGIC, scoring over the full
type-matched space, assigned 0–79% of peaks to the correct residue: it does well
where the space is small and the intensities clean (ubiquitin 79.1%) but collapses
as the space grows and the network sparsens (REC2 1.6%, REC3 0.0%), the near-flat
full-space landscape leaving its multistart in a wrong basin. MAUS, forced to
commit a single residue where its option set does not pin one (it cannot rank
within a set), reached 2.3–84.2%; its envelope held the truth for 100% of peaks on
every simulated target. magicmaus committed a single answer for every peak while
preserving that **100%** envelope, and was the most accurate of the three on every
target whose bounds it could enumerate tightly: **56.5–100%** residue-wise, a
perfect 43/43 on ubiquitin and **60.9%** on MBP — beating full-space scoring on
seven of eight targets and the constraint-forced baseline on all eight. This comes
from the multistart search: on a simulated intensity network built from the true
structure the NOE objective's global optimum is the truth (a truth-seeded ascent
scores ~100% on ubiquitin), and the diverse-seed multistart reaches it inside the
pruned space where a single-seed ascent or annealer stalls 10–20% short. The
optional soft-ambiguous evidence is a genuine but two-sided lever: it helps where
the discarded NOEs carry residue information (IL-2 64.4→86.4, HNH 86.0→93.0, MBP
60.9→64.6) and hurts where they mostly add noise (REC2 71.4→58.7, REC3 56.5→51.8,
TNF-α 30.6→28.2), so it is offered as an option, not a default, and the tables
below report the base configuration.

**The scale limit (MSG).** MSG is reported honestly as the point where the method's
bounding step becomes impractical. Its Leu/Val domains are 138 methyls wide and the
simulated (H)CCH network resolves only ~85 cross peaks to a firm constraint, so
(i) the exact per-candidate option enumeration does not finish within a 20-min
budget even with the propagation pre-filter, and (ii) the tractable fall-back —
scoring over the arc-consistency-pruned domains, which stay wide — leaves the bounds
too loose and the network too sparse for the scoring layer to disambiguate:
magicmaus commits 14.8%, *below* full-space MAGIC's 30.0%, though still under the
100% envelope. MSG thus delimits where the synthesis pays off: magicmaus's advantage
requires either tractable tight bounds or a dense-enough NOE network, and MSG at
this scale offers neither. This is a property of the method, not a defect of the
implementation, and we state it rather than omit the target.

**Residue type shapes difficulty (Table 2).** The labeling scheme shapes the problem
in two ways the type-resolved accuracy makes explicit. First, the number of labeled
types sets the granularity of the candidate partition: a peak competes only against
methyls of its own residue type, so the two AILMTV targets (ubiquitin, MBP) spread
their peaks over six types while a three-type ILV set concentrates competition.
Second, labeling sets the irreducible geminal load: Leu and Val carry prochiral
methyl pairs, whereas Ile contributes a single δ1 methyl. Ile is accordingly the
most reliably assigned type — 100% on four targets and 73–89% on the mid-size ones —
because it is single-methyl and, being typically buried, NOE-rich; Leu and Val are
the bottleneck, at 100% only where the network is dense (ubiquitin) and falling to
~30–70% on the sparser REC3/MBP networks. Being single-methyl is necessary but not
sufficient: on MBP, Ala Cβ (43%) and Thr Cγ2 (60%) are among the least accurate
types despite no prochiral degeneracy, because these often surface-exposed methyls
simply make too few NOEs to be pinned — whereas Ile on the same protein reaches 73%.
Difficulty is therefore the interaction of labeling with NOE information content,
not labeling alone.

**The residual error is a geminal swap (Table 3).** Resolving accuracy to the
individual prochiral methyl confirms that the residue-vs-methyl gap for Leu/Val is a
geminal swap, not random misassignment: on nearly every target the two members of
each pair degrade in near-lockstep — REC2 Leu δ1/δ2 at 58%/62%, REC3 at 48%/44%, HNH
Val γ1/γ2 both 50% — the signature of the achiral network placing the pair on the
right *residue* but the wrong *methyl*. Both members stay inside the MAUS envelope
(the truth is never excluded), so the swap is a calibrated coin flip the confidence
tier flags as `ambiguous`; only a signal able to tell δ1 from δ2 — the intensity
score where the two methyls' contacts differ, or an independent stereospecific
assignment where they do not — can break it. This is exactly why residue-wise is the
meaningful metric: the methyl-level totals (Table 1, "methyl") sit below the
residue-wise ones by precisely the unresolved-swap fraction (ubiquitin 95.3 vs 100,
MBP 56.8 vs 60.9, HNH 80.7 vs 86.0).

**MAUS is more residue-decisive than it looks.** MAUS's decisive fraction counts
only peaks pinned to a single *methyl*, booking every geminal-unresolved pair as an
abstention — yet such a pair is already decisive at the residue level. Collapsing
each option set to its residues, MAUS is residue-decisive (and, by never-exclude,
correct) on far more peaks than its methyl-unique count: 55.8% on ubiquitin, 68.4%
on HNH, 44.4% on REC2, versus low-single-digit methyl-unique fractions. The
complement — option sets still spanning more than one residue — is the *cross-residue*
ambiguity that dominates the Leu-crowded ILV targets (IL-2, MSG) and that an
HMBC-HMQC geminal link collapses by tying each pair to one residue.

**Real-experimental TNF-α.** On the real homotrimer, matched at H±0.02/C±0.10, the
method behaves consistently with the simulated targets while exposing two properties
only genuine data reveals. First, multimer handling is load-bearing: the measured
NOEs include inter-subunit contacts that have no structural home unless the three
chains are treated as symmetry images (contact = minimum distance over subunits),
and the reciprocal (H)CCH rows resolve both NOE ends by full (C, H) so every firm
edge is a real methyl–methyl contact. Second, the never-exclude guarantee is
conditional on the NOEs being consistent with the structure — a condition simulated
data satisfies by construction but a *predicted* fold does not: three truths carry a
firm NOE the AlphaFold3 model does not support at the cutoffs and fall out, so the
envelope is 96.5%, not 100%. No peak is uniquely pinned on this sparse, wide-domain
target, so the constraint-forced baseline is near-random (10.6%); full-space scoring
reaches 20.0%; and magicmaus is the best of the three at **30.6%** residue-wise. Its
gap to the 24.7% methyl-level figure is entirely geminal swaps (Table 3, TNF-α: Ile
75% but Leu/Val ≤31%), the orientation the global objective cannot fix here because,
scored against a predicted fold, its optimum no longer sits exactly at the truth.

**Table 1.** Residue-wise accuracy (correct residue, geminal δ1/δ2 orientation
ignored) on seven structure-simulated targets plus one real-experimental multimer
(TNF-α), all engines scored on the same 1/r^6^ intensity NOESY. MAGIC = full-space
scoring; MAUS = constraint layer forced to commit a single residue (it cannot rank
within an option set); magicmaus = base configuration (HMQC + 3D NOESY + 3D HMBC);
+soft = with soft-ambiguous evidence; methyl = magicmaus base at the atom-exact
level; Envelope = fraction of peaks whose MAUS option set contains the true residue
(never-exclude). MSG† uses the arc-consistency-pruned bounds (exact enumeration
intractable at 257 methyls); TNF-α‡ is the real-experimental, only multimeric target.

| Target | BMRB / PDB | Labeling | Methyls | MAGIC | MAUS | magicmaus | +soft | methyl | Envelope |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Ubiquitin | 6457 / 1UBQ | AILMTV | 43 | 79.1% | 72.1% | **100.0%** | 100.0% | 95.3% | 100% |
| IL-2 | 28104 / 1M47 | ILV | 59 | 54.2% | 39.0% | **64.4%** | 86.4% | 57.6% | 100% |
| HNH (Cas9) | 27949 / 6O56 | AILTV | 57 | 54.4% | 84.2% | **86.0%** | 93.0% | 80.7% | 100% |
| REC2 (Cas9) | 28105 / 4CMP | ILV | 63 | 1.6% | 54.0% | **71.4%** | 58.7% | 68.3% | 100% |
| REC3 (Cas9) | 28110 / 4ZT0 | ILV | 85 | 0.0% | 41.2% | **56.5%** | 51.8% | 50.6% | 100% |
| MBP | 7114 / 1ANF | AILMTV | 192 | 34.4% | 35.4% | **60.9%** | 64.6% | 56.8% | 100% |
| MSG† | SI / 1D8C | ILV | 257 | **30.0%** | 2.3% | 14.8% | 17.9% | 13.6% | 100% |
| TNF-α‡ | real / AF3 trimer | AILTV | 85 | 20.0% | 10.6% | **30.6%** | 28.2% | 24.7% | 96.5% |

†MSG methyl shifts are digitised from the reference assignment table of Pritišanac
*et al.* (2019); its exact option-set enumeration is intractable, so the bounds are
the arc-consistency-pruned domains — still 138 wide — and with only ~85 firm NOEs
the sparse network leaves magicmaus below full-space scoring, the method's scale
limit. ‡TNF-α: genuine methyl-NMR HMQC/NOESY/HMBC peak lists of the homotrimer,
scored against an AlphaFold3 model, matched at H±0.02/C±0.10 like the rest. Its
reciprocal (H)CCH rows resolve both NOE ends by full (C, H); treating the three
chains as symmetry images explains the inter-subunit NOEs; three truths carry an NOE
the *predicted* structure cannot support, so the envelope is 96.5%. No peak is
uniquely forced, so the MAUS column is entirely arbitrary-tiebreak.

**Table 2.** Residue-wise accuracy of the committed magicmaus call (base) resolved
by residue type. Each cell is correct/observed for that type; a dash marks a type
absent from the target's labeling. Ile leads (single-methyl, NOE-rich); Leu/Val
carry the geminal degeneracy; Ala/Thr, though single-methyl, are NOE-poor on MBP.

| Target | Ile | Leu | Val | Ala | Thr | Met | Total |
|---|---|---|---|---|---|---|---|
| Ubiquitin | 100% (7/7) | 100% (18/18) | 100% (8/8) | 100% (2/2) | 100% (7/7) | 100% (1/1) | 100% (43/43) |
| IL-2 | 78% (7/9) | 60% (25/42) | 75% (6/8) | — | — | — | 64% (38/59) |
| HNH | 100% (7/7) | 93% (26/28) | 62% (10/16) | 100% (2/2) | 100% (4/4) | — | 86% (49/57) |
| REC2 | 89% (8/9) | 65% (31/48) | 100% (6/6) | — | — | — | 71% (45/63) |
| REC3 | 54% (7/13) | 50% (25/50) | 73% (16/22) | — | — | — | 56% (48/85) |
| MBP | 73% (16/22) | 62% (37/60) | 72% (29/40) | 43% (19/44) | 60% (12/20) | 67% (4/6) | 61% (117/192) |
| MSG | 22% (9/41) | 14% (18/133) | 13% (11/83) | — | — | — | 15% (38/257) |
| TNF-α | 75% (6/8) | 31% (11/36) | 19% (5/26) | 33% (4/12) | 0% (0/3) | — | 31% (26/85) |

**Table 3.** Methyl-level (atom-exact) accuracy of the committed magicmaus call
(base) resolved to the individual methyl carbon. Geminal partners (Leu δ1/δ2, Val
γ1/γ2) are listed separately; their near-equal columns are the geminal-swap
signature — right residue, wrong prochiral methyl — and the reason the residue-wise
totals of Table 2 exceed these. A dash marks a methyl absent from the labeling.

| Target | Ile δ1 | Leu δ1 | Leu δ2 | Val γ1 | Val γ2 | Ala β | Thr γ2 | Met ε |
|---|---|---|---|---|---|---|---|---|
| Ubiquitin | 100% (7/7) | 89% (8/9) | 89% (8/9) | 100% (4/4) | 100% (4/4) | 100% (2/2) | 100% (7/7) | 100% (1/1) |
| IL-2 | 78% (7/9) | 52% (11/21) | 52% (11/21) | 50% (2/4) | 75% (3/4) | — | — | — |
| HNH | 100% (7/7) | 93% (13/14) | 86% (12/14) | 50% (4/8) | 50% (4/8) | 100% (2/2) | 100% (4/4) | — |
| REC2 | 89% (8/9) | 58% (14/24) | 62% (15/24) | 100% (3/3) | 100% (3/3) | — | — | — |
| REC3 | 54% (7/13) | 48% (12/25) | 44% (11/25) | 55% (6/11) | 64% (7/11) | — | — | — |
| MBP | 73% (16/22) | 47% (14/30) | 63% (19/30) | 65% (13/20) | 60% (12/20) | 43% (19/44) | 60% (12/20) | 67% (4/6) |
| MSG | 22% (9/41) | 10% (7/67) | 14% (9/66) | 7% (3/42) | 17% (7/41) | — | — | — |
| TNF-α | 75% (6/8) | 17% (3/18) | 17% (3/18) | 8% (1/13) | 31% (4/13) | 33% (4/12) | 0% (0/3) | — |

## 4 Conclusion

magicmaus shows that constraint satisfaction and scoring are not competing solutions
to methyl assignment but complementary stages of one: SAT to bound the answer with
certainty, an intensity-weighted NOE score to commit within the bound. The bounding
step is what makes the scoring tractable and accurate — a diverse-seed multistart
reaches the objective's optimum inside the small pruned space where the same search
over the full candidate space cannot — and where the bounds can be enumerated
tightly magicmaus is the most accurate of the three engines, at no cost to the
never-exclude envelope, which it always returns as an explicit ambiguity map. We
report accuracy residue-wise because the residue is the hard part and the geminal
δ1/δ2 orientation is a near-symmetric coin flip an achiral NOE network cannot settle;
the methyl-level gap is exactly that unresolved swap. Two honest boundaries frame the
method: at scale (malate synthase G, 257 methyls) the exact bounds become intractable
and, over the loose fallback, a sparse network leaves the synthesis below plain
scoring; and on real spectra scored against a *predicted* fold, the never-exclude
guarantee — not the committed call — is what survives a structure the data
contradict (a 96.5%, not 100%, envelope on TNF-α). The two clean-room parent
implementations and magicmaus are released together to support reuse and further
hybridization.

## Acknowledgements

The authors gratefully acknowledge Eun-Hee Kim and Dr. Hae-Kap Cheong from the
Korea Basic Science Institute (KBSI) for their valuable support and discussions,
and thank the maintainers of PySAT and the BMRB.

## Conflict of interest

None declared.

## Funding

This work was supported by the Korea Basic Science Institute (KBSI) internal
research programs [C612120, C623200, C623300, C625130, and C539110].

## Data availability

The program underlying this article is available on the project website at
<https://github.com/deepnmr/magicmaus>.

## References

Audemard,G. and Simon,L. (2009) Predicting learnt clauses quality in modern SAT
solvers. In *Proc. 21st Int. Joint Conf. on Artificial Intelligence (IJCAI)*,
pp. 399–404.

Ignatiev,A., Morgado,A. and Marques-Silva,J. (2018) PySAT: a Python toolkit for
prototyping with SAT oracles. In *Theory and Applications of Satisfiability
Testing (SAT 2018)*, LNCS 10929, pp. 428–437.

Monneau,Y.R. *et al.* (2017) Automatic methyl assignment in large proteins by the
MAGIC algorithm. *J. Biomol. NMR*, **69**, 215–227.

Nerli,S., De Paula,V.S., McShan,A.C. and Sgourakis,N.G. (2021)
Backbone-independent NMR resonance assignments of methyl probes in large
proteins. *Nat. Commun.*, **12**, 691.

Pritišanac,I. *et al.* (2017) Automatic assignment of methyl-NMR spectra of
supramolecular machines using graph theory. *J. Am. Chem. Soc.*, **139**,
9523–9533.

Pritišanac,I., Würz,J.M., Alderson,T.R. and Güntert,P. (2019) Automatic
structure-based NMR methyl resonance assignment in large proteins. *Nat.
Commun.*, **10**, 4922.

Pritišanac,I., Alderson,T.R. and Güntert,P. (2020) Automated assignment of methyl
NMR spectra from large proteins. *Prog. Nucl. Magn. Reson. Spectrosc.*,
**118–119**, 54–73.

Rosenzweig,R. and Kay,L.E. (2014) Bringing dynamic molecular machines into focus
by methyl-TROSY NMR. *Annu. Rev. Biochem.*, **83**, 291–315.

Rossi,P., Xia,Y., Khanra,N., Veglia,G. and Kalodimos,C.G. (2016) ^15^N and
^13^C-SOFAST-HMQC editing enhances 3D-NOESY sensitivity in highly deuterated,
selectively [^1^H,^13^C]-labeled proteins. *J. Biomol. NMR*, **66**, 259–271.

Siemons,L., Mackenzie,H.W., Shukla,V.K. and Hansen,D.F. (2019) Intra-residue
methyl–methyl correlations for valine and leucine residues in large proteins from
a 3D-HMBC-HMQC experiment. *J. Biomol. NMR*, **73**, 749–757.

Spurlino,J.C., Lu,G.-Y. and Quiocho,F.A. (1991) The 2.3-Å resolution structure of
the maltose- or maltodextrin-binding protein, a primary receptor of bacterial
active transport and chemotaxis. *J. Biol. Chem.*, **266**, 5202–5219.

Tugarinov,V. and Kay,L.E. (2004) An isotope labeling strategy for methyl TROSY
spectroscopy. *J. Biomol. NMR*, **28**, 165–172.

Tugarinov,V., Hwang,P.M., Ollerenshaw,J.E. and Kay,L.E. (2003) Cross-correlated
relaxation enhanced ^1^H–^13^C NMR spectroscopy of methyl groups in very high
molecular weight proteins and protein complexes. *J. Am. Chem. Soc.*, **125**,
10420–10428.

Ulrich,E.L. *et al.* (2008) BioMagResBank. *Nucleic Acids Res.*, **36**,
D402–D408.

Wen,J., Zhou,P. and Wu,J. (2012) Efficient acquisition of high-resolution 4-D
diagonal-suppressed methyl–methyl NOESY for large proteins. *J. Magn. Reson.*,
**218**, 128–132.

---

![**Fig. 1.** (**A**) The magicmaus pipeline. A SAT layer (MAUS) bounds each peak
to a per-peak option set that provably contains the truth (100% envelope on
simulated data) and prunes candidates from up to ~60 to 1–3; an intensity-weighted
NOE score (MAGIC-style), applied only within those bounds by a diverse-seed
multistart (independent SAT-feasible seeds, each greedy-ascended), commits to a
single coherent injective map with a per-peak confidence tier (unique / scored /
ambiguous). Accuracy is scored residue-wise, the geminal δ1/δ2 orientation being a
near-symmetric coin flip. (**B**) Residue-wise accuracy across the seven
structure-simulated targets (43–257 methyls, ordered by size) plus the
real-experimental TNF-α homotrimer (*real): magicmaus (blue) versus full-space MAGIC
(red) and the MAUS constraint layer forced to a single residue (purple). Green
markers denote the truth-in-envelope guarantee — 100% on every simulated target and
96.5% on real TNF-α data, where three peaks carry an NOE the predicted structure
cannot support. magicmaus is the most accurate engine on every target but MSG (257
methyls), where the exact bounds are intractable and the sparse network leaves the
loose-bounded score below full-space MAGIC.](figure1.png)
