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
constraints and provably never excludes the correct one, but abstains on any
residual degeneracy; and as a scoring problem (MAGIC), which commits to a single
best assignment but, over the full candidate space, does so on a near-flat
objective and is frequently wrong. We present **magicmaus**, a hybrid that uses
the satisfiability layer to bound the search to a certifiably truth-containing
set of per-peak candidates, then applies an intensity-weighted NOE score *within*
those bounds to commit to a single, globally coherent, injective assignment
carrying a per-peak confidence tier. Across seven structure-simulated benchmark
targets (43–257 methyls) built from real shifts and structures, magicmaus commits
a single answer for every peak at 29–100% methyl-level accuracy — up to an order
of magnitude above the scoring method — while retaining the constraint method's
100% never-exclude guarantee as an explicit ambiguity envelope. On
maltose-binding protein (192 methyls) it assigns 87.0% of methyls correctly
(87.5% with ambiguous-NOE evidence) versus 5.7% for scoring alone. A further
real-experimental multimer (the TNF-α homotrimer, scored against an AlphaFold3
model) demonstrates the method on genuine spectra and honestly bounds it: reciprocal
symmetric NOESY pairing reaches a 98.8% never-exclude envelope (one peak carries an
NOE the predicted structure cannot support), while the plain carbon-only match is
UNSAT on this dense trimer — conditional, unlike the simulated targets, on the
measured NOEs being consistent with the (predicted) structure — and treating
the trimer's chains as symmetry images is required to explain the inter-subunit
NOEs, recovering 6 further peaks in the envelope.

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
enumerated with the solver's assumption interface (Glucose, Audemard and Simon,
2009, via PySAT, Ignatiev *et al.*, 2018). Because the correct global assignment is itself a model, the true methyl
is present in every peak's option set — the *never-exclude* guarantee — and the
option sets prune each peak's candidates from up to ~60 to typically 1–3. This is
also why the tempting shortcut of assigning a peak to the last unclaimed methyl of
its type is unsound: the structure routinely carries many more methyls of a type
than are observed (REC3, for instance, has 292 Leu methyls but only 50 Leu peaks),
so an unclaimed methyl may simply be unobserved rather than the peak's answer.
Because MAUS models a methyl as free to go unassigned, its enumeration already
performs every elimination that injectivity licenses — but no more — and thereby
keeps the truth in the set instead of forcing a wrong pin.

**Layer 2 (commit).** The per-peak option sets are enumerated *independently*, so
their product is not jointly realizable and does not identify a single best map.
This is precisely the question a score answers. Restricting every peak's domain
to its option set — which removes only methyls appearing in no satisfying map,
and so preserves the solution space and the truth — magicmaus obtains one
jointly-consistent assignment from the SAT solver and refines it by
feasibility-preserving search on a MAGIC-style objective: each firm NOE
contributes intensity·(1/r^6^) for the structural contact it is placed on, so a
strong NOE is driven onto a close contact. A plain greedy ascent settles into the
nearest local optimum, which on the near-flat MAGIC landscape lies 10–20% below
the truth's objective; magicmaus instead runs simulated annealing over three
feasibility-preserving move classes — relocate, pairwise swap, and three-cycle
rotation — the last being essential, since a swap requires the displaced methyl to
lie in its partner's option set and stalls in the tightly coupled Leu/Val graphs
that rotations traverse. The best annealed state is polished by a final greedy
ascent. Every move stays injective and NOE-consistent, so the output is always a
valid bijection — a property a naïve per-cluster search cannot guarantee on the
single 138-peak degeneracy cluster of the benchmark below. Because this objective
is well-determined only when the NOESY carries real intensities and most optimised
peaks bear a firm NOE, the annealer is gated on both and reduces to the plain ascent
otherwise (a boolean network, or under 75% firm-NOE coverage). Optionally,
the ambiguous NOE cross peaks MAUS discards are folded back in as diluted,
intensity-weighted soft evidence.

Each peak is reported three ways (Fig. 1A): the committed single call; its MAUS
option set, retained as an explicit ambiguity envelope; and a confidence tier
from the local scoring margin — **unique** (forced by hard constraints),
**scored** (resolved by the NOE score), or **ambiguous** (a tied alternative
exists, i.e. a genuine symmetry). Because the scored search runs only over the
tiny pruned domains, the whole pipeline completes in ~0.3 s for 192 methyls.

## 3 Results

**Benchmark construction.** All inputs are built by `make_peaklists.py` from a
PDB structure and the matching BMRB chemical-shift deposition: methyl carbon
shifts are the deposited carbon values, methyl proton shifts the mean of the
three methyl protons, and each observed methyl becomes an anonymous HMQC peak
(P1…Pn) whose structural identity is written only to a separate truth key.
Because BMRB deposits no NOESY peak list, a methyl–methyl (H)CCH NOESY network
(Rossi *et al.*, 2016; Wen *et al.*, 2012) is simulated from the structure (a cross peak for every methyl
pair within 7.9 Å, both directions)
with 1/r^6^ intensities; the identical network is supplied to all three engines.
We benchmark seven targets: ubiquitin (BMRB 6457, PDB 1UBQ; the reference protein
of biomolecular NMR), *E. coli* maltose-binding protein (MBP; BMRB 7114, PDB
1ANF; Ulrich *et al.*, 2008; Spurlino *et al.*, 1991; 192 methyls), the four
de-novo blind targets of the MAUS study (Nerli *et al.*, 2021) — interleukin-2
and the HNH, REC2 and REC3 domains of *S. pyogenes* Cas9 — and malate synthase G
(MSG; PDB 1D8C), at 257 methyls the largest single-chain protein whose Ile/Leu/Val
methyls have been assigned by solution NMR and the classic large-protein methyl
benchmark, spanning 43–257 observed methyls (Table 1). To test the method on
genuine spectra rather than a simulated network we add one real-experimental
target: measured AILTV methyl HMQC, 3D (H)CCH NOESY-HMQC and HMBC-HMQC peak lists
of the tumour-necrosis-factor-α homotrimer (85 assigned methyls), scored against
an AlphaFold3 model of the trimer. This target is also our only multimer; because
the three protomers are magnetically equivalent, each residue gives one HMQC peak
whose NOE partners may lie in the same or a neighbouring subunit, so `parse_structure`
retains all three chains as symmetry images of each methyl and scores every
structural contact by the minimum distance over subunit pairs.

**Results.** MAGIC, scoring over the full type-matched space, assigned 6–12% of
methyls correctly where it converged (Table 1) — consistent with its reported
4–10% on structure-simulated NOESY and with the near-flat-landscape limitation of
full-space scoring; on the two Leu-dense Cas9 domains and on MSG it did not return
within a 15-min budget, itself illustrating the cost of unbounded scoring. MAUS
resolved 2–35% of peaks uniquely (all correct) and abstained on the rest, with the
truth in the option set for 100% of peaks on every target; its result was unchanged by
the intensity column, as its boolean constraints cannot use it. magicmaus
committed a single answer for every peak while preserving that **100%**
never-exclude envelope throughout, at **82–100%** methyl-level accuracy on the
smaller targets (a perfect 43/43 on ubiquitin) and **87.0%** on MBP — up to an
order of magnitude above scoring over the full space. This accuracy comes from the
scoring layer's 3-cycle simulated-annealing search: on an intensity network the
NOE objective's global optimum is the truth (a truth-seeded search scores ~96%),
and annealing reaches it where a plain greedy ascent stalls 10–20% short — the
rotational moves cross the tightly coupled Leu/Val option graphs that pairwise
swaps cannot. That objective is trustworthy only when it is well-determined, which
requires two conditions, and magicmaus withholds the annealer — falling back to the
plain greedy ascent — when either fails: (i) the NOESY must carry real intensities
that pin each contact to its distance (on a boolean network the annealer would
merely overfit to structural-contact density), and (ii) most of the peaks being
optimised must actually carry a firm NOE. The hard cases are the targets whose 3D
(H)CCH network yields few firm NOEs: REC3 (60.0%; 50 of 85 methyls Leu) and above
all MSG, where only 262 cross peaks resolve to a firm constraint so 95 of 257 peaks
carry none — below magicmaus's 75% firm-NOE coverage cut, so the annealer is
withheld and the greedy ascent commits 29.6% (38.5% with HMBC) where MAGIC does not
converge at all, still under the 100% envelope. On these, an achiral NOE network
leaves many geminal pairs and
shift-degenerate peaks genuinely unresolvable, which magicmaus flags as `ambiguous`
and reports as full option sets rather than guessing. Folding the discarded
ambiguous NOEs back in as soft evidence (`--soft-ambiguous`) helped on most targets
(IL-2 +8.5, HNH +3.5, REC2 +1.6 points) but was a wash on the densely degenerate
REC3, so it is offered as an option, not a default. Adding one further optional
input — an HMBC-HMQC experiment (`--hmbc`; Siemons *et al.*, 2019) that identifies
which two HMQC peaks are the geminal pair of one residue — enters as a hard
constraint tying that pair to a single structural residue. It does not by itself
orient δ1 versus δ2 (the constraint admits both orderings); it collapses
*cross-residue* ambiguity by coupling the pair's NOE evidence, and hands the
orientation to the score. Its net effect on accuracy is target-dependent: it helps
MBP (87.5% → 93.2%; Table 1, +HMBC) but on the crowded REC2/REC3 domains the extra
hard links reshape the constrained landscape into a different equal-scoring optimum
that scores lower, so it too is opt-in. Running the same scoring inside the MAUS
bounds thus yields roughly an order of magnitude more single-answer accuracy than
scoring over the full space, at no cost to the certainty guarantee, and improves
with added experimental input where that input is informative.

On the real-experimental TNF-α homotrimer (matched at H±0.01/C±0.05 with reciprocal
symmetric NOESY pairing) the method behaves consistently with the
simulated targets while exposing two properties only genuine data reveals
(Table 1). First, multimer handling is load-bearing: with the same reciprocal
symmetric NOE edges, parsing the structure as a single protomer leaves 10 measured
inter-subunit contacts without a structural home and drops them from the envelope
(75/85 = 88.2%); retaining the three chains as symmetry images (contact = minimum
distance over subunits) explains 9 of them and recovers a 98.8% envelope (84/85).
Second, the never-exclude guarantee is conditional on the NOEs being consistent with
the structure — a condition simulated data satisfies by construction. Measured
against an AlphaFold3 model, one peak (Leu76δ2) carries a symmetric-confirmed NOE the
predicted structure does not support at the 6/10 Å cutoffs and falls out of the
envelope, so the real-world figure is 98.8%, not 100%. (Reciprocal symmetric pairing
is what makes even this bound attainable: the plain carbon-only firm match — resolving
the partner by carbon alone — collects mutually inconsistent hard edges on this dense
trimer and is UNSAT outright, so it commits nothing; pairing each row with its
reciprocal resolves both ends by full (C,H) and every retained edge is a real
methyl–methyl contact.) Committed accuracy is 37.6% methyl (35.3% before soft
evidence) and 55.3% at residue level — short of the simulated targets, as expected for
a boolean-dominated (H)CCH network scored against a predicted fold. The methyl-vs-
residue gap is entirely geminal swaps: right residue, wrong δ1/δ2 or γ1/γ2 (Table 3).
Because the scoring objective's global optimum is *not* the truth here (climbing it —
more annealing restarts, or a normalized objective — lowers accuracy, since real
intensities against a predicted fold do not pin the true distances), the swaps are
broken instead by a deterministic geminal resolver: for each Leu/Val pair whose two
methyls are both committed, the peak with the stronger NOE cross peak to a shared
partner Q takes the methyl structurally closer to Q (I ~ 1/r⁶), a local 2-way flip
that raises Leu/Val methyl accuracy from 22.6% to 27.4% (Leu 25→33%) and overall
methyl from 34.1% to 37.6% without disturbing the residue assignment or the envelope.
The optional HMBC lever does *not* help this run — only one geminal link matches the
peaks at H±0.01/C±0.05, and adding it to the already max-feasible commit set tips the
SAT infeasible rather than orienting the pairs (Table 1, n.r.).

**Table 1.** Methyl-level accuracy on seven structure-simulated targets plus one
real-experimental multimer (TNF-α), all engines scored on the same 1/r^6^
intensity NOESY. Envelope = fraction of peaks whose MAUS option set contains the
truth (never-exclude guarantee). n.c. = did not converge within a 15-min budget;
n.r. = not reported (the single matched HMBC link over-constrains TNF-α's
max-feasible commit set, tipping the SAT infeasible).

| Target | BMRB / PDB | Labeling | Methyls | MAGIC | MAUS | magicmaus | +soft | +HMBC | Envelope |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Ubiquitin | 6457 / 1UBQ | AILMTV | 43 | 9.3% | 34.9% | 100% | 100.0% | 100.0% | 100% |
| IL-2 | 28104 / 1M47 | ILV | 59 | 8.5% | 8.5% | 88.1% | 96.6% | 96.6% | 100% |
| HNH (Cas9) | 27949 / 6O56 | AILTV | 57 | 12.3% | 26.3% | 82.5% | 86.0% | 86.0% | 100% |
| REC2 (Cas9) | 28105 / 4CMP | ILV | 63 | n.c. | 12.7% | 88.9% | 90.5% | 76.2% | 100% |
| REC3 (Cas9) | 28110 / 4ZT0 | ILV | 85 | n.c. | 8.2% | 60.0% | 57.6% | 52.9% | 100% |
| MBP | 7114 / 1ANF | AILMTV | 192 | 5.7% | 26.6% | 87.0% | 87.5% | 93.2% | 100% |
| MSG | SI† / 1D8C | ILV | 257 | n.c. | 1.6% | 29.6% | 33.5% | 38.5% | 100% |
| TNF-α‡ | real / AF3 trimer | AILTV | 85 | 2.4% | 0.0% | 35.3% | 37.6% | n.r. | 98.8% |

†MSG methyl shifts have no BMRB deposit; they are digitised from the reference
assignment table in the open-access Supplementary Information of Pritišanac *et
al.* (2019). ‡TNF-α is the only real-experimental and only multimeric target:
genuine methyl-NMR HMQC/NOESY/HMBC peak lists of the tumour-necrosis-factor-α
homotrimer, scored against an AlphaFold3 trimer model, and matched at H±0.01/C±0.05
like the rest. Its NOESY match uses **reciprocal symmetric pairing**: a 3D (H)CCH row
gives the partner only by carbon, but the reciprocal row supplies the partner's
proton, so pairing the two resolves both NOE ends by full (C,H) and every firm edge
is a correct methyl–methyl contact — the envelope reaches 98.8% with no wrong hard
constraint (the plain carbon-only firm match is UNSAT on this dense trimer). Multimer
handling is load-bearing: treating the three chains as symmetry images (a contact is
the minimum distance over subunits) explains the inter-subunit NOEs; with the same
symmetric edges the trimer parse recovers 84/85 vs 75/85 single-protomer, and the one
remaining peak (Leu76δ2) carries an NOE the *predicted* structure does not support
(the conditional-envelope caveat below). The commit engine grows a max-feasible hard
set — the symmetric seed plus the carbon-only firm edges that keep the SAT feasible —
so it prunes enough to commit (its own envelope drops, but the reported Envelope column
is engE's 98.8%). The MAUS column is the fraction of peaks committed uniquely (all
correct); none is uniquely forced here, so its coverage is entirely the Envelope
column. The +HMBC column would add the optional HMBC-HMQC geminal-link experiment
(`--hmbc`) on top of +soft; only one geminal link matches at this tolerance and adding
it to the already max-feasible commit set tips the SAT infeasible, so it is not reported
(n.r.).

The labeling scheme (Table 1) shapes the assignment problem in two ways that the
residue-resolved accuracy of Table 2 makes explicit. First, because a peak competes
only with methyls of its own residue type, the number of labeled types sets the
granularity of the candidate partition: the two AILMTV targets (ubiquitin, MBP)
spread their peaks over six types, so each peak competes against fewer same-type
methyls than in a three-type ILV set. Second, labeling sets the irreducible geminal
load: Leu (δ1/δ2) and Val (γ1/γ2) carry prochiral methyl pairs that an achiral NOE
network cannot orient, whereas Ile contributes a single δ1 methyl. Table 2 shows the
two effects: Ile is the most reliably assigned type — 100% on five of seven targets
— because it is single-methyl and, being typically buried, NOE-rich; Leu and Val
are the bottleneck, at 100% only where the network is dense (ubiquitin, IL-2) and
falling to ~30–50% on the sparse REC3/MSG networks. This is why the residual
`ambiguous` tier is dominated by Leu/Val pairs — an orientation the intensity score
resolves only when the two methyls make sufficiently different structural contacts,
and reports as a coin flip when they do not.

Being single-methyl, however, is necessary but not sufficient: removing geminal
ambiguity and carrying NOE information are independent. On MBP, Ala Cβ (68%) and Thr
Cγ2 (70%) are the *least* accurately assigned types (Table 2) despite no prochiral
degeneracy, because these often surface-exposed, shift-clustered methyls simply make
too few methyl–methyl NOEs to be pinned — whereas Ile on the same protein is 100%.
Labeling is therefore not the sole determinant of difficulty; it interacts with NOE
information content. IL-2 and REC2 are Leu-rich (71% and 76% Leu) yet reach ~88–90%
because their networks are dense, while REC3 and MSG fall to 60% and 30% not from
their labeling alone but because sparse carbon-only matching leaves too few firm NOEs
to pin the enlarged ILV domains (Section 3, coverage gate). In practice the levers
attack different bottlenecks: extending the labeling beyond ILV partitions the
problem into finer type classes, an HMBC geminal-link experiment collapses the
cross-residue ambiguity by tying each geminal pair to one residue, and the
intensity-weighted score is what finally orients δ1/δ2 within the pair — none is a
universal fix, and the last two are opt-in.

**Table 2.** Methyl-level accuracy of the committed magicmaus call (+soft) resolved
by residue type. Each cell is correct/observed for that type; a dash marks a type
absent from the target's labeling. Ile is near-perfect except on the sparse
REC3/MSG networks; Leu/Val carry the geminal degeneracy; Ala/Thr, though
single-methyl, are NOE-poor on MBP.

| Target | Ile | Leu | Val | Ala | Thr | Met | Total |
|---|---|---|---|---|---|---|---|
| Ubiquitin | 100% (7/7) | 100% (18/18) | 100% (8/8) | 100% (2/2) | 100% (7/7) | 100% (1/1) | 100% (43/43) |
| IL-2 | 100% (9/9) | 100% (42/42) | 75% (6/8) | — | — | — | 97% (57/59) |
| HNH | 100% (7/7) | 86% (24/28) | 75% (12/16) | 100% (2/2) | 100% (4/4) | — | 86% (49/57) |
| REC2 | 100% (9/9) | 88% (42/48) | 100% (6/6) | — | — | — | 90% (57/63) |
| REC3 | 54% (7/13) | 52% (26/50) | 73% (16/22) | — | — | — | 58% (49/85) |
| MBP | 100% (22/22) | 97% (58/60) | 95% (38/40) | 68% (30/44) | 70% (14/20) | 100% (6/6) | 88% (168/192) |
| MSG | 44% (18/41) | 32% (42/133) | 31% (26/83) | — | — | — | 33% (86/257) |
| TNF-α | 100% (8/8) | 33% (12/36) | 19% (5/26) | 50% (6/12) | 33% (1/3) | — | 38% (32/85) |

The real-data TNF-α row (symmetric-NOESY engine at H±0.01/C±0.05, with the geminal
intensity-ratio resolver) makes the type dependence stark even as overall accuracy is
modest on this sparse network: Ile (no geminal partner) is assigned outright at 100%
(8/8), while Leu (33%) and Val (19%) lag — the firm NOEs cannot fully orient the
geminal pairs. The geminal resolver (I ~ 1/r⁶: the stronger NOE to a shared partner
belongs to the closer methyl) lifts Leu from 25 to 33% and Leu/Val together from 22.6
to 27.4%; the residual is where no shared firm partner separates δ1/δ2 (or γ1/γ2),
which stays inside the MAUS envelope as a genuine coin flip. The residue-level accuracy
is 55% (47/85): the gap to the 38% methyl-level figure is entirely these unresolved
geminal swaps — right residue, wrong δ1/δ2 or γ1/γ2.

Resolving the accuracy to the individual prochiral methyl (Table 3) confirms that
the residual Leu/Val error is a geminal swap rather than random misassignment: on
every target the two members of each pair degrade in near-lockstep — HNH Leu δ1 and
δ2 are both 86%, its Val γ1 and γ2 both 75%, REC2 Leu δ1/δ2 both 88% — the exact
signature of the achiral network placing the pair on the right *residue* but the
wrong *methyl*. Both members stay inside the MAUS envelope (the truth is never
excluded), so the swap is a calibrated coin flip that the confidence tier flags as
`ambiguous` and that only a signal able to tell δ1 from δ2 can break — the intensity
score where their structural contacts differ, or an independent stereospecific
assignment where they do not. An HMBC geminal-link experiment does *not* break it:
by construction its constraint admits both orderings of the pair (Section 3), so it
resolves which residue, not which methyl. Ile δ1, having no geminal partner, carries
no such symmetry and is assigned outright wherever the network is dense.

**Table 3.** Accuracy resolved to the individual methyl carbon (magicmaus +soft).
Geminal partners (Leu δ1/δ2, Val γ1/γ2) are listed separately; their near-equal
columns are the geminal-swap signature. A dash marks a methyl absent from the
target's labeling.

| Target | Ile δ1 | Leu δ1 | Leu δ2 | Val γ1 | Val γ2 | Ala β | Thr γ2 | Met ε |
|---|---|---|---|---|---|---|---|---|
| Ubiquitin | 100% (7/7) | 100% (9/9) | 100% (9/9) | 100% (4/4) | 100% (4/4) | 100% (2/2) | 100% (7/7) | 100% (1/1) |
| IL-2 | 100% (9/9) | 100% (21/21) | 100% (21/21) | 75% (3/4) | 75% (3/4) | — | — | — |
| HNH | 100% (7/7) | 86% (12/14) | 86% (12/14) | 75% (6/8) | 75% (6/8) | 100% (2/2) | 100% (4/4) | — |
| REC2 | 100% (9/9) | 88% (21/24) | 88% (21/24) | 100% (3/3) | 100% (3/3) | — | — | — |
| REC3 | 54% (7/13) | 52% (13/25) | 52% (13/25) | 73% (8/11) | 73% (8/11) | — | — | — |
| MBP | 100% (22/22) | 97% (29/30) | 97% (29/30) | 95% (19/20) | 95% (19/20) | 68% (30/44) | 70% (14/20) | 100% (6/6) |
| MSG | 44% (18/41) | 33% (22/67) | 30% (20/66) | 31% (13/42) | 32% (13/41) | — | — | — |
| TNF-α | 100% (8/8) | 28% (5/18) | 39% (7/18) | 15% (2/13) | 23% (3/13) | 50% (6/12) | 33% (1/3) | — |

The same geminal-swap lens revises how the MAUS column should be read. MAUS's
decisive fraction (Table 1, MAUS) counts only peaks pinned to a single *methyl*, so
it books every geminal-unresolved pair as an abstention — yet such a pair is already
decisive at the *residue* level: the truth's residue is fixed, only its δ1/δ2 (or
γ1/γ2) label is open. Collapsing each option set to its residues, MAUS is in fact
residue-decisive — and, by never-exclude, correct — on far more peaks than its
methyl-unique count suggests: 88.4% vs 34.9% on ubiquitin, 53.6% vs 26.6% on MBP,
and 50.9% vs 26.3% on HNH. This residue-vs-methyl gap is a direct readout of how much
of MAUS's ambiguity is *merely geminal* — a δ1/δ2 orientation that only the score (or
a stereospecific measurement) can settle, since it is the one thing the HMBC link
leaves open. Its complement, the fraction still spanning more than one residue, is the
*cross-residue* ambiguity, and it dominates the Leu-crowded ILV targets (IL-2 81%, REC2
68%, MSG 94% of peaks) where many Leu compete for one peak. That is exactly the part
the HMBC link collapses: tying each geminal pair to one residue couples their NOE
evidence and lifts the residue-decisive fraction sharply on those targets (REC2 32% →
68%, HNH 51% → 79%). Whether that structural collapse improves the final methyl-level
call, however, is target-dependent — it does on MBP but reshapes the crowded REC2/REC3
landscapes into a worse-scoring optimum — which, together with the fact that the swap
itself is left to the score, is why HMBC is opt-in rather than a universal lever.

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

## Acknowledgements

The authors gratefully acknowledge Eun-Hee Kim and Dr. Hae-Kap Cheong from the
Korea Basic Science Institute (KBSI) for their valuable support and discussions.

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
to a per-peak option set that provably contains the truth (100% envelope) and
prunes candidates from up to ~60 to 1–3; an intensity-weighted NOE score
(MAGIC-style), applied only within those bounds via a SAT-feasible seed and
feasibility-preserving 3-cycle simulated annealing, commits to a single coherent map with a
per-peak confidence tier (unique / scored / ambiguous). (**B**) Methyl-level
accuracy across the seven structure-simulated benchmark targets (43–257 methyls,
ordered by size) plus the real-experimental TNF-α homotrimer (*real),
all engines scored on the same 1/r^6^ intensity NOESY: magicmaus (+soft, blue)
versus MAUS unique-only calls (purple; the constraint layer's decisive fraction,
abstaining elsewhere) and full-space MAGIC (red; hatched *n.c.* where scoring did
not converge within a 15-min budget). Green markers denote the truth-in-envelope
guarantee — 100% on every simulated target, and 98.8% on real TNF-α data via
reciprocal symmetric NOESY pairing (one peak carries an NOE the predicted structure
cannot support; the plain carbon-only match is UNSAT on the dense trimer).](figure1.png)
