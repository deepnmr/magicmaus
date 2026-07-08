# TNF-α homotrimer — magicmaus example (real experimental data)

The **only real-experimental** target in the benchmark, and the only
**multimer**. Everything else in `examples/` is a structure-simulated NOESY;
this one is genuine methyl-NMR peak lists of the tumour-necrosis-factor-α
homotrimer, so it is where the method meets real spectral noise, real
inter-subunit NOEs, and a *predicted* (not crystallographic) structure.

| source | file |
|---|---|
| ILVAT ¹³C-HMQC peak list (Sparky) | `TNFa_ILVAT_13C_HMQC.list` |
| 3D `(H)CCH` NOESY-HMQC peak list | `TNFa_ILVAT_NOESY_HMQC.list` |
| 3D HMBC-HMQC peak list | `TNFa_ILVAT_HMBC_HMQC.list` |
| structure | `fold_tnfa_trimer_model_0.cif` (AlphaFold3 trimer model, chains 1/2/3) |

The three `.list` files are the raw assigned spectra. `make_tnfa_input.py`
(repo root) strips the assignments into the magicmaus input set and keeps the
answer separately:

| generated file | contents |
|---|---|
| `hmqc.tsv` | 85 methyl HMQC **input** peaks: `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type` (anonymous `P1…P85`) |
| `hmqc_true.tsv` | truth key: `… ⇥ True` (e.g. `V1CG1`), scoring only |
| `noesy.tsv` / `noesy_intensity.tsv` | 300 NOESY cross peaks, boolean and with Sparky peak heights as intensity |
| `hmbc.tsv` | 69 HMBC geminal-link rows |

3 HMQC rows are dropped as un-typable (`?-?` ×2, the ambiguous `vat85C-H`).
`V1` is present in the AlphaFold model, so all 85 assigned peaks have a
structural home.

Regenerate the inputs (they are committed):

```bash
python make_tnfa_input.py examples/TNFa
```

## Multimer is load-bearing

TNF-α is a symmetric homotrimer; `maus.parse_structure` keeps all three chains
as symmetry images per methyl and scores every contact by the **minimum
distance over subunits**, so inter-subunit NOEs count. This is not optional
here — several real NOEs cross the subunit interface:

| structure model (symmetric NOE edges) | truth-in-envelope |
|---|---|
| single protomer (chain 1 only) | 75/85 = 88.2% |
| full trimer (all 3 chains) | **84/85 = 98.8%** |

Parsed as a monomer, 10 real inter-subunit NOEs have no structural explanation
and force those peaks out of the envelope; the trimer parse explains 9 of them and
recovers a 98.8% envelope (the one remaining peak, Leu76δ2, carries an NOE the
predicted structure does not support).

## Run

All runs use the wider **H±0.02 / C±0.1 ppm** match tolerance appropriate to the
real, broader linewidths (the simulated targets use ±0.01/±0.05).

```bash
python magicmaus.py examples/TNFa/fold_tnfa_trimer_model_0.cif \
    examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv \
    --truth examples/TNFa/hmqc_true.tsv --tol-h 0.02 --tol-c 0.1
```

```
methyls(G nodes)=89  HMQC peaks=85  NOESY cross peaks=300
NOE match (tol H±0.02/C±0.1): firm=47 ambiguous(dropped)=229 unmatched=24
truth in MAUS option set = 84/85 = 98.8%
magicmaus single call    = 6/85 = 7.1% correct   (scored 4/28, ambiguous 2/57)
```

`--soft-ambiguous` gives 8/85 = 9.4%; adding `--hmbc` on top lifts it to
9/85 = 10.6% (at the cost of a slightly tighter 96.5% envelope). This plain
carbon-only CLI stays feasible only at the wider ±0.02/±0.1 tolerance; at the
±0.01/±0.05 tolerance the carbon-only firm edges collide and go UNSAT — use the
symmetric run below for the tight-tolerance result. Or the same three-engine
benchmark as the other targets:

```bash
python bench.py examples/TNFa examples/TNFa/fold_tnfa_trimer_model_0.cif
```

## NOESY symmetry — a 98.8% envelope at tight tolerance

A 3D (H)CCH NOESY row `(C1_partner, C2_obs, H2_obs)` gives the partner **only by
carbon**. At the tight ±0.01/±0.05 tolerance the carbon-only firm edges on this
dense trimer collide into mutually inconsistent hard constraints and the SAT goes
**UNSAT** (it commits nothing). The reciprocal row `(C2_partner, C1_obs, H1_obs)`
supplies the partner's **proton**: pairing the two resolves both endpoints by full
`(C,H)`, so every retained edge is a correct methyl–methyl contact.
`run_tnfa_symmetric.py` builds those edges for the **envelope**, then grows a
max-feasible hard set (symmetric seed + the carbon-only firm edges that keep the SAT
feasible) for the **commitment**, with the carbon-only ambiguous rows as soft
evidence, and finally breaks the Leu/Val geminal swaps with an intensity-ratio
resolver:

```bash
python run_tnfa_symmetric.py
```

```
symmetric NOE edges = 76   HMBC gem-links = 1   geminal flips = 5
MAUS envelope (symmetric)   = 84/85 = 98.8%  (never excludes truth)
committed (greedy+geminal)  = methyl 32/85 = 37.6%  residue 47/85 = 55.3%
committed call in envelope  = 84/85 = 98.8%
```

`run_tnfa_symmetric.py` **merges the two ends of the trade-off** into one
assignment (written to `magicmaus_calls_symmetric.tsv`): the symmetric edges give
every peak a **98.8% never-exclude envelope** (only Leu76δ2, on a
structure-unsupported NOE, is excluded), while the greedy max-feasible hard set plus
the geminal resolver supplies the **best single committed call** (37.6% methyl / 55%
residue). Nearly every committed call falls inside the guaranteed envelope, so each
peak reports a best-guess assignment nested inside a bound that (bar the one peak)
contains the truth. This is the strongest combined result on the given peak lists.

**Geminal intensity-ratio resolver.** The scoring objective's global optimum is *not*
the truth on real intensities against a predicted fold — climbing it (more annealing
restarts, or a per-peak-normalized objective) only *lowers* accuracy. So the Leu/Val
δ1↔δ2 (γ1↔γ2) swaps are broken by a deterministic local rule instead: for a geminal
pair whose two methyls are both committed to peaks Pa,Pb, a shared NOE partner Q gives
a vote `sign(I(Pa,Q) − I(Pb,Q)) · sign(1/r⁶(G1,Q) − 1/r⁶(G2,Q))` — the stronger cross
peak belongs to the methyl closer to Q — and the pair is flipped if the firm-edge
majority disagrees with the current orientation. Five flips lift Leu/Val methyl
accuracy from 22.6 % to 27.4 % (Leu 25 → 33 %) and overall methyl from 34.1 % to
37.6 %, without touching the residue assignment or the envelope. Where no shared firm
partner separates the pair, the swap stays a genuine coin flip inside the envelope.

## MAGIC (sibling engine)

MAGIC is now multimer-aware too: `magic/structure.py` keeps every chain as a
symmetry image (keyed by residue, no chain id — same collapse as
`maus.parse_structure`) and scores each contact by the **minimum distance over
subunits**, so inter-subunit NOEs enter its objective. `convert_to_magic.py`
renders the trimer directly (a `.cif` is converted to a multi-chain PDB), so no
protomer step is needed; `../magic/` scores it (committed result:
`magic_assignments.tsv`):

```bash
python convert_to_magic.py examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv \
    examples/TNFa/fold_tnfa_trimer_model_0.cif magic_run_tnfa/ --tol-h 0.02 --tol-c 0.1
# then from ../magic/:
python -m magic run ../magicmaus/magic_run_tnfa/control.txt --output-dir out
```

The trimer parse adds 64 inter-subunit structural contacts (model-matrix sum
642.8 → 715.9), lifting the best objective from 89.7 (protomer) to 156.4. Even
so MAGIC assigns only **2/85 = 2.4%** methyls (residue-level 10.6%): its
full-space scoring sits on a near-flat real-data landscape, so a
physically-correcter objective does not track the truth — the same regime it
hits on the simulated targets, and exactly why magicmaus scores *inside* the
MAUS envelope instead. (The protomer run scored 5/85 methyl / 7/85 residue; the
multimer objective moves calls around without net gain — climbing the correct
landscape is not enough without the MAUS bound.)

## From raw spectra — the fully-automatic pipeline

Everything above starts from the *assigned* Sparky `.list` files. A separate,
harder track starts from the raw **`.ucsf` spectra** (not committed — ~1.3 GB;
`*.ucsf` is git-ignored) and does peak-picking + typing + assignment with no
human input, via `ucsf.py` (a numpy-only UCSF reader), `make_tnfa_peaks.py`, and
`run_tnfa_picked.py`:

```bash
python make_tnfa_peaks.py     # pick + type -> hmqc_picked.tsv, noesy_picked.tsv
python run_tnfa_picked.py     # magicmaus + score vs the known answer by ppm
```

**Typing by label-selective spectra.** The `ILVAT` HMQC is the tag-free master
list; the type-selective HMQCs assign a residue type to each peak and remove tag
peaks (present in `Val`/`Thr` samples, which carry an N-terminal tag, but absent
from the tag-free `ILVAT`):

- `I` = ¹³C **< 17 ppm** — Ile δ1 is the only methyl that low (Ile 12.8–15.8;
  every other type ≥ 18.3), a clean shift separator that the Thr sample's ¹³C
  window would otherwise clip
- **geminal link (3D HMBC)** — the 3D HMBC-HMQC (`C, C, H`) carries both methyl
  carbons and the proton, so the observed methyl is pinned unambiguously by
  `(C,H)` and the other carbon is its geminal partner (far cleaner than a 2D
  `(C_partner, H)` list). A peak with a geminal partner is Leu or Val; a single
  methyl (no partner) is Ile/Ala/Thr. This separates the paired from the
  single-methyl types and lets `V` propagate across the link: if a peak's geminal
  partner is in `Val_Methyl`, so is it (recovers Val the low-SNR `Val_Methyl`
  detects on only one of the pair).
- `V` = in `Val_Methyl` (directly or via its geminal partner)
- `L` = geminal, not Val · `T` = single, in `Thr_Methyl` · `A` = single, in
  `ILVAT` only

Each type is then capped to its structural methyl count (e.g. ≤36 Leu), because
the injective SAT is infeasible — and the C solver hangs (~30 min on the
over-capacity instance) — if a type has more peaks than methyls.

**Result (committed answer scored against the known assignment by ppm):**

| stage | outcome |
|---|---|
| ILVAT peaks picked → typed (capped) | 87 peaks; **80/85** true peaks recovered |
| type correct | **73/80 = 91%** (Ile 8/8 by shift; 3D-HMBC geminal fills Val to 26/26) |
| NOESY cross peaks | **symmetric** filter: 226 of 315 (37 resolve to firm constraints) |
| NOESY intensity | box-integrated **volume** (±3 ¹³C, ±1 ¹H points) |
| **MAUS envelope** | **74/80 = 92.5%** truth-in-option-set |
| **magicmaus committed** | **8.8%** methyl-level (11.2% residue-level) |

The 3D HMBC geminal link fills the Val type to its full 26/26 (each Val whose
partner the low-SNR `Val_Methyl` catches is propagated), lifting the envelope
from 87.5% (2D HMBC) to **92.5%** — near the curated-`.list` symmetric figure of 98.8%.

**Symmetric NOESY.** A real methyl–methyl NOE appears both ways — (C1,C2,·) and
(C2,C1,·) — while one-sided picking noise does not. Keeping only symmetric cross
peaks (226 of 315) removes the false hard constraints that made the SAT **UNSAT**
above ~15 firm edges, so all symmetric edges are used (37 firm) with no arbitrary
top-N cap and the envelope intact.

Shift-based Ile typing plus 3D-HMBC geminal linking lifts the type accuracy to
91% and the MAUS envelope to 92.5% (from 77.5% with the first, cruder typing).
The symmetric NOESY filter removes the SAT fragility (37 firm edges now feasible,
no arbitrary cap), and box-integrated **volumes** (better NOE-intensity estimates
than single-point heights) lift the committed call from ~6% to 8.8% methyl /
11.2% residue. It is still modest: the auto-picked, boolean-ish (H)CCH network is
sparse and the volumes are noisy, so the MAGIC-style scoring can rank the right
methyl within the (correct, 92.5%) option set only sometimes. Lowering the pick
threshold does **not** help — below it the extra "cross peaks" are noise that
pairs up symmetrically by density alone, so firm edges balloon (37 → 82 → 212)
and the SAT goes UNSAT again; the genuine NOEs are already all captured at the
current threshold. This is the honest limit of end-to-end automation on this
dataset — a strong bounded envelope, weak but non-trivial commitment — versus the
assigned-`.list` track above, where curated peaks let the scoring layer commit
outright.

## What the numbers mean

- **The 100% guarantee is conditional.** MAUS never excludes the truth *when the
  NOEs are consistent with the structure*. On real data measured against a
  predicted structure, one peak (Leu76δ2) carries a symmetric-confirmed NOE the
  AlphaFold model does not support at the 6/10 Å cutoffs, so it falls out of the
  envelope. 98.8%, not 100%, is the real-world number.
- **A sparse network caps the committed accuracy.** Only 47 of 300 cross peaks
  resolve to a firm constraint at the wide tolerance (and the tight-tolerance
  carbon-only match is UNSAT), so most peaks land in large option sets and the
  scored/ambiguous calls are near coin flips. The symmetric run's max-feasible
  commit reaches 34% methyl / 55% residue; the residual is honestly flagged
  `ambiguous`, not guessed.
- **HMBC does not help this run.** `--hmbc` turns each matched HMBC-HMQC cross peak
  into a hard geminal link. At the wide ±0.02/±0.1 tolerance it nudges the
  carbon-only commit up (9.4 → 10.6%) at a slightly tighter 96.5% envelope, but at
  the tight ±0.01/±0.05 symmetric setup only one link matches and adding it to the
  already max-feasible commit set tips the SAT infeasible — so it is not a net lever
  on this dataset.
