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
| `noesy.tsv` / `noesy_intensity.tsv` | 220 NOESY cross peaks, boolean and with Sparky peak heights as intensity |
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

| structure model | truth-in-envelope |
|---|---|
| single protomer (chain 1 only) | 75/85 = 88.2% |
| full trimer (all 3 chains) | **81/85 = 95.3%** |

Parsed as a monomer, 6 real inter-subunit NOEs have no structural explanation
and force those peaks out of the envelope; the trimer parse explains them and
recovers a coherent 95.3% envelope.

## Run

All runs use the wider **H±0.02 / C±0.1 ppm** match tolerance appropriate to the
real, broader linewidths (the simulated targets use ±0.01/±0.05).

```bash
python magicmaus.py examples/TNFa/fold_tnfa_trimer_model_0.cif \
    examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv \
    --truth examples/TNFa/hmqc_true.tsv --tol-h 0.02 --tol-c 0.1
```

```
methyls(G nodes)=89  HMQC peaks=85  NOESY cross peaks=220
NOE match (tol H±0.02/C±0.1): firm=37 ambiguous(dropped)=172 unmatched=11
truth in MAUS option set = 81/85 = 95.3%
magicmaus single call    = 12/85 = 14.1% correct   (scored 8/26, ambiguous 4/59)
```

`--soft-ambiguous` gives 10/85 = 11.8%; adding `--hmbc` on top lifts it to
15/85 = 17.6% (the one experimental lever that helps here). Or the same
three-engine benchmark as the other targets:

```bash
python bench.py examples/TNFa examples/TNFa/fold_tnfa_trimer_model_0.cif
```

## MAGIC (sibling engine)

MAGIC has no multimer support, so it is run on a single protomer. `tnfa_protomer.pdb`
is chain 1 of the AlphaFold model in PDB format; `convert_to_magic.py` builds the
control bundle and `../magic/` scores it (committed result: `magic_assignments.tsv`):

```bash
python convert_to_magic.py examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv \
    examples/TNFa/tnfa_protomer.pdb magic_run_tnfa/ --tol-h 0.02 --tol-c 0.1
# edit LABELING to AILTV (no Met), then from ../magic/:
python -m magic run ../magicmaus/magic_run_tnfa/control.txt --output-dir out
```

MAGIC assigns **2/85 = 2.4%** methyls (residue-level 5.9%) — full-space scoring
over a near-flat real-data landscape, the same regime it hits on the simulated
targets.

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
| NOESY cross peaks | top 80 by height (~13 resolve to firm constraints) |
| **MAUS envelope** | **74/80 = 92.5%** truth-in-option-set |
| **magicmaus committed** | single-digit % (NOESY-noise-limited, see below) |

The 3D HMBC geminal link fills the Val type to its full 26/26 (each Val whose
partner the low-SNR `Val_Methyl` catches is propagated), lifting the envelope
from 87.5% (2D HMBC) to **92.5%** — near the curated-`.list` figure of 95.3%.

Shift-based Ile typing plus 3D-HMBC geminal linking lifts the type accuracy to
91% and the MAUS envelope to 92.5% (from 77.5% with the first, cruder typing).
The committed call stays single-digit, capped by the NOESY: an auto-picked,
boolean-ish (H)CCH network with raw peak heights gives the scoring layer little
to grade, and picking noise makes the SAT **UNSAT** once ~16 cross peaks resolve
to firm hard constraints (so the NOESY is held to its top 80). This is the honest
limit of end-to-end automation on this dataset — a strong bounded envelope but
weak commitment — versus the assigned-`.list` track above, where curated peaks
let the scoring layer commit.

## What the numbers mean

- **The 100% guarantee is conditional.** MAUS never excludes the truth *when the
  NOEs are consistent with the structure*. On real data measured against a
  predicted structure, 4/85 peaks (V13γ1, A22β, L29δ1, L94δ1) carry NOEs the
  AlphaFold model does not support at the 6/10 Å cutoffs, so they fall out of the
  envelope. 95.3%, not 100%, is the real-world number.
- **A sparse network caps the committed accuracy.** Only 37 of 220 cross peaks
  resolve to a firm constraint at this tolerance, so most peaks land in large
  option sets and the scored/ambiguous calls are near coin flips — 14.1% committed
  (11.8% +soft). The residual is honestly flagged `ambiguous`, not guessed.
- **HMBC helps here.** `--hmbc` turns each matched HMBC cross peak into a hard
  geminal link. At this wider tolerance only 4 of 69 rows match firmly, but those
  matches are clean enough to couple the geminal pairs' NOE evidence, lifting the
  committed call 11.8 → 17.6% — the one experimental lever that helps on this
  dataset. (At a tighter tolerance the HMBC matches include a wrong geminal pair
  and one bad hard link can render the SAT infeasible, so `--hmbc` is
  tolerance-sensitive, as elsewhere.)
