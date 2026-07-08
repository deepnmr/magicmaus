# magicmaus (`mm/`) — fresh MAUS + MAGIC fusion

A clean-room build that fuses two methyl-NMR assignment ideas. It does **not**
import the repo's older `maus.py` / `magicmaus.py`; it is self-contained.

- **MAUS idea** — cast assignment as subgraph isomorphism, solve as SAT, and
  return for every HMQC peak the *set* of methyls consistent with all hard
  constraints. Provably never excludes the truth; abstains under degeneracy.
- **MAGIC idea** — score a global peak→methyl map with a distance-weighted NOE
  objective (intensity ∝ 1/r⁶) and commit to the best scorer.

magicmaus uses MAUS to bound the space with certainty, then MAGIC-style scoring
to commit one coherent call **within** the residual degeneracy — so the answer
keeps MAUS's never-exclude envelope while delivering MAGIC's single call.

## Inputs

| Experiment | File | Columns |
|---|---|---|
| HMQC (2D) | `hmqc.tsv` | `label  H_ppm  C_ppm  res_type` |
| 3D HMBC-HMQC | `hmbc.tsv` | `label  C2  C1  H1` |
| 3D NOESY | `noesy.tsv` | `label  C2  C1  H1  [intensity]` |
| Structure | `.cif` / `.pdb` | mmCIF preferred; homo-oligomer → multimer |

In the 3D lists the **detected** methyl resonates at `(C1, H1)` (its own carbon
and proton); the **partner** methyl contributes carbon `C2` only — the
irreducible ambiguity source. Structure is mmCIF (or PDB); a homo-oligomer's
subunits become chain images and inter-subunit NOE contacts are the min distance
over images (**multimer-aware**).

## Usage

```bash
# assign from peak lists (+ optional truth key for scoring)
python -m mm assign STRUCTURE.cif hmqc.tsv noesy.tsv --hmbc hmbc.tsv --truth hmqc_true.tsv

# simulate noisy peak lists from a structure + BMRB shifts
#   measurement error: H ~ N(0, 0.02 ppm), C ~ N(0, 0.1 ppm), one draw per methyl
python -m mm simulate STRUCTURE.cif --bmrb bmrXXXX.str --out-dir out/ --seed 1
python -m mm simulate STRUCTURE.cif --shifts-tsv hmqc_true.tsv --out-dir out/   # or a shift table
```

Library entry points: `mm.assign(...)` and `mm.simulate.simulate(...)`.

## Worked examples (this repo)

```bash
# BMRB monomer (ubiquitin): simulate → assign
python -m mm simulate examples/ubq/1UBQ.pdb --bmrb examples/ubq/bmr6457.str --out-dir /tmp/ubq --seed 1
python -m mm assign  examples/ubq/1UBQ.pdb /tmp/ubq/hmqc.tsv /tmp/ubq/noesy.tsv --hmbc /tmp/ubq/hmbc.tsv --truth /tmp/ubq/hmqc_true.tsv
# -> truth in MAUS option set = 43/43 (100%); single call ~86%

# multimer (TNF-α homotrimer, AlphaFold cif): inter-subunit NOEs over 3 chains
python -m mm simulate examples/TNFa/fold_tnfa_trimer_model_0.cif --shifts-tsv examples/TNFa/hmqc_true.tsv --out-dir /tmp/tnfa --seed 2
python -m mm assign  examples/TNFa/fold_tnfa_trimer_model_0.cif /tmp/tnfa/hmqc.tsv /tmp/tnfa/noesy.tsv --hmbc /tmp/tnfa/hmbc.tsv --truth /tmp/tnfa/hmqc_true.tsv
# -> truth in MAUS option set = 85/85 (100%); single call under 3-fold symmetry ~68%
```

The never-exclude guarantee holds **when every firm edge is correct**. Under the
N(0,σ) measurement error two distinct methyls' shifts can collide within the
match tolerance and produce a *wrong* "firm" edge that prunes the truth; the CLI
flags this (`WARNING: N truth pruned by a wrong firm edge under noise`) instead
of silently claiming the envelope is intact. A multimer whose chains use
non-shared residue numbering also warns (inter-subunit contacts would be lost).

`python test_mm.py` runs the self-check (multimer images, reproducible noise,
never-exclude, (C2,C1,H1) matching).

## MAUS vs MAGIC vs magicmaus

`python compare_all.py` runs all examples end to end (simulate → assign) and
scores the three ideas **residue-wise** — correct residue, ignoring the geminal
CD1/CD2 (Leu) / CG1/CG2 (Val) swap, which is a near-symmetric coin flip.

- **MAUS** — hard SAT bounds only, no scoring. It cannot rank within an
  ambiguous option set, so a forced single call must tiebreak arbitrarily. Its
  real strength is the *envelope* (option sets retain the truth) and the subset
  it pins uniquely, not a single call.
- **MAGIC** — the 1/r⁶ NOE score optimised over the **un-pruned** full
  type-matched space, with no hard combinatorial constraint. Over that huge,
  near-flat landscape the multistart lands in a wrong basin.
- **magicmaus** — the same score optimised **inside MAUS's pruned option sets**,
  by a diverse-seed multistart (many independent feasible seeds, each greedy
  ascended). The pruning is what lets the search reach the objective's optimum.

Residue-wise accuracy, one seed, structures driven through the mmCIF path. The
BMRB rows use **simulated** peaks (true shifts + N(0,σ), σ_H=0.02, σ_C=0.1); the
TNF-α row uses the **real experimental** peak lists picked from its `.ucsf`
spectra, against the predicted AlphaFold homotrimer:

| example | n | multimer | data | MAUS | MAGIC | **magicmaus** | envelope |
|---|---|---|---|---|---|---|---|
| ubq  |  43 | –   | sim  | 72.1 | 79.1 | **100.0** | 100.0 |
| hnh  |  57 | yes | sim  | 84.2 | 54.4 | **86.0**  | 100.0 |
| il2  |  59 | –   | sim  | 39.0 | 54.2 | **64.4**  | 100.0 |
| mbp  | 192 | –   | sim  | 35.4 | 34.4 | **60.9**  | 100.0 |
| rec2 |  63 | –   | sim  | 54.0 |  1.6 | **71.4**  | 100.0 |
| rec3 |  85 | –   | sim  | 41.2 |  0.0 | **56.5**  | 100.0 |
| **TNFa** | 85 | yes | **real** | 10.6 | 20.0 | **30.6** | 96.5 |
| **mean** | | | | 48.1 | 34.8 | **67.1** | 99.5 |

magicmaus is best on **every** example, real data included.

On the real TNF-α data every method is far lower than on simulation: real
peak-picking scatter and artifacts, a *predicted* structure whose distances only
approximate the true contacts, and firm NOE edges that are sometimes wrong under
real noise (envelope drops to 96.5% — 3 truths pruned; no peak is uniquely pinned,
so MAUS alone is near-random at 10.6%). magicmaus still wins, because scoring
inside even a loose but truth-retaining envelope beats both a blind tiebreak
(MAUS) and an un-pruned search (MAGIC). Reproduce the whole table with
`python compare_all.py`; the single real case with
`python -m mm compare examples/TNFa/fold_tnfa_trimer_model_0.cif examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv examples/TNFa/hmqc_true.tsv --hmbc examples/TNFa/hmbc.tsv`.

**MSG (257 peaks) is the scale limit, reported honestly.** Its Leu/Val domains
are 138 methyls wide and the simulated data yields only ~85 firm NOE edges, so
(a) MAUS's exact option-set enumeration is intractable (>20 min even with a
unit-propagation pre-filter), and (b) the tractable fall-back — scoring over the
arc-consistency-pruned domains (`enumerate_options=False`) — leaves domains too
wide and the NOE network too sparse for magicmaus's advantage to show (magicmaus
14.8% vs MAGIC 30% vs MAUS 2.3%, envelope still 100%). magicmaus's edge needs
either tractable tight bounds or sufficient NOE density; MSG at this scale has
neither. This is a real property of the method, not a bug.

Single dataset: `python -m mm compare STRUCTURE hmqc.tsv noesy.tsv truth.tsv --hmbc hmbc.tsv`.

## Layout

```
mm/structure.py  mmCIF/PDB parse, multimer chain images, distance-classified graph
mm/peaks.py      HMQC + 3D (C2,C1,H1) parsing and frequency matching
mm/engine.py     MAUS SAT option-sets  +  MAGIC 1/r⁶ scored commitment
mm/simulate.py   BMRB / shift-table → noisy peak lists (N(0,σ)) + truth key
mm/__main__.py   `assign` / `simulate` CLI
```
