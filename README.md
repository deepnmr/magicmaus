# magicmaus

A hybrid methyl-NMR assignment engine that **fuses the two sibling projects**
[`../maus/`](../maus) (MAUS, Nerli et al., *Nat. Commun.* 2021) and
[`../magic/`](../magic) (MAGIC, Monneau et al., *J. Biomol. NMR* 2017) into one
program that is strictly better than either alone.

A short write-up in *Bioinformatics* Application Note form is in
[`manuscript.md`](manuscript.md) / [`magicmaus_AppNote.docx`](magicmaus_AppNote.docx)
(Fig. 1 rendered by [`make_figure.py`](make_figure.py)).

**New here?** Start with the **[8-methyl dummy tutorial](TUTORIAL.md)**
([PDF](TUTORIAL.pdf)) for the mechanics, then the **[maltose-binding protein
tutorial](TUTORIAL_MBP.md)** ([PDF](TUTORIAL_MBP.pdf)) for the real 192-methyl
run and the three-way comparison against MAGIC and MAUS. For the full
five-protein benchmark, see [**Benchmark**](#benchmark--five-real-bmrb-targets)
below and [`examples/*/README.md`](examples/).

Same problem — assign methyl ¹H/¹³C HMQC peaks to the methyls of a known
structure using NOESY contacts. Two opposite philosophies:

| | **MAUS** (SAT) | **MAGIC** (scoring) |
|---|---|---|
| output | *set* of methyls per peak | single best methyl per peak |
| guarantee | never excludes the truth | none |
| weakness | abstains on all degeneracy | commits to wrong answers over a near-flat landscape |

Neither dominates: MAUS **bounds the space with certainty but abstains**; MAGIC
**commits but over too large a space**. `../maus/COMPARISON.md` already spelled
out the synthesis — *use MAUS to bound the space, MAGIC-style scoring to rank
within the residual degeneracy.* **magicmaus is that synthesis, built.**

## The idea

```
HMQC + NOESY peak lists + PDB
        │
        ▼
┌─────────────────────────────┐   layer 1  (MAUS, reused verbatim)
│ SAT subgraph isomorphism     │   per-peak option set O_i
│ → hard constraints           │   truth ∈ O_i  guaranteed
└─────────────────────────────┘
        │  pruned domains (mostly 1–3 candidates, never the truth removed)
        ▼
┌─────────────────────────────┐   layer 2  (MAGIC-style, new)
│ SAT-feasible seed assignment │   one jointly-consistent, injective map
│ + 3-cycle simulated-anneal   │   maximise Σ NOE-contact strength (~1/r⁶)
│ search on NOE-contact score  │   every move stays feasible
└─────────────────────────────┘
        │
        ▼
per peak:  a single committed call · the MAUS option set (ambiguity envelope)
           · a confidence tier (unique / scored / ambiguous)
```

Why it works and why each half is load-bearing:

- **MAUS's option sets are enumerated *independently*** — a methyl is an option
  if *some* satisfying global map uses it. Their product is **not** jointly
  realizable, and they do not say which single coherent assignment is best.
  MAGIC answers exactly that. MAUS alone can't.
- **MAGIC over the full candidate space is near-flat** and commits to
  near-optimal-but-wrong answers (85 % error on this data). Run over MAUS's tiny
  truth-containing domains instead, the same scoring becomes both fast and
  actually able to find the truth. MAGIC alone can't.
- The MAUS option set is still reported as an **honest ambiguity envelope**, so
  the never-exclude guarantee survives even though a single answer is committed.

The scored step is a **feasibility-preserving 3-cycle simulated-annealing search**
seeded by one SAT model of the pruned domains. A plain greedy ascent stalls in a
local optimum 10–20% below the truth's objective on the near-flat MAGIC landscape;
annealing over relocate/swap/3-cycle-rotation moves climbs into the truth's basin
(the rotations cross the tightly coupled Leu/Val option graphs a swap cannot),
then a final greedy ascent locks it in. Starting feasible and keeping every move
injective and NOE-consistent guarantees a valid bijection as output — where a
naïve per-cluster brute force would choke on MBP's single **138-peak** residual
degeneracy cluster and break injectivity. The annealer runs only where its
objective is trustworthy — a NOESY with real intensities *and* ≥75 % firm-NOE
coverage over the optimised peaks; on a boolean network or a sparsely-constrained
target (MSG) it reduces to the plain ascent.

## Install

```bash
python3 -m pip install python-sat numpy      # or: python3 -m venv .venv && .venv/bin/pip install python-sat numpy
```

The SAT layer reuses the MAUS clean-room code: `maus.py` here is a **vendored
copy** of `../maus/maus.py` (kept in sync), so this repo runs standalone. If a
sibling `../maus/` checkout is present it is used instead — the two are
identical. The three-way MAGIC comparison additionally needs a `../magic/`
checkout (plus `psutil`); the core magicmaus pipeline does not.

## Usage

```bash
python magicmaus.py PDB HMQC.tsv NOESY.tsv \
    [--hmbc HMBC.tsv] [--truth TRUTH.tsv] \
    [--tol-h 0.01] [--tol-c 0.05] [--soft-ambiguous] [--out calls.tsv]
```

Input formats are exactly MAUS's (see `../maus/README.md`):

- `HMQC.tsv` — `label ⇥ H_ppm ⇥ C_ppm ⇥ res_type`. `res_type` may be a bare type
  (`L`) or a **tentative anchor** (`L45D2`) that pins one peak and propagates.
- `NOESY.tsv` — `label ⇥ C1 ⇥ C2 ⇥ H2 ⇥ [intensity]` (3D (H)CCH; observed methyl
  `(H2,C2)`, partner carbon `C1`). The **5th `intensity` column is optional** and
  defaults to 1.0, so a plain 4-column MAUS list reads unchanged. When present,
  intensity feeds the score as `intensity × (1/r⁶)` — a strong NOE is pushed onto
  a close structural contact, which is what lets **real intensities break the
  geometric degeneracy** a boolean network cannot (see below).
- `--hmbc` — optional geminal links; `--truth` — key for scoring.
- `--soft-ambiguous` — fold the NOE cross peaks MAUS **discards** (ambiguous
  matches) back in as MAGIC-style diluted soft tie-breakers. Off by default so
  the headline number rests only on hard-consistent data; **on real NOESY with
  intensities it is the intended lever** (see below).

Output (`--out`) is one row per peak: the committed `call`, its `confidence`
tier and `margin`, the full MAUS `options` set, and (with `--truth`) whether the
call and the envelope contain the truth.

## Example — maltose-binding protein

Real BMRB 7114 shifts + PDB 1ANF, structure-simulated NOESY (the shared dataset
of `../maus/COMPARISON.md`):

```bash
python magicmaus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --out mbp_calls.tsv
```

```
--- MAUS envelope (never excludes truth) ---
unique(1 option)      = 44/192
ambiguous(2-3 options)= 73/192
ambiguous(>3 options) = 75/192
--- magicmaus commitment (single coherent call) ---
confidence: unique=44  scored=66  ambiguous=82
--- scored vs truth ---
truth in MAUS option set = 192/192 = 100.0%  (never-exclude guarantee preserved)
magicmaus single call    = 119/192 = 62.0% correct
    unique   : 44/44 = 100.0% correct
    scored   : 42/66 = 63.6% correct
    ambiguous: 33/82 = 40.2% correct
```

### Head-to-head (same MBP inputs)

| | single-answer correct | truth-in-envelope | note |
|---|---|---|---|
| **MAGIC** alone | 15.1 % | — (no envelope) | commits to all, 84.9 % error, ~279 s |
| **MAUS** alone  | 44/192 decisive, rest abstain | **100 %** | never wrong, but no call on degeneracy |
| **magicmaus**   | **62.0 %** | **100 %** | commits on all 192 *and* keeps the envelope, ~0.3 s |

magicmaus turns MAUS's 148 abstentions into 75 additional committed correct
calls **without giving up the 100 % never-exclude guarantee**, and it multiplies
MAGIC's single-answer accuracy 4× (15 % → 62 %) — because the scored search runs
only over the truth-containing MAUS domains. (This is the boolean network; real
NOESY intensities lift it to 87 %, below.)

### Levers

| variant | single-answer correct | why |
|---|---|---|
| default (`--tol-c 0.05`) | 62.0 % | firm, hard-consistent boolean NOE only |
| `--soft-ambiguous`       | 67.2 % | reuses the ambiguous NOEs MAUS throws away |
| tighter `--tol-c 0.02`   | 79.2 % | sharper NOE matching → more unique peaks |
| +24 tentative anchors    | 71.4 % | anchors propagate through **both** layers |

On the **boolean** simulated NOESY the residual coin-flip tier is dominated by
**geminal methyl pairs** (Leu δ1/δ2, Val γ1/γ2) and shift-degenerate peaks —
symmetries an achiral, intensity-free network cannot resolve. magicmaus flags
them `ambiguous` (margin 0) rather than guessing, and reports both members in the
envelope.

### The intensity lever (where MAGIC's scoring earns its keep)

Real NOESY cross peaks carry an intensity ~ 1/r⁶ that says *how close* the
contact is — exactly the signal the boolean network lacks. `make_intensity_noesy.py`
reconstructs that column for the MBP list (from the truth key + structure),
simulating a real experiment:

```bash
python make_intensity_noesy.py examples/mbp/1ANF.pdb examples/mbp/hmqc_true.tsv \
    examples/mbp/noesy.tsv mbp_noesy_intensity.tsv
python magicmaus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv mbp_noesy_intensity.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous
```

| NOESY | single-answer correct | scored-tier correct |
|---|---|---|
| boolean (no intensity)          | 62.0 % | 63.6 % |
| + real intensities              | 87.0 % | 92.3 % |
| + intensities `--soft-ambiguous`| **87.5 %** | **92.6 %** |

Feeding intensities lifts accuracy 62 % → 87 % **with the envelope still at
100 %** — the scoring layer, idle on boolean data, does real work once the data
carry information, precisely as `../maus/COMPARISON.md` predicted. magicmaus
detects the flat boolean network and holds back the annealer there (it would only
overfit to structural-contact density), then unleashes it once intensities pin
each NOE to its true distance.

### Fair three-way, *all three engines on the same intensity NOESY*

`convert_to_magic.py` re-emits this exact intensity dataset as a MAGIC control
bundle, so MAGIC, MAUS and magicmaus are scored on **identical** shifts,
structure and NOE peaks (`score_three.py` reproduces the table):

| engine | methyl-level | residue-level | truth-in-envelope |
|---|---|---|---|
| **MAGIC** (scoring, uses intensity) | 5.7 % | 10.4 % | — (no envelope) |
| **MAUS** (SAT, *ignores* intensity) | 26.6 % unique, rest abstain | — | **100 %** |
| **magicmaus**                       | 87.0 % | 90.1 % | **100 %** |
| **magicmaus** `--soft-ambiguous`    | **87.5 %** | **89.6 %** | **100 %** |

Two facts the shared network makes plain:

- **MAUS structurally cannot use intensity** — its constraints are boolean, so
  the 5th column changes nothing (verified: identical 51 unique / 100 % envelope
  with or without it). It bounds the truth but stays undecided on degeneracy.
- **MAGIC uses intensity but still lands at ~6–10 %** — global scoring over the
  full type-matched space sits on a near-flat landscape (MAGIC's own
  `VALIDATION.md §2` reports the same 4–10 % on simulated NOESY), and it commits
  to a residue-level answer that cannot even resolve geminal pairs.

magicmaus runs the *same* intensity-weighted scoring MAGIC uses, but only inside
MAUS's truth-containing domains — and lands **15× MAGIC's methyl accuracy while
keeping MAUS's 100 % envelope**. That is the synthesis paying off: MAUS keeps the
truth in reach, MAGIC's scoring then extracts every bit of experimental signal to
commit correctly within it. Neither half delivers this alone.

## Benchmark — seven simulated targets + one real multimer

MBP is not a one-off. `make_peaklists.py` builds a dataset from any PDB + BMRB
deposition (see [`examples/*/README.md`](examples/)); `bench.py` scores MAUS and
magicmaus, `convert_to_magic.py` adds MAGIC. Beyond MBP we ship ubiquitin (the
reference protein of biomolecular NMR, where magicmaus assigns all 43 methyls
correctly) and the four de-novo blind targets of the MAUS paper (Nerli *et al.*
2021, Table 1): interleukin-2 and the HNH, REC2 and REC3 domains of Cas9. The
suite also includes **malate synthase G** (257 methyls) — the classic large-protein
methyl benchmark — built from the open-access MethylFLYA supplement (no BMRB
methyl deposit) via `make_peaklists.py --shifts-tsv`.

| target | BMRB / PDB | labeling | methyls | MAGIC | MAUS | magicmaus | +soft | +HMBC | envelope |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| Ubiquitin | 6457 / 1UBQ | AILMTV | 43 | 9.3 % | 34.9 % | **100 %** | 100.0 % | 100.0 % | **100 %** |
| IL-2 | 28104 / 1M47 | ILV | 59 | 8.5 % | 8.5 % | 88.1 % | 96.6 % | 96.6 % | **100 %** |
| HNH  | 27949 / 6O56 | AILTV | 57 | 12.3 % | 26.3 % | 82.5 % | 86.0 % | 86.0 % | **100 %** |
| REC2 | 28105 / 4CMP | ILV | 63 | n.c. | 12.7 % | 88.9 % | 90.5 % | 76.2 % | **100 %** |
| REC3 | 28110 / 4ZT0 | ILV | 85 | n.c. | 8.2 % | 60.0 % | 57.6 % | 52.9 % | **100 %** |
| MBP  | 7114 / 1ANF  | AILMTV | 192 | 5.7 % | 26.6 % | 87.0 % | 87.5 % | 93.2 % | **100 %** |
| MSG  | SI / 1D8C    | ILV | 257 | n.c. | 1.6 % | 29.6 % | 33.5 % | 38.5 % | **100 %** |
| **TNF-α** † | *real expt* / AF3 trimer | AILTV | 85 | — | 7.1 % | 30.6 % | 28.2 % | n.r. | 92.9 % |

† **TNF-α is the one real-experimental, multimeric target** — genuine methyl-NMR
peak lists of the tumour-necrosis-factor-α homotrimer against an AlphaFold3
trimer model, not a structure-simulated NOESY. It is where the guarantees meet
reality (see [`examples/TNFa/`](examples/TNFa/)). Two honest lessons the
simulated targets cannot show: (1) the **100 % envelope is conditional** on the
NOEs being consistent with the structure — measured against a *predicted*
structure, 6/85 peaks carry contacts the model does not support, so the envelope
is 92.9 %, not 100 %; (2) **multimer handling is load-bearing** — parsed as a
single protomer the real inter-subunit NOEs are unexplainable and the SAT
collapses to a 0 % envelope; keeping all three chains as symmetry images (contact
= min distance over subunits) recovers 92.9 %. HMBC is `n.r.` (not recommended)
here: real HMBC shifts resolve poorly, one wrong matched geminal link makes the
global SAT infeasible, and the envelope collapses — the extreme of the same
target-dependent fragility seen on REC2/REC3.

MAUS commits only on the unique peaks (the % shown, all correct) and abstains on
the rest — its coverage is the envelope column. The **100 % envelope holds on
every simulated target** (and 92.9 % on the real TNF-α data, where a predicted
structure cannot support 6 measured NOEs), and magicmaus beats full-space MAGIC by up to ~15× (MBP 87 % vs
5.7 %). The 3-cycle annealer is what closes the gap: on the intensity network the
objective's global optimum *is* the truth (a truth-seeded search scores ~96 %),
and the annealer reaches it where a plain greedy ascent stalls ~10–20 % short.
`+soft-ambiguous` helps on most targets (IL-2 88→97 %, HNH 82→86 %) but is a wash
on the Leu-densest ones (REC3), so it is opt-in. `+HMBC` (an optional geminal-link
experiment, `--hmbc`, on top of +soft) helps where geminal pairs dominate the
residual (MBP 87→93 %) but *hurts* where the degeneracy is symmetric rather than
geminal (REC2, REC3): the extra hard links reshape the landscape into a different
equal-scoring optimum, so it too is opt-in. MSG is the one target where the
objective itself is underdetermined — only 62 % of its optimised peaks carry a firm
NOE (95 of 257 have none), below the 75 % coverage cut — so magicmaus withholds the
annealer and falls back to the safe greedy ascent (29.6 %) rather than chasing an
unreliable optimum. MAGIC did not converge (`n.c.`) within a 15-min budget on the
two Leu-dense Cas9 domains — the scaling cost of scoring the full space. REC3
(50/85 Leu) is the hard case: an achiral NOE network cannot resolve that much
geminal/shift degeneracy, and magicmaus reports the residual as `ambiguous` option
sets rather than guessing.

Regenerate any target (e.g. IL-2):

```bash
python make_peaklists.py examples/il2/1M47.pdb BMRB28104.str --out-dir examples/il2
python make_intensity_noesy.py examples/il2/1M47.pdb examples/il2/hmqc_true.tsv \
    examples/il2/noesy.tsv examples/il2/noesy_intensity.tsv
python bench.py examples/il2 examples/il2/1M47.pdb
```

## Files

```
magicmaus.py             the hybrid engine (MAUS SAT bound + MAGIC scored commit)
maus.py                  vendored MAUS clean-room SAT layer (== ../maus/maus.py)
make_peaklists.py        build hmqc/noesy peak lists from a PDB + BMRB deposition
make_intensity_noesy.py  reconstruct 1/r^6 intensities for a simulated NOESY
make_tnfa_input.py       convert the real TNF-alpha Sparky peak lists -> magicmaus input + truth
convert_to_magic.py      re-emit a dataset as a MAGIC control bundle (fair compare)
bench.py                 MAUS + magicmaus (+MAGIC) on one example dir
score_three.py           MAGIC vs MAUS vs magicmaus on the MBP intensity NOESY
examples/mbp/            MBP dataset (BMRB 7114 shifts + PDB 1ANF)
examples/ubq/            ubiquitin (BMRB 6457 + PDB 1UBQ), the NMR reference protein
examples/{il2,hnh,rec2,rec3}/   MAUS-paper blind targets (BMRB + PDB)
examples/msg/            malate synthase G (MethylFLYA SI shifts + PDB 1D8C), 257 methyls
examples/TNFa/           real-experimental TNF-alpha homotrimer (Sparky lists + AlphaFold3 model)
```

Each `examples/<target>/` is self-contained: the PDB, the BMRB deposition
(`bmr<id>.str`), the generated peak lists (`hmqc*.tsv`, `noesy*.tsv`, `hmbc.tsv`),
the committed magicmaus output (`magicmaus_calls.tsv`), the MAGIC output
(`magic_assignments.tsv`, where MAGIC converged), and a `README.md` with the exact
regenerate/benchmark commands.

```
magicmaus.py             the hybrid engine (MAUS SAT bound + MAGIC scored commit)
maus.py                  vendored MAUS clean-room SAT layer (== ../maus/maus.py)
make_intensity_noesy.py  reconstruct 1/r^6 intensities for a simulated NOESY
convert_to_magic.py      re-emit a dataset as a MAGIC control bundle (fair compare)
score_three.py           MAGIC vs MAUS vs magicmaus on one intensity NOESY
examples/mbp/            MBP dataset (BMRB 7114 shifts + PDB 1ANF)
```
