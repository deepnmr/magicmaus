# magicmaus

A hybrid methyl-NMR assignment engine that **fuses the two sibling projects**
[`../maus/`](../maus) (MAUS, Nerli et al., *Nat. Commun.* 2021) and
[`../magic/`](../magic) (MAGIC, Monneau et al., *J. Biomol. NMR* 2017) into one
program that is strictly better than either alone.

**New here? Start with the [step-by-step tutorial](TUTORIAL.md)** (also
[`TUTORIAL.pdf`](TUTORIAL.pdf)) — a tiny 8-methyl walkthrough that shows every
output field and the intensity lever in action.

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
│ + coordinate ascent on the   │   maximise Σ NOE-contact strength (~1/r⁶)
│ NOE-contact objective        │   every move stays feasible
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

The scored step is a **feasibility-preserving coordinate ascent** seeded by one
SAT model of the pruned domains. Starting feasible and keeping every move
injective and NOE-consistent guarantees a valid bijection as output — where a
naïve per-cluster brute force would choke on MBP's single **138-peak** residual
degeneracy cluster and break injectivity.

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
unique(1 option)      = 51/192
ambiguous(2-3 options)= 81/192
ambiguous(>3 options) = 60/192
--- magicmaus commitment (single coherent call) ---
confidence: unique=51  scored=67  ambiguous=74
--- scored vs truth ---
truth in MAUS option set = 192/192 = 100.0%  (never-exclude guarantee preserved)
magicmaus single call    = 121/192 = 63.0% correct
    unique   : 51/51 = 100.0% correct
    scored   : 44/67 = 65.7% correct
    ambiguous: 26/74 = 35.1% correct
```

### Head-to-head (same MBP inputs)

| | single-answer correct | truth-in-envelope | note |
|---|---|---|---|
| **MAGIC** alone | 15.1 % | — (no envelope) | commits to all, 84.9 % error, ~279 s |
| **MAUS** alone  | 51/192 decisive, rest abstain | **100 %** | never wrong, but no call on degeneracy |
| **magicmaus**   | **63.0 %** | **100 %** | commits on all 192 *and* keeps the envelope, ~0.6 s |

magicmaus turns MAUS's 141 abstentions into 70 additional committed correct
calls **without giving up the 100 % never-exclude guarantee**, and it multiplies
MAGIC's single-answer accuracy 4× (15 % → 63 %) at ~400× the speed — because the
scored search runs only over the truth-containing MAUS domains.

### Levers

| variant | single-answer correct | why |
|---|---|---|
| default (`--tol-c 0.05`) | 63.0 % | firm, hard-consistent boolean NOE only |
| `--soft-ambiguous`       | 67.2 % | reuses the ambiguous NOEs MAUS throws away |
| tighter `--tol-c 0.02`   | 81.2 % | sharper NOE matching → more unique peaks |
| +24 tentative anchors    | 70.8 % | anchors propagate through **both** layers |

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
| boolean (no intensity)          | 63.0 % | 65.7 % |
| + real intensities              | 72.9 % | 79.7 % |
| + intensities `--soft-ambiguous`| **79.7 %** | **85.8 %** |

Feeding intensities lifts accuracy 63 % → 80 % **with the envelope still at
100 %** — the scoring layer, idle on boolean data, does real work once the data
carry information, precisely as `../maus/COMPARISON.md` predicted.

### Fair three-way, *all three engines on the same intensity NOESY*

`convert_to_magic.py` re-emits this exact intensity dataset as a MAGIC control
bundle, so MAGIC, MAUS and magicmaus are scored on **identical** shifts,
structure and NOE peaks (`score_three.py` reproduces the table):

| engine | methyl-level | residue-level | truth-in-envelope |
|---|---|---|---|
| **MAGIC** (scoring, uses intensity) | 5.7 % | 10.4 % | — (no envelope) |
| **MAUS** (SAT, *ignores* intensity) | 26.6 % unique, rest abstain | — | **100 %** |
| **magicmaus**                       | 72.9 % | 79.2 % | **100 %** |
| **magicmaus** `--soft-ambiguous`    | **79.7 %** | **85.4 %** | **100 %** |

Two facts the shared network makes plain:

- **MAUS structurally cannot use intensity** — its constraints are boolean, so
  the 5th column changes nothing (verified: identical 51 unique / 100 % envelope
  with or without it). It bounds the truth but stays undecided on degeneracy.
- **MAGIC uses intensity but still lands at ~6–10 %** — global scoring over the
  full type-matched space sits on a near-flat landscape (MAGIC's own
  `VALIDATION.md §2` reports the same 4–10 % on simulated NOESY), and it commits
  to a residue-level answer that cannot even resolve geminal pairs.

magicmaus runs the *same* intensity-weighted scoring MAGIC uses, but only inside
MAUS's truth-containing domains — and lands **13× MAGIC's methyl accuracy while
keeping MAUS's 100 % envelope**. That is the synthesis paying off: MAUS keeps the
truth in reach, MAGIC's scoring then extracts every bit of experimental signal to
commit correctly within it. Neither half delivers this alone.

## Files

```
magicmaus.py             the hybrid engine (MAUS SAT bound + MAGIC scored commit)
maus.py                  vendored MAUS clean-room SAT layer (== ../maus/maus.py)
make_intensity_noesy.py  reconstruct 1/r^6 intensities for a simulated NOESY
convert_to_magic.py      re-emit a dataset as a MAGIC control bundle (fair compare)
score_three.py           MAGIC vs MAUS vs magicmaus on one intensity NOESY
examples/mbp/            MBP dataset (BMRB 7114 shifts + PDB 1ANF)
```
