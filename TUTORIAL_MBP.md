# magicmaus tutorial — maltose-binding protein (real data)

The [dummy walkthrough](TUTORIAL.md) showed the mechanics on 8 methyls. This one
runs the same pipeline on a **real, full-size dataset**: the 370-residue
*E. coli* maltose-binding protein (MBP), **192 methyls**, using experimental
chemical shifts from **BMRB entry 7114** and the crystal structure **PDB 1ANF**.
By the end you will have reproduced the headline result — magicmaus committing a
single answer for all 192 peaks while never dropping the truth from its envelope
— and the three-way comparison against MAGIC and MAUS on identical data.

Prerequisites: `.venv` with `python-sat` and `numpy` (see [`README.md`](README.md)).
All commands run from the repo root.

---

## 1. The dataset

`examples/mbp/` ships the real inputs:

| file | what it is |
|---|---|
| `1ANF.pdb` | MBP crystal structure |
| `hmqc.tsv` | 192 methyl HMQC peaks (anonymous ids `P1…P192`) — **the input** |
| `hmqc_true.tsv` | answer key (BMRB assignment of each peak) |
| `noesy.tsv` | simulated 3D (H)CCH NOE network (nearest-8 within 8 Å), boolean |
| `hmqc_tentative.tsv` | same HMQC but with 24 peaks pre-anchored (used in §4) |

192 methyls across 6 labelled types: Ala 44, Ile 22, Leu 60, Met 6, Thr 20,
Val 40. That is a realistic large-protein methyl assignment problem.

---

## 2. First run — boolean NOESY

```bash
.venv/bin/python magicmaus.py \
    examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 \
    --out mbp_calls.tsv
```

```
methyls(G nodes)=192  HMQC peaks=192  NOESY cross peaks=1650
NOE match (tol H+-0.01/C+-0.05): firm=502 ambiguous(dropped)=1148 unmatched=0
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

Read it top to bottom:

- The SAT layer resolves **51/192 peaks uniquely** (all correct) and leaves 141
  as multi-option sets — but **the truth is in the option set for all 192**
  (`truth in MAUS option set = 100 %`). That envelope guarantee never breaks.
- On top of that, magicmaus **commits a single answer for all 192** peaks:
  **63.0 % correct**, tiered by confidence (51 `unique`, 67 `scored`, 74
  `ambiguous`). That is the number to beat with better data (§3–4).

### The per-peak table

```bash
head -9 mbp_calls.tsv | column -t -s $'\t'
```

```
label  res_type  call    confidence  n_options  options        truth   call_correct  truth_in_set
P1     L         L7CD2   scored      2          L7CD1,L7CD2     L7CD1   0             1
P2     L         L7CD1   scored      2          L7CD1,L7CD2     L7CD2   0             1
P3     V         V8CG1   ambiguous   2          V8CG1,V8CG2     V8CG1   1             1
P4     V         V8CG2   scored      2          V8CG1,V8CG2     V8CG2   1             1
P5     I         I9CD1   unique      1          I9CD1           I9CD1   1             1
P6     I         I11CD1  unique      1          I11CD1          I11CD1  1             1
P7     L         L20CD1  scored      2          L20CD1,L20CD2   L20CD1  1             1
P8     L         L20CD2  ambiguous   2          L20CD1,L20CD2   L20CD2  1             1
```

`P5`/`P6` are pinned uniquely (`unique`, certain). `P1`/`P2` are a Leu geminal
pair the boolean network cannot orient, so magicmaus swaps them here — but note
`truth_in_set = 1` for both: the right answer is still inside each envelope.
This is exactly the boolean coin-flip the intensity column fixes next.

---

## 3. Add intensities

Real NOE cross peaks carry an intensity ~ 1/r⁶ — the "how close" signal a boolean
edge lacks. `make_intensity_noesy.py` reconstructs it for this dataset (from the
truth key + structure), simulating a real experiment:

```bash
.venv/bin/python make_intensity_noesy.py examples/mbp/1ANF.pdb \
    examples/mbp/hmqc_true.tsv examples/mbp/noesy.tsv mbp_noesy_intensity.tsv

.venv/bin/python magicmaus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv \
    mbp_noesy_intensity.tsv --truth examples/mbp/hmqc_true.tsv \
    --tol-h 0.01 --tol-c 0.05 --soft-ambiguous
```

| NOESY | single-answer correct | scored-tier correct |
|---|---|---|
| boolean                          | 63.0 % | 65.7 % |
| + intensities                    | 87.0 % | 92.3 % |
| + intensities `--soft-ambiguous` | **87.5 %** | **92.6 %** |

`--soft-ambiguous` folds in the ambiguous NOE cross peaks MAUS discards as extra,
intensity-weighted tie-breakers. Intensities lift the run **63 % → 87 %** — and
the envelope stays at **100 %** the whole way. The scoring layer, idle on boolean
data, does real work once the data carry information: with the intensity column in
hand the NOE objective's optimum is the truth, and magicmaus's 3-cycle annealer
climbs to it (a plain greedy ascent stalls ~15 % short).

---

## 4. Two more levers

**Tighter matching.** Sharper shift tolerances make more NOEs resolve firmly, so
more peaks land in the `unique` tier:

```bash
.venv/bin/python magicmaus.py examples/mbp/1ANF.pdb examples/mbp/hmqc.tsv \
    examples/mbp/noesy.tsv --truth examples/mbp/hmqc_true.tsv \
    --tol-h 0.01 --tol-c 0.02
# unique(1 option) = 70/192 ;  single call = 156/192 = 81.2% correct
```

**Tentative anchors.** If you already know a few assignments, write them into the
`res_type` cell of the HMQC file (e.g. `L45D2` instead of `L`). `hmqc_tentative.tsv`
carries 24 such anchors; they propagate through the NOE network and sharpen the
rest:

```bash
.venv/bin/python magicmaus.py examples/mbp/1ANF.pdb \
    examples/mbp/hmqc_tentative.tsv examples/mbp/noesy.tsv \
    --truth examples/mbp/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05
# unique(1 option) = 79/192 ;  single call = 136/192 = 70.8% correct
```

24 anchors lift the `unique` tier from 51 to 79 — a few knowns go a long way.

---

## 5. The payoff — MAGIC vs MAUS vs magicmaus on identical data

`convert_to_magic.py` re-emits the intensity dataset as a MAGIC control bundle,
so all three engines are scored on the **same** shifts, structure and NOE peaks.
`score_three.py` runs the comparison:

```bash
.venv/bin/python convert_to_magic.py examples/mbp/hmqc.tsv \
    mbp_noesy_intensity.tsv examples/mbp/1ANF.pdb magic_run_intensity
# then run MAGIC on magic_run_intensity/control.txt (needs a ../magic checkout),
# and:
.venv/bin/python score_three.py
```

```
                            methyl-level     residue-level   truth-in-envelope
MAGIC (scoring)           11/192 =  5.7%    20/192 = 10.4%     — (no envelope)
MAUS (SAT)            unique 51/192 = 26.6%                      192/192 = 100.0%
magicmaus                167/192 = 87.0%   173/192 = 90.1%    192/192 = 100.0%
magicmaus +soft-amb      168/192 = 87.5%   172/192 = 89.6%    192/192 = 100.0%
```

Two facts the shared network makes plain:

- **MAUS cannot use intensity** — its constraints are boolean, so the 5th column
  changes nothing (identical 51 unique / 100 % envelope). It bounds the truth but
  stays undecided on the degeneracy.
- **MAGIC uses intensity but still lands at ~6–10 %** — global scoring over the
  full candidate space sits on a near-flat landscape and commits to a
  residue-level answer that cannot even resolve geminal pairs.

magicmaus runs the *same* intensity-weighted scoring MAGIC uses, but only inside
MAUS's truth-containing domains — landing **15× MAGIC's methyl accuracy while
keeping MAUS's 100 % envelope**. That is the synthesis: MAUS keeps the truth in
reach, MAGIC's scoring extracts every bit of experimental signal to commit
correctly within it. Neither half delivers this alone.

---

## 6. What to keep from a real run

On data without a truth key, drop `--truth` and read `mbp_calls.tsv`:

- **`unique`** calls — forced by hard constraints; trust them.
- **`scored`** calls — the NOE geometry's best pick; strong leads, verify where
  it matters.
- **`ambiguous`** calls — carry the *whole option set*, not the single call;
  these are genuine symmetries (mostly geminal pairs) awaiting stereo labelling,
  a 4D experiment, or more anchors.

The one invariant across every setting above: **the true methyl is never dropped
from a peak's option set.** magicmaus commits when it can and abstains — visibly —
when it cannot.

## 7. Beyond MBP — the five-protein benchmark

MBP is one of five targets. `make_peaklists.py` builds the same peak lists from
any PDB + BMRB deposition, so the four de-novo blind targets of the MAUS paper
(IL-2 and the Cas9 HNH/REC2/REC3 domains) drop straight in — see
[`../README.md#benchmark--five-real-bmrb-targets`](../README.md) and each
`examples/<target>/README.md`. Across all five the **100 % never-exclude envelope
holds**, and magicmaus beats full-space MAGIC by up to ~15× while committing on
every peak.
