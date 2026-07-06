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
here — many real NOEs cross the subunit interface:

| structure model | firm NOE edges | truth-in-envelope | unique |
|---|---|---|---|
| single protomer (chain 1 only) | fewer | **0/85** (SAT collapses) | 0 |
| full trimer (all 3 chains) | more | **79/85 = 92.9%** | 6 |

Parsed as a monomer the real inter-subunit NOEs are structurally
unexplainable, the hard constraints become jointly infeasible, and the option
sets collapse. The trimer parse recovers a coherent 92.9% envelope.

## Run

```bash
python magicmaus.py examples/TNFa/fold_tnfa_trimer_model_0.cif \
    examples/TNFa/hmqc.tsv examples/TNFa/noesy_intensity.tsv \
    --truth examples/TNFa/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05
```

```
methyls(G nodes)=89  HMQC peaks=85  NOESY cross peaks=220
NOE match (tol H±0.01/C±0.05): firm=65 ambiguous(dropped)=111 unmatched=44
truth in MAUS option set = 79/85 = 92.9%
magicmaus single call    = 26/85 = 30.6% correct   (unique 6/6, scored 14/38, ambiguous 6/41)
```

Or the same three-engine benchmark as the other targets:

```bash
python bench.py examples/TNFa examples/TNFa/fold_tnfa_trimer_model_0.cif
```

## What the numbers mean

Two honest facts real data exposes that the simulated targets cannot:

- **The 100% guarantee is conditional.** MAUS never excludes the truth *when the
  NOEs are consistent with the structure*. On real data measured against a
  predicted structure, 6/85 peaks (A33, L75δ2, V91γ1/γ2, L94δ1/δ2) carry NOEs
  the AlphaFold model does not support at the 6/10 Å cutoffs, so they fall out
  of the envelope. 92.9%, not 100%, is the real-world number.
- **HMBC hard-constraints are fragile to noise.** `--hmbc` turns each matched
  HMBC cross peak into a hard geminal link. Real HMBC shifts resolve poorly at
  these tolerances (65 of 69 rows unmatched); the few that match include a
  *wrong* geminal pair, and one bad hard constraint makes the global SAT
  infeasible → the envelope collapses to 0. On this dataset HMBC is **not
  recommended** — consistent with its target-dependent role elsewhere, taken to
  its extreme by real spectral noise.

Accuracy (30.6%) is lower than the simulated targets because the input is a
boolean-ish 3D `(H)CCH` network with raw peak-height intensities, matched
against a predicted structure — the hardest, most realistic setting in the
benchmark.
