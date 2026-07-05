# Malate synthase G (MSG) — magicmaus example

The landmark large-protein methyl-assignment target: the 723-residue, 82 kDa
enzyme MSG, the first protein of this size whose Ile/Leu/Val methyls were
assigned by solution NMR (Tugarinov, Ollerenshaw & Kay, landmark ILV methyl
assignment). This is the biggest, most degenerate target in the suite.

## Data provenance

Unlike the other examples, MSG's methyl shifts are **not** in a BMRB deposition
(BMRB 5471 is backbone-only). They are digitised from the reference methyl
assignment table published, open access (CC-BY), as Supplementary Table 1 of:

> Pritišanac, Würz, Alderson & Güntert, *Automatic structure-based NMR methyl
> resonance assignment in large proteins*, Nat. Commun. 10, 4922 (2019).
> <https://doi.org/10.1038/s41467-019-12837-8>

| file | contents |
|---|---|
| `msg_methyl_shifts.tsv` | paired ¹H/¹³C methyl shifts (`resnum res_type atom H_ppm C_ppm`), 268 ILV methyls parsed from the SI |
| `1D8C.pdb` | MSG crystal structure (chain-A ATOM records) |
| `hmqc.tsv`, `hmqc_true.tsv`, `noesy.tsv`, `noesy_intensity.tsv`, `hmbc.tsv` | generated peak lists |

The parsed shifts align to `1D8C` numbering with **zero residue-type mismatches**;
257 of 268 methyls map onto a structural methyl (the rest are disordered/missing
in the crystal). MSG is ILV-labelled: 41 Ile, 133 Leu, 83 Val.

Because MSG's methyl shifts are a paired ¹H/¹³C table (not a BMRB `.str`),
`make_peaklists.py` reads them with `--shifts-tsv`:

```bash
python make_peaklists.py examples/msg/1D8C.pdb \
    --shifts-tsv examples/msg/msg_methyl_shifts.tsv \
    --out-dir examples/msg --noe-cut 7.9 --labeling "I;L;V"
python make_intensity_noesy.py examples/msg/1D8C.pdb examples/msg/hmqc_true.tsv \
    examples/msg/noesy.tsv examples/msg/noesy_intensity.tsv
python msg_run.py
```

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | did not converge (>15 min) | — |
| MAUS (unique only, rest abstain) | 1.6% | 100.0% |
| magicmaus | 29.6% | 100.0% |
| magicmaus +soft-ambiguous | 33.5% | 100.0% |
| magicmaus +soft +HMBC | 38.5% | 100.0% |

MSG is the hard case of the suite: only 262 of its 3D (H)CCH cross peaks resolve
to a firm NOE (carbon-only partner matching is highly degenerate), so 95 of 257
peaks carry no firm constraint and the true option sets are large. The
never-exclude envelope still holds at 100%, and magicmaus commits on every peak
where full-space MAGIC does not even converge. This is also the one target where
the scoring objective is underdetermined: with only 62% of the free peaks pinned by
a firm NOE (below the 75% cut), the objective's global optimum is no longer the
truth, so the 3-cycle annealer would climb to a higher-scoring but *less* accurate
map. magicmaus detects the low coverage and falls back to the safe greedy ascent
here — the numbers above are that fallback. The enumeration is compute-heavy at this
scale (weak constraints + 133-Leu symmetry); `msg_run.py` reuses a single
option-set enumeration for the plain and +soft calls to keep it practical.
