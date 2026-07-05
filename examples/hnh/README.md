# Cas9 HNH nuclease domain — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [27949](https://bmrb.io/data_library/summary/index.php?bmrbId=27949) |
| PDB structure | [6O56](https://www.rcsb.org/structure/6O56) (15.7 kDa) |

**57 observed methyls** (2 Ala, 7 Ile, 28 Leu, 4 Thr, 16 Val). The HMQC peaks are anonymised (`P1…P57`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr27949.str`), the generated peak lists, the committed magicmaus
output (`magicmaus_calls.tsv`: per-peak call, confidence tier, option set, truth),
and the MAGIC assignment output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/hnh/6O56.pdb examples/hnh/bmr27949.str --out-dir examples/hnh
python make_intensity_noesy.py examples/hnh/6O56.pdb examples/hnh/hmqc_true.tsv \
    examples/hnh/noesy.tsv examples/hnh/noesy_intensity.tsv
```

Run magicmaus (the benchmark protocol: firm NOE + `--soft-ambiguous`, matching
`magicmaus_calls.tsv` and the table below), or the full head-to-head with
`bench.py`:

```bash
python magicmaus.py examples/hnh/6O56.pdb examples/hnh/hmqc.tsv \
    examples/hnh/noesy_intensity.tsv \
    --truth examples/hnh/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous \
    --out examples/hnh/magicmaus_calls.tsv
python bench.py examples/hnh examples/hnh/6O56.pdb \
    --magic examples/hnh/magic_assignments.tsv
```

`hmbc.tsv` is an optional geminal-link lever: adding `--hmbc examples/hnh/hmbc.tsv`
resolves Leu/Val geminal pairs and raises accuracy further (not part of the table).

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | 12.3% | — |
| magicmaus | 73.7% | 100.0% |
| magicmaus +soft-ambiguous | 57.9% | 100.0% |

magicmaus commits a single call for all 57 peaks while preserving the MAUS
never-exclude envelope (100.0%). Soft ambiguous evidence does not help on this target (dense Leu degeneracy), so the plain call is preferred.
