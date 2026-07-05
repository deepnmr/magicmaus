# Ubiquitin — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [6457](https://bmrb.io/data_library/summary/index.php?bmrbId=6457) |
| PDB structure | [1UBQ](https://www.rcsb.org/structure/1UBQ) (8.6 kDa) |

**43 observed methyls** (2 Ala, 7 Ile, 18 Leu, 1 Met, 7 Thr, 8 Val). The HMQC peaks are anonymised (`P1…P43`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr6457.str`), the generated peak lists, the committed magicmaus
output (`magicmaus_calls.tsv` for the +soft protocol, `magicmaus_calls_hmbc.tsv`
with the optional HMBC lever: per-peak call, confidence tier, option set, truth),
and the MAGIC assignment output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/ubq/1UBQ.pdb examples/ubq/bmr6457.str --out-dir examples/ubq
python make_intensity_noesy.py examples/ubq/1UBQ.pdb examples/ubq/hmqc_true.tsv \
    examples/ubq/noesy.tsv examples/ubq/noesy_intensity.tsv
```

Run magicmaus (the benchmark protocol: firm NOE + `--soft-ambiguous`, matching
`magicmaus_calls.tsv` and the table below), or the full head-to-head with
`bench.py`:

```bash
python magicmaus.py examples/ubq/1UBQ.pdb examples/ubq/hmqc.tsv \
    examples/ubq/noesy_intensity.tsv \
    --truth examples/ubq/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous \
    --out examples/ubq/magicmaus_calls.tsv
python bench.py examples/ubq examples/ubq/1UBQ.pdb \
    --magic examples/ubq/magic_assignments.tsv
```

`hmbc.tsv` is an optional geminal-link lever: adding `--hmbc examples/ubq/hmbc.tsv`
resolves Leu/Val geminal pairs and raises accuracy further (not part of the table).

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | 9.3% | — |
| MAUS (unique only, rest abstain) | 34.9% | 100.0% |
| magicmaus | 100.0% | 100.0% |
| magicmaus +soft-ambiguous | 90.7% | 100.0% |
| magicmaus +soft +HMBC | 90.7% | 100.0% |

magicmaus commits a single call for all 43 peaks while preserving the MAUS
never-exclude envelope (100.0%). Soft ambiguous evidence does not help on this target, so the plain call is preferred. Adding the optional HMBC geminal-link
experiment (`--hmbc`) on top of +soft reaches 90.7%.
