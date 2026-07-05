# Interleukin-2 (IL-2) — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [28104](https://bmrb.io/data_library/summary/index.php?bmrbId=28104) |
| PDB structure | [1M47](https://www.rcsb.org/structure/1M47) (15.4 kDa) |

**59 observed methyls** (9 Ile, 42 Leu, 8 Val). The HMQC peaks are anonymised (`P1…P59`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr28104.str`), the generated peak lists, the committed magicmaus
output (`magicmaus_calls.tsv` for the +soft protocol, `magicmaus_calls_hmbc.tsv`
with the optional HMBC lever: per-peak call, confidence tier, option set, truth),
and the MAGIC assignment output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/il2/1M47.pdb examples/il2/bmr28104.str --out-dir examples/il2
python make_intensity_noesy.py examples/il2/1M47.pdb examples/il2/hmqc_true.tsv \
    examples/il2/noesy.tsv examples/il2/noesy_intensity.tsv
```

Run magicmaus (the benchmark protocol: firm NOE + `--soft-ambiguous`, matching
`magicmaus_calls.tsv` and the table below), or the full head-to-head with
`bench.py`:

```bash
python magicmaus.py examples/il2/1M47.pdb examples/il2/hmqc.tsv \
    examples/il2/noesy_intensity.tsv \
    --truth examples/il2/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous \
    --out examples/il2/magicmaus_calls.tsv
python bench.py examples/il2 examples/il2/1M47.pdb \
    --magic examples/il2/magic_assignments.tsv
```

`hmbc.tsv` is an optional geminal-link lever: adding `--hmbc examples/il2/hmbc.tsv`
resolves Leu/Val geminal pairs and raises accuracy further (not part of the table).

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | 8.5% | — |
| MAUS (unique only, rest abstain) | 8.5% | 100.0% |
| magicmaus | 88.1% | 100.0% |
| magicmaus +soft-ambiguous | 89.8% | 100.0% |
| magicmaus +soft +HMBC | 89.8% | 100.0% |

magicmaus commits a single call for all 59 peaks while preserving the MAUS
never-exclude envelope (100.0%). Soft ambiguous evidence helps here. Adding the optional HMBC geminal-link
experiment (`--hmbc`) on top of +soft reaches 89.8%.
