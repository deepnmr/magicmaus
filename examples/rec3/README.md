# Cas9 REC3 domain — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [28110](https://bmrb.io/data_library/summary/index.php?bmrbId=28110) |
| PDB structure | [4ZT0](https://www.rcsb.org/structure/4ZT0) (24.5 kDa) |

**85 observed methyls** (13 Ile, 50 Leu, 22 Val). The HMQC peaks are anonymised (`P1…P85`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr28110.str`), the generated peak lists, and the MAGIC assignment
output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/rec3/4ZT0.pdb examples/rec3/bmr28110.str --out-dir examples/rec3
python make_intensity_noesy.py examples/rec3/4ZT0.pdb examples/rec3/hmqc_true.tsv \
    examples/rec3/noesy.tsv examples/rec3/noesy_intensity.tsv
```

Run magicmaus, or the full benchmark (`bench.py` adds MAGIC via the committed
assignments):

```bash
python magicmaus.py examples/rec3/4ZT0.pdb examples/rec3/hmqc.tsv \
    examples/rec3/noesy_intensity.tsv --hmbc examples/rec3/hmbc.tsv \
    --truth examples/rec3/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous
python bench.py examples/rec3 examples/rec3/4ZT0.pdb
```

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | did not converge (>15 min) | — |
| magicmaus | 32.9% | 100.0% |
| magicmaus +soft-ambiguous | 28.2% | 100.0% |

magicmaus commits a single call for all 85 peaks while preserving the MAUS
never-exclude envelope (100.0%). Soft ambiguous evidence does not help on this target (dense Leu degeneracy), so the plain call is preferred.
