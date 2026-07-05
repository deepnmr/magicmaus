# Cas9 REC2 domain — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [28105](https://bmrb.io/data_library/summary/index.php?bmrbId=28105) |
| PDB structure | [4CMP](https://www.rcsb.org/structure/4CMP) (15.6 kDa) |

**63 observed methyls** (9 Ile, 48 Leu, 6 Val). The HMQC peaks are anonymised (`P1…P63`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr28105.str`), the generated peak lists, and the MAGIC assignment
output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/rec2/4CMP.pdb examples/rec2/bmr28105.str --out-dir examples/rec2
python make_intensity_noesy.py examples/rec2/4CMP.pdb examples/rec2/hmqc_true.tsv \
    examples/rec2/noesy.tsv examples/rec2/noesy_intensity.tsv
```

Run magicmaus, or the full benchmark (`bench.py` adds MAGIC via the committed
assignments):

```bash
python magicmaus.py examples/rec2/4CMP.pdb examples/rec2/hmqc.tsv \
    examples/rec2/noesy_intensity.tsv --hmbc examples/rec2/hmbc.tsv \
    --truth examples/rec2/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous
python bench.py examples/rec2 examples/rec2/4CMP.pdb
```

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | did not converge (>15 min) | — |
| magicmaus | 74.6% | 100.0% |
| magicmaus +soft-ambiguous | 76.2% | 100.0% |

magicmaus commits a single call for all 63 peaks while preserving the MAUS
never-exclude envelope (100.0%). Soft ambiguous evidence helps here.
