"""Regenerate per-example README.md from the benchmark table below.
Run from repo root: python examples/gen_readmes.py
"""
from pathlib import Path

# key: (name, bmrb, pdb, mw_kda, comp, MAGIC%, MAUS_unique%, mm%, mm_soft%, mm_hmbc%, env%, n)
DATA = {
  'il2':  ('Interleukin-2 (IL-2)',            28104, '1M47', 15.4, '9 Ile, 42 Leu, 8 Val',
           8.5, 8.5, 88.1, 89.8, 89.8, 100.0, 59),
  'hnh':  ('Cas9 HNH nuclease domain',        27949, '6O56', 15.7, '2 Ala, 7 Ile, 28 Leu, 4 Thr, 16 Val',
           12.3, 26.3, 73.7, 57.9, 64.9, 100.0, 57),
  'rec2': ('Cas9 REC2 domain',                28105, '4CMP', 15.6, '9 Ile, 48 Leu, 6 Val',
           None, 12.7, 74.6, 76.2, 82.5, 100.0, 63),
  'rec3': ('Cas9 REC3 domain',                28110, '4ZT0', 24.5, '13 Ile, 50 Leu, 22 Val',
           None, 8.2, 32.9, 28.2, 45.9, 100.0, 85),
}

TEMPLATE = """# {name} — magicmaus example

Real benchmark target from the MAUS study (Nerli *et al.*, 2021, Table 1), built
from experimental data with `make_peaklists.py`:

| source | id |
|---|---|
| BMRB chemical shifts | [{bmrb}](https://bmrb.io/data_library/summary/index.php?bmrbId={bmrb}) |
| PDB structure | [{pdb}](https://www.rcsb.org/structure/{pdb}) ({mw} kDa) |

**{n} observed methyls** ({comp}). The HMQC peaks are anonymised (`P1…P{n}`); the
truth key (`hmqc_true.tsv`) is used only for scoring. The NOESY is a simulated 3D
`(H)CCH` network (`noesy.tsv` boolean, `noesy_intensity.tsv` with `1/r⁶`
intensities), re-matched to HMQC peaks by frequency at run time.

Everything to reproduce this target lives in this directory: the PDB, the BMRB
deposition (`bmr{bmrb}.str`), the generated peak lists, the committed magicmaus
output (`magicmaus_calls.tsv` for the +soft protocol, `magicmaus_calls_hmbc.tsv`
with the optional HMBC lever: per-peak call, confidence tier, option set, truth),
and the MAGIC assignment output (`magic_assignments.tsv`) where MAGIC converged.

Regenerate the peak lists:

```bash
python make_peaklists.py examples/{k}/{pdb}.pdb examples/{k}/bmr{bmrb}.str --out-dir examples/{k}
python make_intensity_noesy.py examples/{k}/{pdb}.pdb examples/{k}/hmqc_true.tsv \\
    examples/{k}/noesy.tsv examples/{k}/noesy_intensity.tsv
```

Run magicmaus (the benchmark protocol: firm NOE + `--soft-ambiguous`, matching
`magicmaus_calls.tsv` and the table below), or the full head-to-head with
`bench.py`:

```bash
python magicmaus.py examples/{k}/{pdb}.pdb examples/{k}/hmqc.tsv \\
    examples/{k}/noesy_intensity.tsv \\
    --truth examples/{k}/hmqc_true.tsv --tol-h 0.01 --tol-c 0.05 --soft-ambiguous \\
    --out examples/{k}/magicmaus_calls.tsv
python bench.py examples/{k} examples/{k}/{pdb}.pdb{magic_arg}
```

`hmbc.tsv` is an optional geminal-link lever: adding `--hmbc examples/{k}/hmbc.tsv`
resolves Leu/Val geminal pairs and raises accuracy further (not part of the table).

## Benchmark (same intensity NOESY for all engines)

| engine | methyl-level correct | truth-in-envelope |
|---|---|---|
| MAGIC (scoring) | {magic} | — |
| MAUS (unique only, rest abstain) | {maus}% | {env}% |
| magicmaus | {mm}% | {env}% |
| magicmaus +soft-ambiguous | {mms}% | {env}% |
| magicmaus +soft +HMBC | {mmh}% | {env}% |

magicmaus commits a single call for all {n} peaks while preserving the MAUS
never-exclude envelope ({env}%). {soft_note} Adding the optional HMBC geminal-link
experiment (`--hmbc`) on top of +soft reaches {mmh}%.
"""


def main():
  for k, (name, bmrb, pdb, mw, comp, magic, maus, mm, mms, mmh, env, n) in DATA.items():
    magic_s = f'{magic}%' if magic is not None else 'did not converge (>15 min)'
    soft_note = ('Soft ambiguous evidence helps here.' if mms >= mm
                 else 'Soft ambiguous evidence does not help on this target '
                      '(dense Leu degeneracy), so the plain call is preferred.')
    magic_arg = f' \\\n    --magic examples/{k}/magic_assignments.tsv' if magic is not None else ''
    txt = TEMPLATE.format(name=name, bmrb=bmrb, pdb=pdb, mw=mw, comp=comp, n=n,
                          k=k, magic=magic_s, maus=maus, mm=mm, mms=mms, mmh=mmh,
                          env=env, soft_note=soft_note, magic_arg=magic_arg)
    Path(f'examples/{k}/README.md').write_text(txt)
    print(f'wrote examples/{k}/README.md')


if __name__ == '__main__':
  main()
