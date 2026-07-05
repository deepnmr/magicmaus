"""Guard for the simulated-annealing commitment layer. The search must (1) never
break feasibility -- the committed map stays injective and every firm NOE edge
sits on a structural contact -- and (2) never score below a single greedy ascent,
since it only keeps strictly-better restarts (and falls back to that ascent when
the NOESY carries no intensity signal).

Run: python test_magicmaus_ils.py
"""
from pathlib import Path

import maus
import magicmaus as mm


def _engine(exdir, pdb):
  lab = maus.parse_labeling('A;I;L;M;T;V')
  meth = maus.parse_structure(Path(f'{exdir}/{pdb}').read_text().splitlines(), lab)
  peaks = maus.load_hmqc(f'{exdir}/hmqc.tsv')
  cross = mm.load_noesy_rows(f'{exdir}/noesy_intensity.tsv')
  gem, sg, lg = maus.build_structure_graph(meth, 6.0, 10.0)
  noe, ei, _amb, _ = mm.match_noe_intensity(peaks, cross, 0.01, 0.05)
  return mm.MagicMaus(meth, peaks, gem, sg, lg, noe, edge_intensity=ei)


def check(exdir, pdb):
  eng = _engine(exdir, pdb)
  options = eng.option_sets()
  free = [p.index for p in eng.peaks if len(options[p.index]) > 1]
  single = eng._improve(eng._feasible_model(options), options, free)
  chosen, opts = eng.solve()

  # (1) injective: no methyl assigned to two peaks
  used = [g for g in chosen.values() if g is not None]
  assert len(used) == len(set(used)), 'ILS output not injective'
  # (1) NOE-consistent: every firm edge on an allowed structural contact
  for (i, j) in eng.noe_edges:
    gi, gj = chosen[i], chosen[j]
    if gi is not None and gj is not None and gi != gj:
      assert (gi, gj) in eng.allowed, f'firm NOE {i},{j} off-contact'
  # (1) every call stays inside the MAUS envelope (never-exclude preserved)
  for p in eng.peaks:
    assert chosen[p.index] in opts[p.index], 'call left MAUS option set'
  # (2) ILS never worse than a single ascent
  assert eng._total_obj(chosen) >= eng._total_obj(single) - 1e-12, 'ILS below single ascent'
  print(f'OK: {exdir} obj single={eng._total_obj(single):.4g} ils={eng._total_obj(chosen):.4g}')


if __name__ == '__main__':
  check('examples/dummy', 'model.pdb')
  check('examples/ubq', '1UBQ.pdb')
