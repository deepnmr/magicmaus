"""Regression guard for the MAUS speed optimisations (arc-consistency domain
prune + sequential-counter at-most-one). They must not change any result: same
option sets, same unique count, and the never-exclude guarantee intact.

Run: python test_maus_fast.py
"""
from pathlib import Path

import maus


def _mbp():
  lab = maus.parse_labeling('A;I;L;M;T;V')
  meth = maus.parse_structure(Path('examples/mbp/1ANF.pdb').read_text().splitlines(), lab)
  peaks = maus.load_hmqc('examples/mbp/hmqc.tsv')
  cross = maus.load_noesy('examples/mbp/noesy_intensity.tsv')
  gem, sg, lg = maus.build_structure_graph(meth, 6.0, 10.0)
  noe, _ = maus.match_noe(peaks, cross, 0.01, 0.05)
  truth = maus.load_truth('examples/mbp/hmqc_true.tsv')
  m = maus.MAUS(meth, peaks, gem, sg, lg, noe)
  opts = m.solve_options()
  lbi = {x.index: x.label for x in meth}
  return peaks, opts, lbi, truth


def test_mbp_unchanged():
  peaks, opts, lbi, truth = _mbp()
  unique = sum(1 for p in peaks if len(opts[p.index]) == 1)
  in_set = sum(1 for p in peaks
               if truth.get(p.peak_id) in [lbi[g] for g in opts[p.index]])
  assert len(peaks) == 192, len(peaks)
  assert unique == 51, unique                       # optimisation must not change this
  assert in_set == 192, in_set                      # never-exclude guarantee holds
  # arc-consistency must never drop the true methyl from a peak's domain
  for p in peaks:
    t = truth.get(p.peak_id)
    if t:
      assert t in [lbi[g] for g in opts[p.index]], (p.peak_id, t)
  print(f'OK: MBP 192 peaks, {unique} unique, {in_set}/192 truth-in-envelope')


if __name__ == '__main__':
  test_mbp_unchanged()
